import logging
import math
import os
import pickle
import re
import sys
from abc import abstractmethod
from collections import OrderedDict, defaultdict
from enum import Enum
from itertools import starmap
from operator import methodcaller
from typing import List, Union

import casadi as ca
from casadi import (
    Function,
    MX,
    SX,
    det,
    hessian,
    jacobian,
    substitute,
    symvar,
    trace,
    vertcat,
)

import numpy as np
from numpy import inf

from rtctools.optimization.goal_programming_mixin import Goal
from rtctools.optimization.optimization_problem import OptimizationProblem
from rtctools.optimization.timeseries import Timeseries

from .polygon_enclosure import DeadEndError, enclosing_segments
from .util import _ObjectParameterWrapper

logger = logging.getLogger(__name__)


class CommonStructureSwitchFunctions:
    def StructureHistory(
        self, structure, times, status_sym, pump_history=None, station_history=None
    ):
        hist_needed = max(structure.minimum_off, structure.minimum_on)
        hist_earliest = np.min(times - hist_needed)

        history_map = {Pump: pump_history, PumpingStation: station_history}

        if type(structure) in history_map:
            history = history_map[type(structure)]
            hist_start_ind = np.searchsorted(history.times, hist_earliest)
            if history.times[hist_start_ind] != hist_earliest:
                # We need an earlier time step to be sure
                hist_start_ind -= 1
        else:
            raise Exception("Unknown structure type for history processing.")

        hist_end_ind = np.searchsorted(history.times, times[0])
        hist_times = history.times[hist_start_ind : hist_end_ind + 1]
        hist_status = history.values[hist_start_ind : hist_end_ind + 1]

        if any(np.isnan(hist_status[hist_start_ind:hist_end_ind])):
            logger.info(
                "Missing values in history of {}, skipping status history constraints.".format(
                    structure.symbol
                )
            )
            raise KeyError

        rev_hist_times = np.abs(hist_times[::-1])
        rev_hist_status = hist_status[::-1]

        if not all((hist_status == 0) | (hist_status == 1)):
            raise Exception(
                "Invalid values in history of {} {}".format(
                    type(structure), structure.symbol
                )
            )

        # Force the initial state to match the history
        min_vals = np.full(len(times), 0.0)
        max_vals = np.full(len(times), 1.0)
        min_vals[0] = rev_hist_status[0]
        max_vals[0] = rev_hist_status[0]
        # ToDo: user might want hard constraint on history so need to look further back here as an option

        if self.pumpingstation_history_constraints == "hard":
            hist_bounds = (Timeseries(times, min_vals), Timeseries(times, max_vals))
            if isinstance(structure, Pump):
                cur_bounds = self.__pump_status_bounds[status_sym]
                self.__pump_status_bounds[status_sym] = self.merge_bounds(
                    hist_bounds, cur_bounds
                )
            if isinstance(structure, PumpingStation):
                cur_bounds = self._PumpingStationMixin__station_status_bounds[
                    status_sym
                ]
                self.__station_status_bounds[status_sym] = self.merge_bounds(
                    hist_bounds, cur_bounds
                )
        else:
            if isinstance(structure, Pump):
                self._PumpingStationMixin__psmixin_initial_pump_status[
                    structure.symbol
                ] = rev_hist_status[0]
            if isinstance(structure, PumpingStation):
                self._PumpingStationMixin__psmixin_initial_station_status[
                    structure.symbol
                ] = rev_hist_status[0]

        hist_on = next(
            (rev_hist_times[i] for i, s in enumerate(rev_hist_status) if s == 0),
            rev_hist_times[-1],
        )
        hist_off = next(
            (rev_hist_times[i] for i, s in enumerate(rev_hist_status) if s == 1),
            rev_hist_times[-1],
        )

        return hist_on, hist_off

    def GenerateStatusSymbols(
        self, structure, need_on_sym, need_off_sym, status_sym=None, power_sym=None
    ):
        if need_on_sym > 0:
            sw_on_sym = "{}__switched_on".format(structure.symbol)
        else:
            sw_on_sym = None
        if need_off_sym > 0:
            sw_off_sym = "{}__switched_off".format(structure.symbol)
        else:
            sw_off_sym = None

        # Store variable names
        if isinstance(structure, Pump):
            self._PumpingStationMixin__pump_discrete_symbols.append(status_sym)
            if need_on_sym:
                self._PumpingStationMixin__pump_discrete_symbols.append(sw_on_sym)
            if need_off_sym:
                self._PumpingStationMixin__pump_discrete_symbols.append(sw_off_sym)
            # Store bounds
            self._PumpingStationMixin__pump_status_bounds[status_sym] = (0, 1)
            # If using continuous approach, force status to be one
            constant_inputs = self.constant_inputs(0)
            if structure.semi_continuous is not None:
                continuous = constant_inputs[structure.semi_continuous]
                self._PumpingStationMixin__pump_status_bounds[status_sym] = (
                    continuous,
                    1,
                )
            if need_on_sym:
                self._PumpingStationMixin__pump_status_bounds[sw_on_sym] = (0, 1)
            if need_off_sym:
                self._PumpingStationMixin__pump_status_bounds[sw_off_sym] = (0, 1)
        elif isinstance(structure, PumpingStation):
            self._PumpingStationMixin__station_discrete_symbols.append(status_sym)
            if need_on_sym:
                self._PumpingStationMixin__station_discrete_symbols.append(sw_on_sym)
            if need_off_sym:
                self._PumpingStationMixin__station_discrete_symbols.append(sw_off_sym)
            self._PumpingStationMixin__station_status_bounds[status_sym] = (0, 1)
            if need_on_sym:
                self._PumpingStationMixin__station_status_bounds[sw_on_sym] = (0, 1)
            if need_off_sym:
                self._PumpingStationMixin__station_status_bounds[sw_off_sym] = (0, 1)

        # Generate and store MX variables
        status_mx = MX.sym(status_sym)
        if isinstance(structure, Pump):
            power_mx = MX.sym(power_sym)
        if need_on_sym > 0:
            sw_on_mx = MX.sym(sw_on_sym)
        if need_off_sym > 0:
            sw_off_mx = MX.sym(sw_off_sym)

        self._PumpingStationMixin__pumping_station_mx_path_variables.append(status_mx)
        if isinstance(structure, Pump):
            self._PumpingStationMixin__pumping_station_mx_path_variables.append(
                power_mx
            )
        if need_on_sym:
            self._PumpingStationMixin__pumping_station_mx_path_variables.append(
                sw_on_mx
            )
        if need_off_sym:
            self._PumpingStationMixin__pumping_station_mx_path_variables.append(
                sw_off_mx
            )

        if isinstance(structure, Pump):
            return status_sym, sw_on_sym, sw_off_sym, power_sym
        elif isinstance(structure, PumpingStation):
            return status_sym, sw_on_sym, sw_off_sym


class Pump(_ObjectParameterWrapper):
    """
    Python Pump object as an interface to the
    :cpp:class:`~Deltares::HydraulicStructures::PumpingStation::Pump` object
    in the model.
    """

    def __init__(
        self,
        optimization_problem,
        symbol,
        energy_price_symbol,
        semi_continuous=None,
        status_history=None,
    ):
        super().__init__(optimization_problem)

        self.optimization_problem = optimization_problem
        self.symbol = symbol
        self.energy_price_symbol = energy_price_symbol
        self.semi_continuous = semi_continuous
        self.status_history = status_history

    def discharge(self):
        """
        Get the state corresponding to the pump discharge.

        :returns: `MX` expression of the pump discharge.
        """

        # TODO: We would rather use self.symbol + ".Q" as the control
        # variable, but only top level input variables are allowed. We
        # therefore use the convention that a symbol exists where all dots are
        # replaced with underscores.
        return self.optimization_problem.state(self.symbol.replace(".", "_") + "_Q")

    def head(self):
        """
        Get the state corresponding to the pump head. This depends on the
        ``head_option`` that was specified by the user.

        :returns: `MX` expression of the pump head.
        """
        return self.optimization_problem.state(self.symbol + ".dH")

    @property
    def _need_switched_on_symbol(self):
        return (
            self.minimum_on > 0.0
            or self.start_up_energy > 0.0
            or self.start_up_cost > 0.0
        )

    @property
    def _need_switched_off_symbol(self):
        return (
            self.minimum_off > 0.0
            or self.shut_down_energy > 0.0
            or self.shut_down_cost > 0.0
        )


class Resistance(_ObjectParameterWrapper):
    """
    Python Resistance object as an interface to the
    :cpp:class:`~Deltares::HydraulicStructures::PumpingStation::Resistance`
    object in the model.
    """

    def __init__(self, optimization_problem, symbol):
        super().__init__(optimization_problem)

        self.optimization_problem = optimization_problem
        self.symbol = symbol

    def discharge(self):
        """
        Get the state corresponding to the discharge through the resistance.

        :returns: `MX` expression of the discharge.
        """
        return self.optimization_problem.state(self.symbol + ".HQUp.Q")

    def head_loss(self):
        """
        Get the state corresponding to the head loss over the resistance.

        :returns: `MX` expression of the head loss.
        """

        # Can't we use the dot notation instead, as the two are equated in the
        # Modelica model anyway?
        return self.optimization_problem.state(self.symbol.replace(".", "_") + "_dH")


class PumpingStation(_ObjectParameterWrapper):
    """
    Python PumpingStation object as an interface to the
    :cpp:class:`~Deltares::HydraulicStructures::PumpingStation::PumpingStation` object in the model.
    """

    def __init__(
        self,
        optimization_problem: OptimizationProblem,
        symbol: str,
        pump_symbols: List[str] = None,
        energy_price_symbols: Union[str, List[str]] = None,
        semi_continuous: Union[str, List[str]] = None,
        status_history: Union[str, List[str]] = None,
        **kwargs,
    ):
        """
        Initialize the pumping station object.

        :param optimization_problem:
               :py:class:`~rtctools.optimization.optimization_problem.OptimizationProblem` instance.
        :param symbol: Symbol name of the pumping station in the model.
        :param pump_symbols: Symbol names of the pumps in the pumping station.
        :param energy_price_symbols: String or list of names of the energy price's time series of the pumps
                                     in the pumping station.
        :param semi_continuous: String or list of names of the constant input
                                indicating use of the semi-continuous approach
        :param status_history: String or list of names of the pump status
                               history time series. If string, one can use e.g. "{pump}_status_hist" to map
                               to "pumpingstation1.pump1_status_hist", with "pumpingstation1.pump1" the pump
                               symbol.
        """
        super().__init__(optimization_problem)

        self.optimization_problem = optimization_problem
        self.symbol = symbol
        self.semi_continuous = semi_continuous
        self.status_history = status_history

        # NOTE: We use pump symbols to guarantee the order in which we process
        # the pumps. This is important for e.g. the pump switching matrix,
        # where we need to know what row represents what pump.
        self.pump_symbols = pump_symbols

        self._pumps = None
        self._resistances = None

        if energy_price_symbols is None:
            self.energy_price_symbols = (
                optimization_problem.pumpingstation_energy_price_symbol
            )
        else:
            self.energy_price_symbols = energy_price_symbols

    def discharge(self):
        """
        Get the state corresponding to the station discharge.

        :returns: `MX` expression of the station discharge.
        """

        # TODO: We would rather use self.symbol + ".Q" as the control
        # variable, but only top level input variables are allowed. We
        # therefore use the convention that a symbol exists where all dots are
        # replaced with underscores.
        return self.optimization_problem.state(self.symbol + ".HQUp.Q")

    def pumps(self) -> List[Pump]:
        """
        Get a list of :py:class:`Pump` objects that are part of this pumping station
        in the model.

        :returns: List of :py:class:`Pump` objects.
        """
        if self._pumps is None:
            if self.pump_symbols is None:
                # TODO: Until we are able to guarantee an order in Modelica, we can only come here
                # if the pump switching matrix is all zeros.
                matrix = self.pump_switching_matrix
                if not np.all(matrix == 0):
                    raise Exception(
                        "Automatic finding of pumps not allowed with non-zero switching matrix."
                    )

                _pump_symbols = set()
                for x in self.optimization_problem.parameters(0).keys():
                    m = re.search(r"({}\..+?)\.working_area\[".format(self.symbol), x)
                    if m is None:
                        continue
                    else:
                        _pump_symbols.add(m.group(1))

                self.pump_symbols = sorted(_pump_symbols)

            if isinstance(self.energy_price_symbols, list):
                if len(self.energy_price_symbols) != len(self.pump_symbols):
                    raise Exception(
                        "Each pump in {} must have a corresponding energy price".format(
                            self.symbol
                        )
                    )
                energy_price_symbols = self.energy_price_symbols
            else:
                energy_price_symbols = [self.energy_price_symbols] * len(
                    self.pump_symbols
                )

            if isinstance(self.semi_continuous, list):
                if len(self.semi_continuous) != len(self.pump_symbols):
                    raise Exception(
                        "Each pump in {} must have a corresponding semi continuous series".format(
                            self.symbol
                        )
                    )
                semi_continuous_symbols = self.semi_continuous
            else:
                semi_continuous_symbols = [self.semi_continuous] * len(
                    self.pump_symbols
                )

            if isinstance(self.status_history, list):
                if len(self.status_history) != len(self.pump_symbols):
                    raise Exception(
                        "Each pump in {} must have a corresponding status history series".format(
                            self.symbol
                        )
                    )
                status_history_ts_names = self.status_history
            elif isinstance(self.status_history, str):
                status_history_ts_names = [
                    self.status_history.format(pump=x) for x in self.pump_symbols
                ]
            elif self.status_history is None:
                status_history_ts_names = [None] * len(self.pump_symbols)
            else:
                raise ValueError(
                    "Pump status history must be a list or string, not {}".format(
                        type(self.status_history)
                    )
                )

            self._pumps = [
                Pump(self.optimization_problem, x, e, s, h)
                for x, e, s, h in zip(
                    self.pump_symbols,
                    energy_price_symbols,
                    semi_continuous_symbols,
                    status_history_ts_names,
                )
            ]

        return self._pumps

    def resistances(self) -> List[Resistance]:
        """
        Get a list of :py:class:`Resistance` objects that are part of this pumping station
        in the model.

        :returns: List of :py:class:`Resistance` objects.
        """

        if self._resistances is None:
            _resist_symbols = set()
            for x in self.optimization_problem.parameters(0).keys():
                # TODO: Isn't there a better way to find these components
                # instead of looking for some type of signature (which can
                # change).
                m = re.search(r"({}\..+?)\.C".format(self.symbol), x)
                if m is None:
                    continue
                else:
                    _resist_symbols.add(m.group(1))

            _resist_symbols = sorted(_resist_symbols)

            self._resistances = [
                Resistance(self.optimization_problem, x) for x in _resist_symbols
            ]

        return self._resistances

    @property
    def pump_switching_matrix(self) -> np.ndarray:
        # TODO: Move default values to Modelica, and delete this property
        # method (i.e. let it be handled automatically by __getattr__)
        # TODO: For some reason using super() does not work. Why?
        matrix = _ObjectParameterWrapper.__getattr__(self, "pump_switching_matrix")

        # FIXME: Detect placeholder array for JModelica workaround
        if np.all(matrix == -999):
            matrix = np.tril(np.ones(matrix.shape), -1)
            for i in range(matrix.shape[0]):
                matrix[i, i] = -1 * sum(matrix[i, :])

        # Only lower triangle matrices are allowed
        if not (np.tril(matrix) == matrix).all():
            raise Exception(
                "Switching matrices may only contain a non-zeros in the lower triangle."
            )

        return matrix

    @property
    def pump_switching_constraints(self) -> np.ndarray:
        # TODO: Move default values to Modelica, and delete this property
        # method (i.e. let it be handled automatically by __getattr__)
        constraints = _ObjectParameterWrapper.__getattr__(
            self, "pump_switching_constraints"
        )

        # FIXME: Detect placeholder array for JModelica workaround
        if np.all(constraints == -999):
            constraints = np.transpose(
                [np.zeros(self.n_pumps), list(range(self.n_pumps))]
            )

        return constraints

    @property
    def n_pumps(self) -> int:
        return int(_ObjectParameterWrapper.__getattr__(self, "n_pumps"))

    @property
    def _need_switched_on_symbol(self):
        return (
            self.minimum_on > 0.0
            or self.start_up_energy > 0.0
            or self.start_up_cost > 0.0
        )

    @property
    def _need_switched_off_symbol(self):
        return (
            self.minimum_off > 0.0
            or self.shut_down_energy > 0.0
            or self.shut_down_cost > 0.0
        )


class _MinimizePumpGoalType(Enum):
    ENERGY = 1
    COST = 2


class _MinimizePumpGoal(Goal):
    priority = sys.maxsize
    order = 1

    _type = None

    # NOTE: Based on tests we can typically minimize to about 10% of the
    # original value. Targeting an objective of ~1.0 at the end, that means
    # multiplying the reference nominal (based on previous priorities) by 0.1.
    _dynamic_nominal_scaling = 0.1

    def __init__(
        self,
        use_dynamic_nominal: bool = True,
        exclude_continuous: bool = False,
        *args,
        **kwargs,
    ):
        """
        :param use_dynamic_nominal: Whether to use a dynamically calculated
            nominal based on results of previous priorities.
        """
        super().__init__(*args, **kwargs)

        self.use_dynamic_nominal = use_dynamic_nominal
        self.exclude_continuous = exclude_continuous

    def function(self, o, ensemble_member):
        if self.use_dynamic_nominal:
            priorities, nominals = list(
                zip(*o._psmixin_pump_minimization_nominal[self._type])
            )
            ind = np.searchsorted(priorities, self.priority) - 1
            nominal = nominals[ind] * self._dynamic_nominal_scaling
            # We store the dynamic nominal that we use for debugging purposes.
            # We check/assert that the value does not change, even when this
            # function is called after the priority completes.
            assert (
                not hasattr(self, "_dynamic_nominal")
                or self._dynamic_nominal == nominal
            )
            self._dynamic_nominal = nominal

            if self.function_nominal != 1.0:
                raise Exception(
                    "The minimization goal's function_nominal has to be 1.0 when using dynamic scaling."
                )
        else:
            self._dynamic_nominal = 1.0

        assert self._dynamic_nominal > 0.0
        assert np.isfinite(self._dynamic_nominal)

        costs = 0.0

        times = o.times()

        constant_inputs = o.constant_inputs(ensemble_member)

        for ps in o.pumping_stations():
            for p in ps.pumps():
                if p.semi_continuous is not None and self.exclude_continuous:
                    assert np.array_equal(
                        constant_inputs[p.semi_continuous].times, times
                    )
                    continuous = constant_inputs[p.semi_continuous].values
                else:
                    continuous = np.full_like(times, 0.0)

                for ts, tf, cont in zip(times[:-1], times[1:], continuous[1:]):
                    tstep = tf - ts

                    if self._type == _MinimizePumpGoalType.COST:
                        price = o.timeseries_at(p.energy_price_symbol, tf)
                        # Zero or negative price values will lead to wrong results
                        if price < 0:
                            raise Exception(
                                "Price for pump {} at t = {} s is negative.".format(
                                    p.symbol, tf
                                )
                            )

                    elif self._type == _MinimizePumpGoalType.ENERGY:
                        price = 1.0
                    else:
                        raise NotImplementedError("Unknown minimization goal type.")

                    # TODO: Not pretty to use the same formatting again
                    # Pump power
                    costs += (
                        tstep
                        * o.state_at("{}__power".format(p.symbol), tf)
                        * price
                        * (1 - cont)
                    )

                    # Start-up energy
                    if p.start_up_energy > 0.0:
                        costs += (
                            o.state_at("{}__switched_on".format(p.symbol), tf)
                            * p.start_up_energy
                            * price
                        )

                    # Shut-down energy
                    if p.shut_down_energy > 0.0:
                        costs += (
                            o.state_at("{}__switched_off".format(p.symbol), tf)
                            * p.shut_down_energy
                            * price
                        )

                    if self._type == _MinimizePumpGoalType.ENERGY:
                        continue
                    else:
                        assert self._type == _MinimizePumpGoalType.COST

                    # Fixed start-up costs (other than energy)
                    if p.start_up_cost > 0.0:
                        costs += (
                            o.state_at("{}__switched_on".format(p.symbol), tf)
                            * p.start_up_cost
                        )

                    # Fixed shut-down costs (other than energy)
                    if p.shut_down_cost > 0.0:
                        costs += (
                            o.state_at("{}__switched_off".format(p.symbol), tf)
                            * p.shut_down_cost
                        )
        return costs / self._dynamic_nominal


class MinimizePumpCostGoal(_MinimizePumpGoal):
    """
    Goal that minimizes overall energy costs.

    Loops over all pumps in all pumping stations, integrating all
    instantaneous pump operating costs (and any start-up and shut-down
    costs/energy) in the optimization window into one objective value.

    :cvar function_nominal: Nominal value of needed for scaling. Guessed
       automatically based on the power range of all pumps.
    :cvar priority: Priority of this goal. Default is ``sys.maxsize``.
    """

    _type = _MinimizePumpGoalType.COST


class MinimizePumpEnergyGoal(_MinimizePumpGoal):
    """
    Goal that minimizes overall energy consumption.

    Loops over all pumps in all pumping stations, integrating all
    instantaneous pump powers (and any start-up and shut-down energy) in the
    optimization window into one objective value.

    :cvar function_nominal: Nominal value of needed for scaling. Guessed
       automatically based on the power range of all pumps.
    :cvar priority: Priority of this goal. Default is ``sys.maxsize``.
    """

    _type = _MinimizePumpGoalType.ENERGY


class StructureStatusGoal(Goal):
    order = 1

    def __init__(
        self,
        optimization_problem,
        structure,
        minimum_on=None,
        minimum_off=None,
        initial_status=True,
        horizon_status=False,
        priority=1,
        weight=1,
        function_nominal=1,
    ):
        if minimum_on is None:
            minimum_on = structure.minimum_on
        if minimum_off is None:
            minimum_off = structure.minimum_off

        if minimum_on > 0 and not structure._need_switched_on_symbol:
            raise Exception(
                "StructureStatusGoal: A switched on symbol is needed when using minimum on time. "
                "For example, set the parameter {}.minimum_on to a small, but non-zero, value.".format(
                    structure.symbol
                )
            )

        if minimum_off > 0 and not structure._need_switched_off_symbol:
            raise Exception(
                "StructureStatusGoal: A switched off symbol is needed when using minimum off time. "
                "For example, set the parameter {}.minimum_off to a small, but non-zero, value.".format(
                    structure.symbol
                )
            )

        self.structure = structure
        self.minimum_on = minimum_on
        self.minimum_off = minimum_off
        self.initial_status = initial_status
        self.horizon_status = horizon_status
        self.priority = priority
        self.weight = weight
        self.function_nominal = function_nominal

        initial_constraints = []
        horizon_constraints = []

        if initial_status:
            initial_constraints = (
                optimization_problem._psmixin_initial_status_constraints(
                    structure, minimum_on, minimum_off
                )
            )
        if horizon_status:
            horizon_constraints = (
                optimization_problem._psmixin_horizon_status_constraints(
                    structure, minimum_on, minimum_off
                )
            )

        total_size = 0
        self._constraints = []

        if initial_status:
            total_size += sum(len(x[1]) for x in initial_constraints)
            self._constraints.extend(initial_constraints)

            if isinstance(structure, Pump):
                try:
                    init_val = optimization_problem._PumpingStationMixin__psmixin_initial_pump_status[
                        structure.symbol
                    ]
                    status_sym, *_ = (
                        optimization_problem._PumpingStationMixin__pump_status_pairs[
                            structure.symbol
                        ]
                    )

                    def f(o, status_sym=status_sym):
                        return o.state_at(status_sym, o.initial_time)

                    self._constraints.append((f, init_val, init_val, -1, 2))
                    total_size += 1
                except KeyError:
                    pass
            if isinstance(structure, PumpingStation):
                try:
                    init_val = optimization_problem._PumpingStationMixin__psmixin_initial_station_status[
                        structure.symbol
                    ]
                    status_sym, *_ = (
                        optimization_problem._PumpingStationMixin__station_status_pairs[
                            structure.symbol
                        ]
                    )

                    def f(o, status_sym=status_sym):
                        return o.state_at(status_sym, o.initial_time)

                    self._constraints.append((f, init_val, init_val, -1, 2))
                    total_size += 1
                except KeyError:
                    pass

        if horizon_status:
            total_size += len(horizon_constraints)
            self._constraints.extend(horizon_constraints)

        if total_size > 0:
            self.size = total_size

            self.target_min = np.hstack([x[1] for x in self._constraints])
            self.target_max = np.hstack([x[2] for x in self._constraints])

            fr_min = np.hstack([x[3] for x in self._constraints])
            fr_max = np.hstack([x[4] for x in self._constraints])

            if self.size == 1:
                self.target_min = float(self.target_min[0])
                self.target_max = float(self.target_max[0])

                fr_min = float(fr_min[0])
                fr_max = float(fr_max[0])

            self.function_range = (fr_min, fr_max)

    def function(self, optimization_problem, ensemble_member):
        if self._constraints:
            return ca.vertcat(*[x[0](optimization_problem) for x in self._constraints])
        else:
            return ca.MX(0)


class PumpStatusGoal(StructureStatusGoal):
    pass


class _InvalidCacheError(Exception):
    pass


class PumpingStationMixin(OptimizationProblem, CommonStructureSwitchFunctions):
    """
    Adds handling of PumpingStation objects in your model to your optimization
    problem.

    Relevant parameters and variables are read from the model, and from this
    data a set of constraints and objectives are automatically generated to
    minimize cost.

    If historical data regarding the status of the pumps is provided, this
    information is used to ensure that the minimum amount of time a pump must
    be on / off is respected.
    """

    # TODO: Vijzels are different in that the pump head is just the upstream
    # head, and the discharge/working area is a non-smooth function that does
    # not fit a polynomial well. How to handle these?
    _pumping_station_mx_path_variables = []

    # In the post() routine we check if the non-linear equality constraints
    # (e.g. pump power and resistance head loss) minimized to their equality
    # constraint. A warning is raised if the relative error for any particular
    # time step exceeds this value:
    _pumpingstation_ineq_relative_error = 1e-4

    # In addition to the relative error check, we also check the absolute
    # error. Note that the absolute tolerance applies to the scaled power. In
    # other words, the value entered here will be multiplied with the pump's
    # maximum power before checking.
    _pumpingstation_absolute_error = 1e-5

    #: Use pickle to cache the HQ subproblems that are solved.
    pumpingstation_cache_hq_subproblem = True

    #: Energy price symbol to use if no symbol specified per pumping station
    #: or per pump.
    pumpingstation_energy_price_symbol = "energy_price"

    #: How pump status history constraints should be enforced. Either as a
    #: hard constraint ('hard'), or by using the PumpStatusGoal ('soft' /
    #: None). Note that the latter requires 'keep_soft_constraints' to be set.
    pumpingstation_history_constraints = "hard"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # We also want to output a bunch of variables that are calculated in
        # post(), which we return in the extract_results() call. We store them
        # in this dictionary.
        self.__additional_results = OrderedDict()
        self.__additional_output_variables = set()

        assert "model_folder" in kwargs
        self.__model_folder = kwargs["model_folder"]

        self.__hq_subproblem_cache = OrderedDict()
        self.__hq_subproblem_cache_path = os.path.join(
            kwargs["model_folder"], "_hq_subproblem_cache.pickle"
        )

    def pre(self):
        super().pre()

        # TODO: Is a reset necessary? It is a bit double with the statements above.
        self.__pumping_station_mx_path_variables = []
        self.__pump_discrete_symbols = []
        self.__station_discrete_symbols = []
        self.__pump_status_bounds = OrderedDict()
        self.__station_status_bounds = OrderedDict()
        self.__psmixin_initial_pump_status = OrderedDict()
        self.__psmixin_initial_station_status = OrderedDict()
        self.__pump_power_bounds = OrderedDict()
        self._psmixin_pump_discharge_bounds = OrderedDict()
        self.__pump_status_pairs = OrderedDict()
        self.__station_status_pairs = OrderedDict()
        self.__pump_power_range_on = OrderedDict()

        self._psmixin_pump_working_area_head_range = OrderedDict()
        self._psmixin_pump_extended_working_area_head_range = OrderedDict()

        self.__head_range_up = OrderedDict()
        self.__head_range_down = OrderedDict()
        self.__head_range_diff = OrderedDict()

        self._psmixin_pump_minimization_nominal = OrderedDict()
        self.__pump_minimization_function_values = []

        self.__pump_status_hist = {}
        self.__station_status_hist = {}

        # Uses the same mapping as Pump.head_option for easy access
        self._psmixin_head_range = {
            -1: self.__head_range_up,
            0: self.__head_range_diff,
            1: self.__head_range_down,
            2: self.__head_range_diff,
        }

        # Check validity of HQ subproblem cache. We check if _any_ file other
        # than ourself in the model folder is newer.
        try:
            if not self.pumpingstation_cache_hq_subproblem:
                raise _InvalidCacheError("Caching disabled")

            # Calling getmtime() on a file that does not exist returns a vague
            # platform varying Exception. Try to prevent that from happening.
            if not os.path.exists(self.__hq_subproblem_cache_path):
                raise _InvalidCacheError("Cache file does not exist")

            cache_mtime = os.path.getmtime(self.__hq_subproblem_cache_path)
            cache_abspath = os.path.abspath(self.__hq_subproblem_cache_path)
            for root, _dir, files in os.walk(self.__model_folder):
                for f in files:
                    f_abspath = os.path.abspath(os.path.join(root, f))
                    if (
                        cache_abspath != f_abspath
                        and os.path.getmtime(f_abspath) > cache_mtime
                    ):
                        raise _InvalidCacheError("Cache no longer valid")

            with open(self.__hq_subproblem_cache_path, "rb") as file:
                # Load subproblem.
                try:
                    self.__hq_subproblem_cache = pickle.load(file)
                except RuntimeError as error:
                    logger.warning(f"Loading cache file {file} failed: {error}.")

        except _InvalidCacheError:
            self.__hq_subproblem_cache = OrderedDict()

        # Automatic deriviation of maximum head over pumping station. Note
        # that it is possible that different pumps use different head
        # definitions (up, down, diff) for the maximum head of the enclosing
        # pumping station.
        for ps in self.pumping_stations():
            # dependent on head option
            # assume all pumps have the same head option if option 2 is used
            for p in ps.pumps():
                head_option = p.head_option
                break
            if head_option in [-1, 0, 1]:
                symbol_up = ps.symbol + ".HQUp.H"
                symbol_down = ps.symbol + ".HQDown.H"
            elif head_option == 2:
                symbol_up = ps.symbol + ".HW"
                symbol_down = ps.symbol + ".HQDown.H"

            # Use bounds if specified, otherwise try deriving bounds from time series
            for hr, s in [
                (self.__head_range_up, symbol_up),
                (self.__head_range_down, symbol_down),
            ]:
                m, M = self.bounds().get(s, [None, None])

                if isinstance(m, Timeseries):
                    m = min(m.values)
                if isinstance(M, Timeseries):
                    M = min(M.values)

                try:
                    ts = self.get_timeseries(s)
                    canonical_state, sign = self.alias_relation.canonical_signed(s)

                    # Discard history values for head estimation.
                    ts_values = np.array(
                        [v for t, v in zip(ts.times, ts.values) if t in self.times()]
                    )

                    # Attempt to avoid using time series used as initial conditions only by checking for NaN
                    if m is None or not np.isfinite(m) and not any(np.isnan(ts_values)):
                        m = min(ts_values)
                        logger.info(
                            'Using {} value "{}" in time series "{}" as lower bound for "{}".'.format(
                                "minimum" if sign == 1 else "maximum",
                                m,
                                canonical_state,
                                s,
                            )
                        )
                    if M is None or not np.isfinite(M) and not any(np.isnan(ts_values)):
                        M = max(ts_values)
                        logger.info(
                            'Using {} value "{}" in time series "{}" as upper bound for "{}".'.format(
                                "maximum" if sign == 1 else "minimum",
                                M,
                                canonical_state,
                                s,
                            )
                        )
                except KeyError:
                    # Time series does not exist
                    pass

                if m is None or M is None or not np.isfinite(m) or not np.isfinite(M):
                    raise Exception(
                        "Specify (finite) bounds or time series for '{}', currently found {}".format(
                            s, (m, M)
                        )
                    )

                hr[ps.symbol] = (m, M)

            self.__head_range_diff[ps.symbol] = [
                a - b
                for a, b in zip(
                    self.__head_range_down[ps.symbol],
                    reversed(self.__head_range_up[ps.symbol]),
                )
            ]

        # Automatic derivation of discharge and power if not specified by user
        for ps in self.pumping_stations():
            for p in ps.pumps():
                discharge_sym = p.symbol.replace(".", "_") + "_Q"
                head_sym = p.symbol + "_head"
                power_sym = "{}__power".format(p.symbol)

                # Find the maximum discharge and head in the working area
                Q = SX.sym("Q")
                H = SX.sym("H")

                hr = self._psmixin_head_range[p.head_option][ps.symbol]

                _, (_, q_max) = self.__solve_working_area_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, -1 * Q
                )
                self._psmixin_pump_discharge_bounds[discharge_sym] = (0.0, q_max)

                _, (h_max, _) = self.__solve_working_area_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, -1 * H
                )
                self._psmixin_pump_working_area_head_range[head_sym] = (0.0, h_max)
                # TODO: determine h_min, or do we count on it being 0.0?

                _, (h_max, _) = self.__solve_working_area_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, -1 * H, 0
                )
                _, (h_min, _) = self.__solve_working_area_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, H, 0
                )
                self._psmixin_pump_extended_working_area_head_range[head_sym] = (
                    h_min,
                    h_max,
                )

                # Recalculate the maximum discharge, but now for the
                # _extended_ working area, so that we can use it for the
                # maximum power calculation below.
                _, (_, q_max) = self.__solve_working_area_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, -1 * Q, 0
                )

                # Lower power bound (when on)
                coeffs = p.power_coefficients
                power_functions = self.__power_functions(H, Q, coeffs)

                min_power, _ = self.__solve_power_subproblem(
                    p.working_area, p.working_area_direction, hr, H, Q, power_functions
                )

                # For a monotonically increasing convex function we can find
                # the maximum by checking all vertices of the polygon. We do
                # this not on the working area (which is not a polygon), but
                # on an enclosing square. Note that we do this on the H and Q
                # ranges of the _extended_ working area.
                max_power = 0.0

                for h, q in [(h_min, 0), (h_max, 0), (h_min, q_max), (h_max, q_max)]:
                    powers = self.__power_functions(h, q, coeffs)
                    max_power = np.maximum(powers, max_power)

                self.__pump_power_bounds[power_sym] = (0.0, max(max_power))
                self.__pump_power_range_on[power_sym] = (min_power, max_power)

        # Convexity and increasing-with-H check of power coefficients
        Q = SX.sym("Q")
        H = SX.sym("H")
        X = vertcat(Q, H)

        for ps in self.pumping_stations():
            for p in ps.pumps():
                discharge_sym = p.symbol.replace(".", "_") + "_Q"
                head_sym = p.symbol + "_head"

                coeffs = p.power_coefficients
                powers = self.__power_functions(H, Q, coeffs)

                for power_i, power in enumerate(powers):
                    # Check if power is increasing with H (only necessary if there are resistances)
                    if ps.resistances():
                        sx_jac = jacobian(power, H)

                        sx_hess = hessian(sx_jac, X)[0]
                        sx_determinant = det(sx_hess)
                        sx_trace = trace(sx_hess)

                        # CasADi returns NaN if the expression is still a function of H and/or Q
                        determinant = float(sx_determinant)
                        trace_calculated = float(sx_trace)

                        if np.isnan(determinant):
                            logger.warning(
                                "Cannot determine monotonicity in H of power equation "
                                'at index {} of pump "{}".'.format(
                                    power_i + 1, p.symbol
                                )
                            )
                        elif determinant < 0.0 or trace_calculated < 0.0:
                            # Concave function of which we are trying to find the minimum --> use enclosing rectangle
                            h_min, h_max = self._psmixin_pump_working_area_head_range[
                                head_sym
                            ]
                            q_min, q_max = self._psmixin_pump_discharge_bounds[
                                discharge_sym
                            ]

                            min_jac = np.inf

                            for h, q in [
                                (h_min, q_min),
                                (h_max, q_min),
                                (h_min, q_max),
                                (h_max, q_max),
                            ]:
                                cur_jac = float(
                                    substitute(substitute(sx_jac, H, h), Q, q)
                                )
                                min_jac = min(cur_jac, min_jac)

                            if min_jac < 0.0:
                                logger.warning(
                                    'Power equation at index {} of pump "{}" likely not increasing '
                                    "with H in working area.".format(power_i, p.symbol)
                                )
                        else:
                            hr = self._psmixin_head_range[p.head_option][ps.symbol]
                            # We require convexity on the working area (i.e. when pump is on)
                            minimum_jac, _ = self.__solve_working_area_subproblem(
                                p.working_area,
                                p.working_area_direction,
                                hr,
                                H,
                                Q,
                                sx_jac,
                            )

                            if minimum_jac < 0.0:
                                logger.error(
                                    'Power equation at index {} of pump "{}" is not increasing '
                                    "with H in working area.".format(power_i, p.symbol)
                                )

                    # Convexity check:
                    sx_hess = hessian(power, X)[0]
                    sx_determinant = det(sx_hess)
                    sx_trace = trace(sx_hess)

                    # CasADi returns NaN if the expression is still a function of H and/or Q
                    determinant = float(sx_determinant)
                    trace_calculated = float(sx_trace)

                    if (not np.isnan(determinant) and determinant < 0.0) or (
                        not np.isnan(trace_calculated) and trace_calculated < 0.0
                    ):
                        logger.error(
                            'Non-convex power relationship specified for pump "{}" '
                            "at power equation index {}.".format(p.symbol, power_i)
                        )
                    elif np.isnan(determinant):
                        # The determinant is an expression of H and Q. Check if it
                        # is a convex expression, and if so, find the minimum.
                        sx_det_hess = hessian(sx_determinant, X)[0]
                        sx_det_determinant = det(sx_det_hess)

                        det_determinant = float(sx_det_determinant)

                        if not np.isnan(det_determinant):
                            hr = self._psmixin_head_range[p.head_option][ps.symbol]
                            # We require convexity on the expanded working area, i.e when pump is off
                            minimum_determinant, _ = (
                                self.__solve_working_area_subproblem(
                                    p.working_area,
                                    p.working_area_direction,
                                    hr,
                                    H,
                                    Q,
                                    sx_determinant,
                                    0,
                                )
                            )
                            minimum_trace, _ = self.__solve_working_area_subproblem(
                                p.working_area,
                                p.working_area_direction,
                                hr,
                                H,
                                Q,
                                sx_trace,
                                0,
                            )

                            if minimum_determinant < 0.0 or minimum_trace < 0.0:
                                logger.error(
                                    'Power for pump "{}" is not convex over extended working area'.format(
                                        p.symbol
                                    )
                                )

                            # For correct maximum power estimation, we also require
                            # convexity on the rectangular region:
                            # Q \in [0, max_q]
                            # H \in [0, max_h]
                            # where  max_q and max_h are the maximum discharge
                            # and pump head when the pump is _on_.
                            bnds = {
                                Q.name(): self._psmixin_pump_discharge_bounds[
                                    discharge_sym
                                ],
                                H.name(): self._psmixin_pump_working_area_head_range[
                                    head_sym
                                ],
                            }
                            minimum_determinant, _ = self.__solve_hq_subproblem(
                                H, Q, sx_determinant, bounds=bnds
                            )
                            minimum_trace, _ = self.__solve_hq_subproblem(
                                H, Q, sx_trace, bounds=bnds
                            )

                            if minimum_determinant < 0.0 or minimum_trace < 0.0:
                                logger.error(
                                    'Power equation at index {} for pump "{}" is not convex '
                                    "over max. head/discharge range.".format(
                                        power_i, p.symbol
                                    )
                                )
                        else:
                            logger.warning(
                                "Cannot detect convexity of power coefficients "
                                'at index {} of pump "{}".'.format(power_i, p.symbol)
                            )
                    else:
                        # Positive determinant --> convex
                        continue

        for ps in self.pumping_stations():
            # Add discharge symbols for each of the pumps in this pumping station
            # TODO: Add on/off + switched on/off symbols for each pump as well.
            #       Make sure that this is accompanied by also adding constraints on
            #       pump 2 only being able to switch on (or be on?) when pump 1 is on.
            station_status_sym = "{}__status".format(ps.symbol)
            station_need_on_sym = False
            station_need_off_sym = False
            for p in ps.pumps():
                if p.status_history is not None:
                    station_hist_values = np.full_like(
                        self.get_timeseries(p.status_history).values, 0.0
                    )
            for p in ps.pumps():
                # Define symbol names
                status_sym = "{}__status".format(p.symbol)
                power_sym = "{}__power".format(p.symbol)

                need_on_sym = p._need_switched_on_symbol
                station_need_on_sym += need_on_sym
                need_off_sym = p._need_switched_off_symbol
                station_need_off_sym += need_off_sym
                # Generate status symbols
                (status_sym, sw_on_sym, sw_off_sym, power_sym) = (
                    CommonStructureSwitchFunctions.GenerateStatusSymbols(
                        self,
                        p,
                        need_on_sym,
                        need_off_sym,
                        status_sym=status_sym,
                        power_sym=power_sym,
                    )
                )

                # Store all symbols together
                self.__pump_status_pairs[p.symbol] = (
                    status_sym,
                    sw_on_sym,
                    sw_off_sym,
                    power_sym,
                )

                # Pump status history on/off
                times = self.times()
                if p.status_history is not None:
                    try:
                        pump_history = self.get_timeseries(p.status_history)
                        station_hist_values += pump_history.values
                    except KeyError as e:
                        raise KeyError(
                            "History for pump '{}' specified, but Timeseries '{}' was not found.".format(
                                p.symbol, p.status_history
                            )
                        ) from e
                    hist_on, hist_off = CommonStructureSwitchFunctions.StructureHistory(
                        self, p, times, status_sym, pump_history=pump_history
                    )
                else:
                    hist_on = 0.0
                    hist_off = 0.0

                self.__pump_status_hist[p.symbol] = (hist_on, hist_off)

            # station level
            (
                station_status_sym,
                station_sw_on_sym,
                station_sw_off_sym,
            ) = CommonStructureSwitchFunctions.GenerateStatusSymbols(
                self,
                ps,
                station_need_on_sym,
                station_need_off_sym,
                status_sym=station_status_sym,
            )
            # Store all symbols together
            self.__station_status_pairs[ps.symbol] = (
                station_status_sym,
                station_sw_on_sym,
                station_sw_off_sym,
            )
            if ps.status_history is not None:
                # ensure station history status is 0 or 1.
                if not all(v in (0, 1) for v in station_hist_values):
                    if all(v in (0, 1) or np.isnan(v) for v in station_hist_values):
                        # if all values are 0, 1 or nan, we can just fill nan with 0
                        logger.info(
                            f"Status history for pumping station '{ps.symbol}' contains NaN values. "
                            "Filling NaNs with 0."
                        )
                        station_hist_values = np.array(
                            [0 if np.isnan(v) else v for v in station_hist_values]
                        )
                    else:
                        logger.warning(
                            f"Status history for pumping station '{ps.symbol}' contains "
                            "values other than 0 or 1. Converting to binary status history "
                            "by converting all values > 0 to 1 and all others to 0."
                        )
                        station_hist_values = np.array(
                            [1 if v > 0 else 0 for v in station_hist_values]
                        )
                try:
                    hist_on, hist_off = CommonStructureSwitchFunctions.StructureHistory(
                        self,
                        ps,
                        times,
                        station_status_sym,
                        station_history=Timeseries(
                            self.io.times_sec, station_hist_values
                        ),
                    )
                except AttributeError:
                    continue

            else:
                hist_on = 0.0
                hist_off = 0.0

            self.__station_status_hist[ps.symbol] = (hist_on, hist_off)

        # If a minimization goal is used at the first priority, we still need
        # a value for the dynamic function nominal. This is done based on the
        # power range.
        minimization_nominals = defaultdict(int)  # Default of zero

        for ps in self.pumping_stations():
            for p in ps.pumps():
                avg_pump_power = (
                    max(self.__pump_power_range_on[p.symbol + "__power"][1]) / 2
                )

                # Energy
                minimization_nominals[_MinimizePumpGoalType.ENERGY] += avg_pump_power

                # Cost
                ts = self.get_timeseries(p.energy_price_symbol)
                price = self.interpolate(self.times(), ts.times, ts.values)

                # price_key = "{}_{}".format(p.symbol, "price")
                # self.__additional_output_variables.add(price_key)
                # self.__additional_results[price_key] = price

                avg_pump_energy_price = np.average(price)
                pump_nominal = avg_pump_energy_price * avg_pump_power
                minimization_nominals[_MinimizePumpGoalType.COST] += pump_nominal

        for k, v in minimization_nominals.items():
            self._psmixin_pump_minimization_nominal[k] = [
                (-np.inf, v * (self.times()[-1] - self.times()[0]))
            ]

        # Store cache to disk
        if self.pumpingstation_cache_hq_subproblem:
            with open(self.__hq_subproblem_cache_path, "wb") as f:
                pickle.dump(self.__hq_subproblem_cache, f)

    def __power_functions(self, head, discharge, coeffs, status=1):
        powers = []
        for i in range(coeffs.shape[0]):
            power = self.__power_function(head, discharge, coeffs[i], status)
            powers.append(power)
        return powers

    def __power_function(self, head, discharge, coeffs, status=1):
        power = 0.0
        for i in range(coeffs.shape[0]):
            for j in range(coeffs.shape[1]):
                power += coeffs[i, j] * head**i * discharge**j * status
        return power

    def spec_energy_functions(self, head, discharge, coeffs, status=1):
        speces = []
        for i in range(coeffs.shape[0]):
            spece = self.spec_energy_function(head, discharge, coeffs[i], status)
            speces.append(spece)
        return speces

    def spec_energy_function(self, head, discharge, coeffs, status=1):
        spece = 0.0
        for i in range(coeffs.shape[0]):
            for j in range(coeffs.shape[1]):
                # unit: kWh/1000m3
                spece += (
                    coeffs[i, j]
                    * head**i
                    * discharge**j
                    * status
                    / (discharge * 3600)
                    * 1000
                )
        return spece

    def _psmixin_initial_status_constraints(
        self, structure, minimum_on=0, minimum_off=0
    ):
        initial_constraints = []

        if isinstance(structure, Pump):
            hist_on, hist_off = self.__pump_status_hist[structure.symbol]
            status_sym, sw_on_sym, sw_off_sym, _ = self.__pump_status_pairs[
                structure.symbol
            ]
        elif isinstance(structure, PumpingStation):
            hist_on, hist_off = self.__station_status_hist[structure.symbol]
            status_sym, sw_on_sym, sw_off_sym = self.__station_status_pairs[
                structure.symbol
            ]
        else:
            logger.error(
                "We currently do not support structures of type {}".format(
                    type(structure)
                )
            )

        times = self.times(status_sym)

        # Force status at start of horizon based on history
        if minimum_on > 0.0:
            if hist_on > 0.0 and minimum_on > hist_on:
                min_on = minimum_on - hist_on

                max_ind = min(np.searchsorted(times, min_on) + 1, len(times))

                def _f(o, status_sym=status_sym, max_ind=max_ind):
                    return o.state_vector(status_sym)[1:max_ind]

                initial_constraints.append(
                    (
                        _f,
                        np.ones(max_ind - 1),
                        np.full(max_ind - 1, np.inf),
                        np.zeros(max_ind - 1),
                        np.ones(max_ind - 1),
                    )
                )

        if minimum_off > 0.0:
            if hist_off > 0.0 and minimum_off > hist_off:
                min_off = minimum_off - hist_off

                max_ind = min(np.searchsorted(times, min_off) + 1, len(times))

                def _f(o, status_sym=status_sym, max_ind=max_ind):
                    return o.state_vector(status_sym)[1:max_ind]

                initial_constraints.append(
                    (
                        _f,
                        np.full(max_ind - 1, -np.inf),
                        np.zeros(max_ind - 1),
                        np.zeros(max_ind - 1),
                        np.ones(max_ind - 1),
                    )
                )
        return initial_constraints

    def _psmixin_horizon_status_constraints(
        self, structure, minimum_on=0, minimum_off=0
    ):
        horizon_constraints = []
        if isinstance(structure, Pump):
            status_sym, sw_on_sym, sw_off_sym, _ = self.__pump_status_pairs[
                structure.symbol
            ]
        elif isinstance(structure, PumpingStation):
            status_sym, sw_on_sym, sw_off_sym = self.__station_status_pairs[
                structure.symbol
            ]
        else:
            logger.error(
                "We currently do not support structures of type {}".format(
                    type(structure)
                )
            )
        times = self.times(status_sym)
        num_tsteps = len(times)

        # Force pump status minimum on/off time throughout horizon
        if minimum_on > 0.0:
            # Figure out at what times the structure status symbol should be on
            # when the structure switches on at a particular timestep
            status_end_inds = np.searchsorted(times, times[:-1] + minimum_on).tolist()
            status_range_inds = [
                np.arange(i + 1, min((e + 1), len(times))).tolist()
                for i, e in enumerate(status_end_inds)
            ]
            assert isinstance(status_range_inds[0][0], int), (
                "Indexing CasADi symbols only fast with list of _Python_ ints"
            )

            for i in range(1, num_tsteps):
                cur_status_inds = status_range_inds[i - 1]
                if len(cur_status_inds) <= 1:
                    # Trivial constraint
                    continue

                def _f(
                    o,
                    status_sym=status_sym,
                    cur_status_inds=cur_status_inds,
                    sw_on_sym=sw_on_sym,
                    i=i,
                ):
                    return (
                        ca.sum1(o.state_vector(status_sym)[cur_status_inds])
                        - len(cur_status_inds) * o.state_vector(sw_on_sym)[i]
                    )

                horizon_constraints.append(
                    (_f, 0.0, np.inf, -len(cur_status_inds), 0.0)
                )

        if minimum_off > 0.0:
            # Figure out at what times the pump status symbol should be off
            # when the pump switches off at a particular timestep
            status_end_inds = np.searchsorted(times, times[:-1] + minimum_off).tolist()
            status_range_inds = [
                np.arange(i + 1, min((e + 1), len(times))).tolist()
                for i, e in enumerate(status_end_inds)
            ]
            assert isinstance(status_range_inds[0][0], int), (
                "Indexing CasADi symbols only fast with list of _Python_ ints"
            )

            for i in range(1, num_tsteps):
                cur_status_inds = status_range_inds[i - 1]
                if len(cur_status_inds) <= 1:
                    # Trivial constraint
                    continue

                def _f(
                    o,
                    status_sym=status_sym,
                    cur_status_inds=cur_status_inds,
                    sw_off_sym=sw_off_sym,
                    i=i,
                ):
                    return (
                        ca.sum1(1 - o.state_vector(status_sym)[cur_status_inds])
                        - len(cur_status_inds) * o.state_vector(sw_off_sym)[i]
                    )

                horizon_constraints.append(
                    (_f, 0.0, np.inf, -len(cur_status_inds), 0.0)
                )

        return horizon_constraints

    @abstractmethod
    def pumping_stations(self) -> List[PumpingStation]:
        """
        User problem returns list of :class:`PumpingStation` objects.

        :returns: A list of pumping stations.
        """
        raise NotImplementedError()

    def constraints(self, ensemble_member):
        constraints = super().constraints(ensemble_member)

        for ps in self.pumping_stations():
            pumps_index = 0.0
            for p in ps.pumps():
                status_sym, sw_on_sym, sw_off_sym, _ = self.__pump_status_pairs[
                    p.symbol
                ]

                d_tm1 = self.state_vector(status_sym)[:-1]
                d_t = self.state_vector(status_sym)[1:]
                d_diff = d_t - d_tm1

                if sw_on_sym is not None:
                    # On/off symbols should be collocated at same times as
                    # status symbols. It does not make sense to have more of
                    # one than the others.
                    assert np.array_equal(self.times(status_sym), self.times(sw_on_sym))

                    x = self.state_vector(sw_on_sym)[1:]

                    # x is 1 if and only if pump switched on (else 0)
                    constraints.append((d_t - x, 0, inf))
                    constraints.append((x - d_diff, 0, inf))
                    constraints.append((1 - d_tm1 - x, 0, inf))

                if sw_off_sym is not None:
                    assert np.array_equal(
                        self.times(status_sym), self.times(sw_off_sym)
                    )

                    y = self.state_vector(sw_off_sym)[1:]

                    # y is 1 if and only if pump switched off (else 0)
                    constraints.append((d_tm1 - y, 0, inf))
                    constraints.append((y + d_diff, 0, inf))
                    constraints.append((1 - d_t - y, 0, inf))

                if self.pumpingstation_history_constraints == "hard":
                    constraints.extend(
                        [
                            (x[0](self), *x[1:3])
                            for x in self._psmixin_initial_status_constraints(
                                p, p.minimum_on, p.minimum_off
                            )
                        ]
                    )
                try:
                    ii = [y[0] for y in self.PumpsCardinality].index(self.priority)
                    matrix_indicator = self.PumpsCardinality[ii][1]
                except Exception:
                    try:
                        matrix_indicator = self.PumpsCardinality[0][1]
                    except Exception:
                        matrix_indicator = True
                if pumps_index == 0.0 or matrix_indicator:
                    # Always force constraints for minimum on/off time throughout the horizon
                    constraints.extend(
                        [
                            (x[0](self), *x[1:3])
                            for x in self._psmixin_horizon_status_constraints(
                                p, p.minimum_on, p.minimum_off
                            )
                        ]
                    )

                pumps_index += 1.0
            # constraints for station here.
            if ps.status_history is not None:
                status_sym, sw_on_sym, sw_off_sym = self.__station_status_pairs[
                    ps.symbol
                ]

                d_tm1 = self.state_vector(status_sym)[:-1]
                d_t = self.state_vector(status_sym)[1:]
                d_diff = d_t - d_tm1

                if sw_on_sym is not None:
                    assert np.array_equal(self.times(status_sym), self.times(sw_on_sym))

                    x = self.state_vector(sw_on_sym)[1:]
                    # x is 1 if and only if station switched on (else 0)
                    constraints.append((d_t - x, 0, inf))
                    constraints.append((x - d_diff, 0, inf))
                    constraints.append((1 - d_tm1 - x, 0, inf))

                if sw_off_sym is not None:
                    assert np.array_equal(
                        self.times(status_sym), self.times(sw_off_sym)
                    )
                    y = self.state_vector(sw_off_sym)[1:]
                    # y is 1 if and only if station switched off (else 0)
                    constraints.append((d_tm1 - y, 0, inf))
                    constraints.append((y + d_diff, 0, inf))
                    constraints.append((1 - d_t - y, 0, inf))
        return constraints

    def __solve_hq_subproblem(self, H, Q, f, constraints=None, bounds=None):
        # Caching of the results of this function. SWIG object cannot be
        # pickled directly, so instead we just convert everything to a string.
        pickle_key = (str(H), str(Q), str(f), str(constraints), str(bounds))

        if self.pumpingstation_cache_hq_subproblem:
            try:
                return self.__hq_subproblem_cache[pickle_key]
            except KeyError:
                pass

        # Discharge of pump is always positive
        if bounds is None:
            bounds = {Q.name(): (0.0, np.inf), H.name(): (-np.inf, np.inf)}

        if constraints is None:
            g, lbg, ubg = list(starmap(vertcat, ([], [], [])))
        else:
            g, lbg, ubg = list(starmap(vertcat, list(zip(*constraints))))

        # State vector
        X = [H, Q]

        additional_vars = set(symvar(vertcat(g, f)))
        additional_vars = additional_vars.difference(X)
        X.extend(additional_vars)

        # Bounds
        lbx = vertcat(
            *(bounds[x.name()][0] for x in X[:2]), *([-inf] * len(additional_vars))
        )
        ubx = vertcat(
            *(bounds[x.name()][1] for x in X[:2]), *([inf] * len(additional_vars))
        )

        X = vertcat(*X)

        nlp = {"f": f, "g": g, "x": X}

        # Use the same solver and solver settings as for the overall
        # optimization problem
        options = self.solver_options().copy()

        my_solver = options.pop("solver")

        # Delete unused entries
        del options["optimized_num_dir"]
        options.pop("expand", None)  # Everything is SX already

        casadi_solver = options.pop("casadi_solver")
        if isinstance(casadi_solver, str):
            casadi_solver = getattr(ca, casadi_solver)

        # Remove ipopt and bonmin defaults if they are not used
        if my_solver != "ipopt":
            options.pop("ipopt", None)
        if my_solver != "bonmin":
            options.pop("bonmin", None)

        options.pop("export_lp", None)
        solver = casadi_solver("nlp", my_solver, nlp, options)

        results = solver(x0=np.zeros(X.shape), lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)

        # Sanity check that the solve actually succeeded
        success, _ = self.solver_success(solver.stats(), True)
        if not success:
            raise Exception("Solve of HQ-subproblem failed.")

        objective_value = float(results["f"])
        solver_output = np.array(results["x"])[:, 0]

        self.__hq_subproblem_cache[pickle_key] = objective_value, solver_output[:2]

        return objective_value, solver_output[:2]

    def __solve_working_area_subproblem(
        self, working_area, working_area_direction, head_range, H, Q, f, pump_status=1
    ):
        constraints = self._psmixin_working_area_constraints(
            working_area, working_area_direction, head_range, H, Q, pump_status
        )

        return self.__solve_hq_subproblem(H, Q, f, constraints)

    def __solve_power_subproblem(
        self,
        working_area,
        working_area_direction,
        head_range,
        H,
        Q,
        power_functions,
        pump_status=1,
    ):
        constraints = self._psmixin_working_area_constraints(
            working_area, working_area_direction, head_range, H, Q, pump_status
        )

        power_expr = ca.SX(-np.inf)

        X = SX.sym("X")
        for power_function in power_functions:
            constraints.append((X - power_function, 0.0, inf))
            power_expr = ca.fmax(power_expr, power_function)

        obj_power, (h, q) = self.__solve_hq_subproblem(H, Q, X, constraints)

        if obj_power < 0.0:
            raise Exception("Power inside working area cannot be negative")

        # Check if the inequality resolved to (approximately) an equality
        power_eq = ca.Function("power", [H, Q], [power_expr]).expand()(h, q)

        if not np.allclose(
            power_eq,
            obj_power,
            rtol=self._pumpingstation_ineq_relative_error,
            atol=self._pumpingstation_absolute_error
            * min(abs(power_eq), abs(obj_power)),
        ):
            raise Exception("H-Q subproblem failed to solve power subproblem")

        return obj_power, (h, q)

    def _psmixin_working_area_constraints(
        self, working_area, working_area_direction, head_range, head, discharge, status
    ):
        constraints = []

        for poly, direction in zip(working_area, working_area_direction):
            # When the pump is off, we increasing the working area to
            # include all possible points on the H-axis. We calculate
            # the minimum needed offset to accomplish this.
            offset_min_h = 0.0
            offset_max_h = 0.0

            constr_f = 0.0

            for i in range(poly.shape[0]):
                offset_min_h += head_range[0] ** i * poly[i, 0]
                offset_max_h += head_range[1] ** i * poly[i, 0]

                for j in range(poly.shape[1]):
                    constr_f += poly[i, j] * head**i * discharge**j

            # TODO: Maybe compensate with more than exactly what is
            # needed, e.g. with 130% of the calculated offset. For
            # example, when we have a straight vertical line saying Q
            # > 0.2 m3/s, we might not want to shift the line to
            # exactly Q = 0 m3/s when the pump is off. To make it
            # easier for the mixed integer optimizer to find a good
            # solution, we probably want Q = 0 m3/s to already be an
            # acceptable solution for something like status = 0.2.

            # Apply the working area changes, but only if it increases
            # the working area size. We do not want to shrink it, even
            # if we hypothetically could, as we would rather keep the
            # constraints constant in that case.
            if np.sign(direction) != np.sign(offset_min_h) and np.sign(
                direction
            ) != np.sign(offset_max_h):
                # We are violating both the lowest H-value as well as
                # the highest H-value when the pump is off. We have to
                # compensate only for the largest difference.
                max_offset = np.sign(offset_min_h) * np.max(
                    np.abs([offset_min_h, offset_max_h])
                )
                constr_f -= (1 - status) * max_offset
            elif np.sign(direction) != np.sign(offset_min_h):
                constr_f -= (1 - status) * offset_min_h
            elif np.sign(direction) != np.sign(offset_max_h):
                constr_f -= (1 - status) * offset_max_h

            if direction == -1:
                constraints.append((constr_f, -inf, 0.0))
            elif direction == 1:
                constraints.append((constr_f, 0.0, inf))
            else:
                raise Exception(
                    "Working area polynomial needs a direction of 1 or -1, but got {}".format(
                        direction
                    )
                )

        return constraints

    def path_constraints(self, ensemble_member):
        constraints = super().path_constraints(ensemble_member)

        constant_inputs = self.constant_inputs(ensemble_member)

        for ps in self.pumping_stations():
            sum_q_max = 0
            station_power = 0.0
            station_m = None
            station_Ms = None
            for p in ps.pumps():
                status_sym, _, _, power_sym = self.__pump_status_pairs[p.symbol]
                discharge_sym = p.symbol.replace(".", "_") + "_Q"

                status = self.state(status_sym)

                continuous = 0.0
                if p.semi_continuous is not None:
                    if p.semi_continuous not in constant_inputs:
                        raise TypeError(
                            "Semi-continuous series of pump '{}' not found in constant inputs".format(
                                p.symbol
                            )
                        )

                    continuous_values = constant_inputs[p.semi_continuous].values

                    if not np.all((continuous_values == 0) | (continuous_values == 1)):
                        raise ValueError(
                            "Semi-continuous series of pump '{}' should consist of only 0 and 1".format(
                                p.symbol
                            )
                        )
                    if not np.array_equal(
                        constant_inputs[p.semi_continuous].times, self.times()
                    ):
                        raise ValueError(
                            "Semi-continuous series of pump '{}' should use optimization times".format(
                                p.symbol
                            )
                        )

                    # We need a symbol to easily multiply Timeseries with path variables
                    continuous = self.variable(p.semi_continuous)

                hr = self._psmixin_head_range[p.head_option][ps.symbol]

                constraints.extend(
                    self._psmixin_working_area_constraints(
                        p.working_area,
                        p.working_area_direction,
                        hr,
                        p.head(),
                        p.discharge(),
                        (1 - continuous) * status,
                    )
                )

                # Power calculation which we need for optimization/minimization and constraints.
                coeffs = p.power_coefficients
                # TODO: Is it better if we simplify here directly (e.g. x^0 = 1, x^1 = x)
                powers = self.__power_functions(p.head(), p.discharge(), coeffs)

                m, Ms = self.__pump_power_range_on[power_sym]
                assert len(powers) == len(Ms)
                station_power += self.state(power_sym)

                if len(powers) == 1 and powers[0].is_constant():
                    # Power is constant when pump is on. Can use an equality constraint.
                    constraints.append(
                        (self.state(power_sym) - powers[0] * status, 0.0, 0.0)
                    )
                    station_m += float(powers[0])
                    station_Ms += float(powers[0])
                else:
                    constraints.append((self.state(power_sym) - m * status, 0.0, inf))
                    constraints.append(
                        (self.state(power_sym) - max(Ms) * status, -inf, 0.0)
                    )
                    if station_m is None:
                        station_m = m
                    else:
                        # m is the minimum station power when any pump is on
                        if isinstance(station_m, float):
                            station_m = min(m, station_m)
                        elif isinstance(m, float):
                            station_m = min(m, station_m)
                        else:
                            station_m = np.minimum(m, station_m)
                    if station_Ms is None:
                        station_Ms = Ms
                    else:
                        if isinstance(station_Ms, float):
                            station_Ms = Ms + station_Ms
                        else:
                            station_Ms += Ms

                    for power, M in zip(powers, Ms):
                        # NOTE: Inequality constraint for power, as an equality constraint would have to be affine
                        constraints.append(
                            (
                                self.state(power_sym) - (power - M * (1 - status)),
                                0.0,
                                inf,
                            )
                        )

                # Pump needs to always have a positive discharge
                constraints.append((p.discharge(), 0.0, inf))

                # Pump needs to have zero discharge when off (when not using continuous approach)
                _, q_max = self._psmixin_pump_discharge_bounds[discharge_sym]

                constraints.append((p.discharge() - (status * q_max), -inf, 0.0))
                sum_q_max += q_max
            status_sym, _, _ = self._PumpingStationMixin__station_status_pairs[
                ps.symbol
            ]

            status = self.state(status_sym)
            q_max = sum_q_max
            M = ps.n_pumps
            constraints.append((station_power - station_m * status, 0.0, inf))
            constraints.append((station_power - max(station_Ms) * status, -inf, 0.0))
            constraints.append((ps.discharge() - (status * q_max), -inf, 0.0))
            # For each pumping station, we use the switching matrix only when
            # all the pumps have the same energy price
            if len({p.energy_price_symbol for p in ps.pumps()}) == 1:
                # To handle pump switching constraints easily, we make a
                # vector of the status symbols of all pumps.
                try:
                    ii = [y[0] for y in self.PumpsCardinality].index(self.priority)
                    matrix_indicator = self.PumpsCardinality[ii][1]
                except Exception:
                    try:
                        matrix_indicator = self.PumpsCardinality[0][1]
                    except Exception:
                        matrix_indicator = True

                if matrix_indicator:
                    switch_matrix = ps.pump_switching_matrix
                    if np.all(switch_matrix == 0.0):
                        logger.info("Pumps are independent")
                    else:
                        logger.info("Pumps are dependent according to switching matrix")
                else:
                    switch_matrix = np.zeros_like(ps.pump_switching_matrix)
                    logger.info("Pumps are independent")
                switch_constraints = ps.pump_switching_constraints
                pump_status_vector = np.empty(len(ps.pumps()), dtype=object)
                continuous_variables = []
                for i, p in enumerate(ps.pumps()):
                    if p.semi_continuous is not None:
                        continuous_variables.append(self.variable(p.semi_continuous))
                    else:
                        continuous_variables.append(0.0)

                    pump_status_vector[i] = self.state(
                        self.__pump_status_pairs[p.symbol][0]
                    )

                for i in range(switch_matrix.shape[0]):
                    if any(switch_matrix[i, :] != 0.0):
                        # A pump can only be on when it is allowed to be according to the pump switching matrix
                        # If any pump is continuous, we disable the constraint.
                        any_pump_continuous = ca.mmax(ca.vertcat(*continuous_variables))
                        constraints.append(
                            (
                                (
                                    sum(
                                        np.multiply(
                                            switch_matrix[i, :], pump_status_vector
                                        )
                                    )
                                    * (1 - any_pump_continuous)
                                    + switch_constraints[i, 0] * any_pump_continuous
                                ),
                                switch_constraints[i, 0],
                                switch_constraints[i, 1],
                            )
                        )

            for r in ps.resistances():
                C = r.C
                if C > 0.0:
                    constraints.append(
                        (r.head_loss() - C * r.discharge() ** 2, 0.0, inf)
                    )

                    # To force the head loss to go to zero, we need an upper bound as well.
                    _, max_head_loss = self._psmixin_head_range[0][ps.symbol]
                    q_max_dh = (max_head_loss / C) ** 0.5
                    constraints.append(
                        (
                            max_head_loss / q_max_dh * r.discharge() - r.head_loss(),
                            0.0,
                            inf,
                        )
                    )
                elif C == 0.0:
                    # Force the head loss to zero in case of zero resistance
                    constraints.append((r.head_loss(), 0.0, 0.0))
                else:
                    # Resistance cannot have a negative value
                    raise Exception(
                        'Resistance has a negative value of "{}"'.format(r.C)
                    )

        return constraints

    @property
    def path_variables(self):
        variables = super().path_variables.copy()
        variables.extend(self.__pumping_station_mx_path_variables)
        return variables

    def bounds(self):
        bounds = super().bounds()
        bounds.update(self.__pump_status_bounds)
        bounds.update(self.__pump_power_bounds)
        bounds.update(self._psmixin_pump_discharge_bounds)
        return bounds

    def variable_is_discrete(self, variable):
        if variable in self.__pump_discrete_symbols:
            return True
        if variable in self.__station_discrete_symbols:
            return True
        else:
            return super().variable_is_discrete(variable)

    def priority_completed(self, priority):
        super().priority_completed(priority)

        if not self.pumping_stations():
            return

        times = self.times()
        results = self.extract_results()

        # Round the status. We use this a few times later on to check whether
        # a pump is on, but also to force small errors to zero.
        for ps in self.pumping_stations():
            for p in ps.pumps():
                status_rounded = np.abs(np.around(results[p.symbol + "__status"]))
                self.__additional_results[p.symbol + "_status"] = status_rounded
            if ps.status_history is not None:
                status_rounded = np.ceil(np.around(results[ps.symbol + "__status"], 2))
                self.__additional_results[ps.symbol + "_status"] = status_rounded

        # Extract the pump head and pump discharge from the results.
        # NOTE: We do this here instead of in post(), because in post() we do
        # not have a good solver_output when a priority fails.
        path_expressions = []

        for ps in self.pumping_stations():
            for p in ps.pumps():
                path_expressions.append(p.head())
                path_expressions.append(p.discharge())

        expression = self.map_path_expression(vertcat(*path_expressions), 0)
        f = Function("f", [self.solver_input], [expression])
        evaluated_path_expressions = f(self.solver_output)
        evaluated_path_expressions = np.array(evaluated_path_expressions)

        # Append pump head and discharge to the results dictionary and
        # time series export.
        idx = 0
        for ps in self.pumping_stations():
            for p in ps.pumps():
                head_ts = Timeseries(times, evaluated_path_expressions[:, idx])
                head_key = "{}_{}".format(p.symbol, "head")

                self.__additional_output_variables.add(head_key)
                self.__additional_results[head_key] = head_ts.values
                idx += 1

                status_realised = self.__additional_results[p.symbol + "_status"]

                discharge_ts = Timeseries(
                    times, evaluated_path_expressions[:, idx] * status_realised
                )
                discharge_key = "{}_{}".format(p.symbol, "discharge")

                self.__additional_output_variables.add(discharge_key)
                self.__additional_results[discharge_key] = discharge_ts.values
                idx += 1

        # We estimate the nominal value to be used in the pump energy/cost
        # minimization using results of earlier priorities. Note that we
        # assume that start-up/shut-down costs are small compared to actual
        # running costs, and can therefore be excluded from the nominal
        # calculation.
        tsteps = np.diff(times)

        minimization_nominals = defaultdict(int)  # Default of zero
        # Make sure the dynamic nominal is not too small, e.g. when the
        # pumps happen to all be off
        minimal_nominals = {}
        minimal_nominals[_MinimizePumpGoalType.ENERGY] = min(
            m for m, M in self.__pump_power_range_on.values()
        ) * min(tsteps)
        minimal_nominals[_MinimizePumpGoalType.COST] = np.inf  # To be calculated

        for ps in self.pumping_stations():
            for p in ps.pumps():
                # The power symbol is typically free (with a big-M upper
                # bound). We therefore calculate the  power using the head and
                # discharge results. Note that this is also why we cannot just
                # evaluate a _MinimizePumpGoal's function() on the solver
                # output.
                head_realised = results[p.symbol + ".dH"]
                discharge_realised = results[p.symbol + ".Q"]
                status_realised = np.abs(np.around(results[p.symbol + "__status"]))
                powers_calculated = self.__power_functions(
                    head_realised,
                    discharge_realised,
                    p.power_coefficients,
                    status_realised,
                )

                power_calculated = np.amax(powers_calculated, axis=0)
                results[p.symbol + "_power"] = power_calculated * status_realised

                # Energy
                minimization_nominals[_MinimizePumpGoalType.ENERGY] += power_calculated[
                    1:
                ]

                # Cost
                ts = self.get_timeseries(p.energy_price_symbol)
                prices = self.interpolate(self.times(), ts.times, ts.values)
                minimization_nominals[_MinimizePumpGoalType.COST] += (
                    power_calculated[1:] * prices[1:]
                )
                min_power, max_power = self.__pump_power_range_on[p.symbol + "__power"]

                assert min_power >= 0.0
                if min_power == 0.0:
                    min_power = 0.2 * max_power

                inds = prices[1:] > 0.0
                try:
                    min_cost = min(tsteps[inds] * prices[1:][inds]) * min_power
                    assert min_cost > 0.0
                except Exception:
                    min_cost = 0

                minimal_nominals[_MinimizePumpGoalType.COST] = min(
                    min_cost, minimal_nominals[_MinimizePumpGoalType.COST]
                )

        for k, v in minimization_nominals.items():
            integral_nominal = sum(v * tsteps)
            integral_nominal = max(integral_nominal, minimal_nominals[k])
            self._psmixin_pump_minimization_nominal[k].append(
                (priority, integral_nominal)
            )

        # Store the objective value the minimization goals evaluated to, so
        # they can be checked for proper scaling in post().
        assert self.ensemble_size == 1
        goals = [
            g
            for g in self.goals()
            if g.priority == priority and isinstance(g, _MinimizePumpGoal)
        ]
        for g in goals:
            obj_val = Function("tmp", [self.solver_input], [g.function(self, 0)])(
                self.solver_output
            )
            self.__pump_minimization_function_values.append((g, float(obj_val)))

    def post(self):
        if not self.pumping_stations():
            super().post()
            return

        times = self.times()
        results = self.extract_results()

        # We want to check that the function nominal used for the minimization
        # goal is appropriate. If this is not the case, a warning is raised
        # providing a suitable suggestion for the function nominal.
        any_pump_on = False
        any_costs_made = False

        for ps in self.pumping_stations():
            for p in ps.pumps():
                ts = self.get_timeseries(p.energy_price_symbol)
                price = self.interpolate(times[1:], ts.times, ts.values)

                pump_on = results[p.symbol + "_status"][1:] == 1.0
                price_positive = price > 0.0

                costs_made = pump_on & price_positive

                any_pump_on |= np.any(pump_on)
                any_costs_made |= np.any(costs_made)

        if any_pump_on:
            for g, obj_val in self.__pump_minimization_function_values:
                # print('hi')
                if isinstance(g, MinimizePumpCostGoal) and not any_costs_made:
                    continue
                elif (obj_val > 100 or obj_val < 0.01) and obj_val != 0.0:
                    cur_func_nom = (
                        g._dynamic_nominal
                        if g.use_dynamic_nominal
                        else g.function_nominal
                    )
                    # Reduce to one significant digit
                    unscaled_obj_val = cur_func_nom * obj_val
                    new_func_nom = np.round(
                        unscaled_obj_val,
                        -int(np.floor(np.log10(np.abs(unscaled_obj_val)))),
                    )

                    logger.warning(
                        "Solution may be unstable. Use {} instead of {} as a function nominal of {}".format(
                            new_func_nom, cur_func_nom, g
                        )
                    )

        # TODO: If we put the calculated pump head and discharge in the
        # extract_results() dictionary, should we maybe move the calculation
        # of these time series to priority_completed() in case that exists (or
        # invalidate cached results with every optimize() call)? It might be
        # useful for debugging of intermediate results.

        # Check if the inequality constraints of pump power and resistance
        # head loss have been succesfully minimized to equality.
        results = self.extract_results()

        for ps in self.pumping_stations():
            for p in ps.pumps():
                head_realised = results[p.symbol + "_head"][1:]
                discharge_realised = results[p.symbol + "_discharge"][1:]
                power_realised = results[p.symbol + "__power"][1:]
                power_calculated = results[p.symbol + "_power"][1:]

                if np.any(isinstance(g, MinimizePumpCostGoal) for g in self.goals()):
                    # Minimized for cost. We will skip any time steps where
                    # the price is zero in the check below.
                    ts = self.get_timeseries(p.energy_price_symbol)
                    price = self.interpolate(times[1:], ts.times, ts.values)
                    inds = price > 0.0

                    # price_key = "{}_{}".format(p.symbol, "price")
                    # self.__additional_output_variables.add(price_key)
                    # self.__additional_results[price_key] = price
                else:
                    inds = np.full(power_realised.shape, True)

                max_power = self.__pump_power_bounds[p.symbol + "__power"][1]

                if not np.allclose(
                    power_calculated[inds],
                    power_realised[inds],
                    rtol=self._pumpingstation_ineq_relative_error,
                    atol=self._pumpingstation_absolute_error * max_power,
                ):
                    logger.error(
                        'Relative/absolute power error exceedence in pump "{}"'.format(
                            p.symbol
                        )
                    )

            for r in ps.resistances():
                C = r.C
                head_loss_realised = results[r.symbol + ".dH"][1:]
                discharge_realised = results[r.symbol + ".HQUp.Q"][1:]

                head_loss_target = C * discharge_realised**2

                if not np.allclose(
                    head_loss_target,
                    head_loss_realised,
                    rtol=self._pumpingstation_ineq_relative_error,
                    atol=self._pumpingstation_absolute_error,
                ):
                    logger.error(
                        'Relative/absolute head loss error exceedence in resistance "{}"'.format(
                            r.symbol
                        )
                    )

            # Append pump speed to results and output timeseries
            for ps in self.pumping_stations():
                for p in ps.pumps():
                    coeffs = p.speed_coefficients

                    # Only calculate and output pump speed if non-zero
                    # coefficients are specified
                    if np.all(coeffs == 0):
                        continue

                    head_realised = results[p.symbol + "_head"]
                    discharge_realised = results[p.symbol + "_discharge"]
                    status_realised = results[p.symbol + "_status"]

                    speed = discharge_realised * 0.0

                    for i in range(coeffs.shape[0]):
                        for j in range(coeffs.shape[1]):
                            speed += (
                                coeffs[i, j]
                                * head_realised**i
                                * discharge_realised**j
                                * status_realised
                            )

                    speed_key = "{}_{}".format(p.symbol, "speed")

                    self.__additional_output_variables.add(speed_key)
                    self.__additional_results[speed_key] = speed

            # Append pump power and status to output timeseries
            for ps in self.pumping_stations():
                for p in ps.pumps():
                    power_calculated = results[p.symbol + "_power"]
                    status_realised = results[p.symbol + "_status"]

                    power_key = "{}_{}".format(p.symbol, "power")
                    self.__additional_output_variables.add(power_key)

                    status_key = "{}_{}".format(p.symbol, "status")
                    self.__additional_output_variables.add(status_key)

        # NOTE: If we call super() first, adding output time series with
        # set_time series has no effect, as e.g. PIMIxin/CSVMixin have already
        # written their export file. That is why we do it at the end instead.
        super().post()

    def extract_results(self, *args, **kwargs):
        results = super().extract_results(*args, **kwargs)
        results.update(self.__additional_results)
        return results

    @property
    def output_variables(self):
        variables = super().output_variables.copy()
        variables.extend([ca.MX.sym(v) for v in self.__additional_output_variables])
        return variables


def plot_operating_points(
    optimization_problem,
    output_folder=None,
    plot_expanded_working_area=True,
    plot_specific_energy=False,
    include_prices=False,
):
    """
    Plot the working area of each pump with its operating points.
    """
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    plots = {}

    for ps in optimization_problem.pumping_stations():
        for p in ps.pumps():
            f = plt.figure(figsize=(8, 6))

            # For the head range, we take the extremes of the head over the
            # pump encountered during optimization, and the maximum head
            # inside the working area.
            hr = optimization_problem._psmixin_head_range[p.head_option][ps.symbol]
            hr = [float(x) for x in hr]  # Convert DMatrix to float

            head_sym = p.symbol + "_head"
            if plot_expanded_working_area:
                hrange_wa = (
                    optimization_problem._psmixin_pump_extended_working_area_head_range[
                        head_sym
                    ]
                )
            else:
                hrange_wa = optimization_problem._psmixin_pump_working_area_head_range[
                    head_sym
                ]
            hrange_wa = [float(x) for x in hrange_wa]  # Convert DMatrix to float

            hrange = [min(hr[0], hrange_wa[0]), max(hr[1], hrange_wa[1])]

            discharge_sym = p.symbol.replace(".", "_") + "_Q"
            qrange = optimization_problem._psmixin_pump_discharge_bounds[discharge_sym]
            qrange = [float(x) for x in qrange]  # Convert DMatrix to float

            # For the lines, use a little bit wider range for both H and Q
            extra_space = 0.25 * (qrange[1] - qrange[0])
            qs_range = (qrange[0] - extra_space, qrange[1] + extra_space)

            extra_space = 0.25 * (hrange[1] - hrange[0])
            hs_range = (hrange[0] - extra_space, hrange[1] + extra_space)

            qs = np.linspace(*qs_range)
            hs = np.linspace(*hs_range)[:, None]

            # For the x and y limits we use slightly less extra space. This is
            # to make sure that the contour lines go all the way to the edge
            # of our plots.
            extra_space = 0.1 * (qrange[1] - qrange[0])
            qplot_range = (qrange[0] - extra_space, qrange[1] + extra_space)

            extra_space = 0.1 * (hrange[1] - hrange[0])
            hplot_range = (hrange[0] - extra_space, hrange[1] + extra_space)

            plt.xlim(*qplot_range)
            plt.ylim(*hplot_range)

            # Plot lines for the horizontal and vertical axes
            plt.axhline(0, color="black", zorder=1)
            plt.axvline(0, color="black", zorder=1)

            wa = p.working_area
            wa_dir = p.working_area_direction

            wa_lines = []

            inner_points = qs * hs * 0.0

            # Plot the working area
            for w in range(len(wa)):
                constraints = optimization_problem._psmixin_working_area_constraints(
                    wa[w : w + 1], wa_dir[w : w + 1], hr, hs, qs, 1
                )

                C = plt.contour(
                    qs, hs.ravel(), constraints[0][0], [0], colors="b", zorder=2
                )

                inner_points += ((constraints[0][0] * wa_dir[w]) > 0).astype(int)
                if len(C.allsegs[0]) > 0:
                    wa_lines.append([tuple(x) for x in C.allsegs[0][0]])
                else:
                    wa_lines.append([])

                if plot_expanded_working_area:
                    constraints = (
                        optimization_problem._psmixin_working_area_constraints(
                            wa[w : w + 1], wa_dir[w : w + 1], hr, hs, qs, 0
                        )
                    )

                    plt.contour(
                        qs,
                        hs.ravel(),
                        constraints[0][0],
                        [0],
                        colors="g",
                        linestyles="dashed",
                        zorder=2,
                    )

            plt.plot(
                [0, 0], list(hr), "yo", ms=6, mec="k", label="Head range", zorder=3
            )

            results = optimization_problem.extract_results()

            # Check if we found any point inside the working area, so that we
            # can color it. We typically will not have found such a point for
            # constant speed pumps (or close to constant speed pumps), in which
            # case we skip the filling. If for some other reason we cannot find
            # and enclosing polyline, we just skip the color fill.
            h_inds, q_inds = np.where(inner_points == len(wa))
            if h_inds.size > 0 and q_inds.size > 0:
                try:
                    point = (qs[q_inds[0]], hs[h_inds[0]])
                    wa_segments = enclosing_segments(point, wa_lines)
                    x, y = list(zip(*(s[0] for s in wa_segments)))
                    poly = plt.fill_between(x, y, alpha=0.25, color="none")
                    verts = np.vstack([k.vertices for k in poly.get_paths()])

                    # Add the specific energy as a gradient to the working area of the pump
                    if plot_specific_energy:
                        # Determine coordinates of the working area of the pump to plot the specific energy
                        minh = 999.0
                        maxh = -999.0
                        minq = 999.0
                        maxq = -999.0
                        rows_points, cols_points = np.shape(inner_points)
                        for row in range(rows_points):
                            for col in range(cols_points):
                                if inner_points[row, col] == len(wa):
                                    if qs[col] < minq:
                                        minq = qs[col]
                                    if hs[row] < minh:
                                        minh = hs[row]
                                    if qs[col] > maxq:
                                        maxq = qs[col]
                                    if hs[row] > maxh:
                                        maxh = hs[row]

                        qplot_range = (minq, maxq)
                        hplot_range = (minh, maxh)

                        step_number = 50
                        h_step_length = (hplot_range[1] - hplot_range[0]) / step_number
                        q_step_length = (qplot_range[1] - qplot_range[0]) / step_number
                        h_step_length = h_step_length.item()
                        grid_h, grid_q = np.mgrid[
                            hplot_range[0].item() : hplot_range[1].item()
                            + h_step_length : h_step_length,
                            qplot_range[0].item() : qplot_range[1].item()
                            + q_step_length : q_step_length,
                        ]

                        # Compute specific energy using the defined coeffs in the working area of the pump
                        coeffs = p.power_coefficients
                        spece_array = optimization_problem.spec_energy_functions(
                            grid_h, grid_q, coeffs
                        )
                        # Compute maximum for given set of coeffs
                        spece_array_max = max(spece_array, key=methodcaller("tolist"))

                        # Plot the specific energy in the working area of the pump
                        gradient = plt.imshow(
                            spece_array_max,
                            cmap="cool",
                            aspect="auto",
                            origin="lower",
                            extent=[
                                verts[1 : len(x), 0].min(),
                                verts[1 : len(x), 0].max(),
                                verts[1 : len(y), 1].min(),
                                verts[1 : len(y), 1].max(),
                            ],
                        )
                        gradient.set_clip_path(
                            poly.get_paths()[0], transform=plt.gca().transData
                        )
                        # Add the colorbar to the plot
                        cb = plt.colorbar(
                            gradient,
                            shrink=0.9,
                            label=r"Spec. energy  [$\mathdefault{kWh/1000m^3}$]",
                        )
                        cb.outline.set_color("black")
                        plt.clim(0, math.ceil(spece_array_max.max()))

                except DeadEndError:
                    pass

            # Plot the operating points of the pump
            if include_prices:
                # Add the cost (energy price) as a color gradient to the operating points.
                colors = results[p.energy_price_symbol][1:]
                sc = plt.scatter(
                    results[discharge_sym][1:],
                    results[head_sym][1:],
                    zorder=4,
                    s=100,
                    cmap="YlOrRd",
                    c=colors,
                    marker="x",
                    label="Operating points",
                )
                cb = plt.colorbar(sc, shrink=0.9, label=r"Energy price [€/kWh]")
                cb.outline.set_color("black")
                plt.clim(0, math.ceil(max(colors)))
            else:
                plt.plot(
                    results[discharge_sym][1:],
                    results[head_sym][1:],
                    "rx",
                    markeredgecolor="black",
                    markeredgewidth=2,
                    label="Operating points",
                    zorder=4,
                )

            # Manually add legend entries for the working area(s), because
            # contour plots do not handle that automatically
            handles, _ = plt.gca().get_legend_handles_labels()
            handles.append(mlines.Line2D([], [], color="b", label="Working area"))
            if plot_expanded_working_area:
                handles.append(
                    mlines.Line2D(
                        [], [], color="g", linestyle="--", label="Extended working area"
                    )
                )

            plt.legend(handles=handles)
            plt.xlabel(r"Discharge [$\mathdefault{m^3\!/s}$]")
            plt.ylabel(r"Head [$\mathdefault{m}$]")
            plt.grid(True)
            f.tight_layout()

            if output_folder is not None:
                fname = "QHP_{}.png".format(p.symbol.replace(".", "_"))
                fname = os.path.join(output_folder, fname)
                plt.savefig(fname, bbox_inches="tight", pad_inches=0.1)

            plots[p.symbol] = f

    return plots
