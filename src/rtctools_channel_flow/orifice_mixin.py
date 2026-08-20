import logging
import os
from abc import abstractmethod
from collections import OrderedDict

from casadi import MX

import numpy as np
from numpy import inf

from rtctools.optimization.goal_programming_mixin import GoalProgrammingMixin
from rtctools.optimization.optimization_problem import OptimizationProblem
from rtctools.optimization.timeseries import Timeseries

from .util import _ObjectParameterWrapper

logger = logging.getLogger("rtc-hydraulic-structures")
logging.basicConfig(level=logging.INFO)


class Orifice(_ObjectParameterWrapper):
    """
    Python Orifice object as an interface to the Orifice object in the model.
    """

    def __init__(self, optimization_problem, name):
        super().__init__(optimization_problem)
        self.optimization_problem = optimization_problem
        self.symbol = name

    def discharge(self):
        """
        Get the state corresponding to the weir discharge.

        :returns: `MX` expression of the weir discharge.
        """
        return self.optimization_problem.state(self.symbol + '.Q')

    def _head_up(self):
        """
        Get the state corresponding to the pump head. This depends on the
        ``head_option`` that was specified by the user.

        :returns: `MX` expression of the pump head.
        """
        return self.optimization_problem.state(self.symbol + '.HQUp.H')

    def _head_down(self):
        return self.optimization_problem.state(self.symbol + '.HQDown.H')

    def head(self):
        return self._head_down() - self._head_up()

    @property
    def dh_max(self):
        return self.dH_max

    def _calc_q(self, dh):
        assert np.all(dh >= 0.0)
        return self.discharge_coefficient * self.area * (2.0 * 9.81 * dh)**0.5

    @property
    def q_max(self):
        return self._calc_q(self.dh_max)


class OrificeMixin(OptimizationProblem):
    """
    Adds handling of Orifice objects in your model to your optimization problem.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        assert isinstance(self, GoalProgrammingMixin), "Can only use OrificeMixin when using GoalProgrammingMixin"

    @abstractmethod
    def orifices(self):
        """
        User problem returns list of :class:`Orifice` objects.

        :returns: A list of orifices.
        """
        return []

    def pre(self):
        super().pre()
        self._orifice_discrete_symbols = []
        self._orifice_mx_path_variables = []
        self._orifice_status_bounds = OrderedDict()
        self._orifice_status_symbols = OrderedDict()

        for o in self.orifices():
            status_sym = '{}__status'.format(o.symbol)
            self._orifice_discrete_symbols.append(status_sym)
            self._orifice_status_bounds[status_sym] = (0, 1)
            self._orifice_mx_path_variables.append(MX.sym(status_sym))
            self._orifice_status_symbols[o.symbol] = status_sym

    def path_constraints(self, ensemble_member):
        constraints = super().path_constraints(ensemble_member)
        bounds = self.bounds()

        for o in self.orifices():
            status_sym = self._orifice_status_symbols[o.symbol]
            status = self.state(status_sym)

            M = o.dh_max
            q_max = bounds[o.symbol + ".Q"][1]

            # Force downhill symbol to be 1 when head is negative (H_down < H_up),
            # and 0 otherwise.
            constraints.append((o.head() - (1 - status) * M, -inf, 0.0))
            constraints.append((o.head() + status * M, 0.0, inf))

            # Force discharge to be zero if head is not downhill
            constraints.append((o.discharge() + (1 - status) * q_max, 0.0, q_max))

            g = 9.81
            constraints.append(
                (((o.discharge() / (o.area * o.discharge_coefficient)) ** 2) / (2 * g) +
                 o.head() - M * (1 - status),
                 -inf, 0.0))

        return constraints

    @property
    def path_variables(self):
        variables = super().path_variables
        variables.extend(self._orifice_mx_path_variables)
        return variables

    def bounds(self):
        bounds = super().bounds()
        bounds.update(self._orifice_status_bounds)

        for o in self.orifices():
            q_sym_name = o.symbol + ".Q"
            m, M = bounds[q_sym_name]
            assert m == 0.0, "Orifice minimum discharge must be zero"
            M = min(M, o.q_max)
            bounds[q_sym_name] = (m, M)

        return bounds

    def variable_is_discrete(self, variable):
        if variable in self._orifice_discrete_symbols:
            return True
        else:
            return super().variable_is_discrete(variable)

    def priority_completed(self, priority):
        """
        Calculate additional timeseries relevant to the orifice:

        - Percentage open to realize the calculated discharge. Assumes linear
          relationship between opening and discharge.
        """
        super().priority_completed(priority)

        times = self.times()
        results = self.extract_results()

        # Round the status.
        for o in self.orifices():
            dh = results[o.symbol + ".HQUp.H"] - results[o.symbol + ".HQDown.H"]
            dh_clipped = np.clip(dh, 0.0, o.dh_max)

            max_qs = o._calc_q(dh_clipped)
            actual_qs = results[o.symbol + ".Q"]

            inds = (max_qs == 0.0)
            max_qs[inds] = 0.1 * o.q_max  # A reasonable value so we can still when something is wrong
            open_fraction = actual_qs / max_qs

            self.set_timeseries(o.symbol + "_fraction_open", Timeseries(times, open_fraction))


def plot_operating_points(optimization_problem, output_folder):
    """
    Plot the operating area of the orifice with its operating points.
    """
    import matplotlib.pyplot as plt

    results = optimization_problem.extract_results()

    for o in optimization_problem.orifices():
        plt.clf()

        # For the head range, we take the extremes of the head over the
        # orifice encountered during optimization, and the maximum head
        # difference specified by the user.
        head_loss = results[o.symbol + ".HQUp.H"] - results[o.symbol + ".HQDown.H"]
        hrange = (min(head_loss), o.dh_max)
        qrange = (0.0, o.q_max)

        # For the lines, use a little bit wider range for both H
        extra_space = 0.25 * (hrange[1] - hrange[0])
        hs_range = (hrange[0] - extra_space, hrange[1] + extra_space)

        hs = np.linspace(*hs_range)[:, None]

        # For the x and y limits we use slightly less extra space. For the H
        # we use no extra space on the upper limit, to avoid confusion.
        extra_space = 0.1 * (qrange[1] - qrange[0])
        qplot_range = (qrange[0] - extra_space, qrange[1] + extra_space)

        extra_space = 0.1 * (hrange[1] - hrange[0])
        hplot_range = (hrange[0] - extra_space, hrange[1])

        plt.xlim(*hplot_range)
        plt.ylim(*qplot_range)

        # Plot lines for the horizontal and vertical axes
        plt.axhline(0, color='black', zorder=1)
        plt.axvline(0, color='black', zorder=1)

        # Plot the maximum discharge at each head, if the orifice is fully open
        hs_clipped = np.clip(hs, 0.0, o.dh_max)
        plt.plot(hs, o._calc_q(hs_clipped), 'b')

        # Plot the operating points
        discharge_sym = o.symbol + ".Q"
        plt.plot(head_loss[1:], results[discharge_sym][1:], 'r+',
                 markeredgewidth=2, label="Operating points")

        plt.xlabel(r'$\Delta H$ [$\mathdefault{m}$]')
        plt.ylabel(r'Discharge [$\mathdefault{m^3\!/s}$]')

        f = plt.gcf()
        f.set_size_inches(8, 6)
        f.tight_layout()

        plt.grid(True)

        fname = '{}_operating_points.png'.format(o.symbol.replace('.', '_'))
        fname = os.path.join(output_folder, fname)
        plt.savefig(fname, bbox_inches='tight', pad_inches=0.1)
