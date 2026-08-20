import logging
import os
from abc import abstractmethod
from collections import OrderedDict
from typing import List

from casadi import MX

from numpy import inf

from rtctools.optimization.optimization_problem import OptimizationProblem

from .util import _ObjectParameterWrapper

logger = logging.getLogger("rtc-hydraulic-structures")
logging.basicConfig(level=logging.INFO)


class Weir(_ObjectParameterWrapper):
    """
    Python Weir object as an interface to the Weir object in the model.

    In the optimization, the weir flow is implemented as constraints. It means
    that the optimization calculated a flow (not weir height!), that is forced
    by the constraints to be a physically possible weir flow.
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

    def _head(self):
        return self.optimization_problem.state(self.symbol + '.HQUp.H')

    @property
    def c_weir(self):
        # coefficient part of the equation
        return (2 * 9.81)**0.5 * 2.0 / 3.0 * self.weir_coef * self.width

    @property
    def q_nom(self):
        # half of the possible dischage
        return (self.q_max + self.q_min) / 2.0

    @property
    def h_nom(self):
        # H corresponding to half of the possible dischage
        return (self.q_nom / self.c_weir)**(2.0 / 3.0) + self.hw_min

    @property
    def slope(self):
        return self.c_weir * 3.0 / 2.0 * (self.h_nom - self.hw_min)**0.5


class WeirMixin(OptimizationProblem):
    """
    Adds handling of Weir objects in your model to your optimization problem.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # We also want to output a bunch of variables that are calculated in
        # post(), which we return in the extract_results() call. We store them
        # in this dictionary.
        self.__additional_results = OrderedDict()

    @abstractmethod
    def weirs(self) -> List[Weir]:
        """
        User problem returns list of :class:`Weir` objects.

        :returns: A list of weirs.
        """
        return []

    def pre(self):
        super().pre()
        self.__weir_discrete_symbols = []
        self.__weir_mx_path_variables = []
        self.__weir_status_bounds = OrderedDict()
        self.__weir_status_pairs = OrderedDict()

        for w in self.weirs():
            # Define symbol names
            status_sym = '{}__status'.format(w.symbol)
            # Store variable names
            self.__weir_discrete_symbols.append(status_sym)
            # Store bounds
            self.__weir_status_bounds[status_sym] = (0, 1)
            # Generate and store MX variables
            status_mx = MX.sym(status_sym)
            self.__weir_mx_path_variables.append(status_mx)
            # Store all symbols together
            self.__weir_status_pairs[w.symbol] = (status_sym)

    def path_constraints(self, ensemble_member):
        constraints = super().path_constraints(ensemble_member)

        for w in self.weirs():
            status_sym = self.__weir_status_pairs[w.symbol]
            status = self.state(status_sym)

            # calculate minimum possible discharge (weir @ max)
            slope_max = w.q_max / ((w.q_max / w.c_weir) ** (2.0/3.0))
            q_min_h = slope_max * w._head() - slope_max * w.hw_max

            # calculate maximum possible discharge (weir @ min), using
            # linearized curve
            q_max_h = w.slope*(w._head() - w.h_nom-1 * (status-1)) + w.q_nom

            # flow should be lower than physical maximum and bigger then zero
            constraints.append((w.discharge() - w.q_max*(status), -inf, 0))
            epsilon = 0.00001
            # Here a small value, epsilon is used, write down the equation
            # status * epsilon < Q, with other words if the weir is "on",
            # this is the minimum flow. If the weir is "off" this flow can still be
            # present. So it should be small...
            constraints.append((-w.discharge() + epsilon*(status), -inf, 0))

            # flow should be lower than max related to water level and weir
            # height
            constraints.append((w.discharge() - q_max_h, -inf, 0))

            # flow should be higher than min related to water level and weir
            # height
            constraints.append((w.discharge() - q_min_h, 0, inf))

            # flow should higher than physical minimum
            constraints.append((-w.discharge(), -inf, -w.q_min))

        return constraints

    @property
    def path_variables(self):
        variables = super().path_variables
        variables.extend(self.__weir_mx_path_variables)
        return variables

    def bounds(self):
        bounds = super().bounds()
        bounds.update(self.__weir_status_bounds)
        return bounds

    def variable_is_discrete(self, variable):
        if variable in self.__weir_discrete_symbols:
            return True
        else:
            return super().variable_is_discrete(variable)

    def post(self):

        results = self.extract_results()

        # Calculating the weir height
        for w in self.weirs():
            weir_wl_up = results[w.symbol + ".HQUp.H"]
            weir_q = results[w.symbol + ".Q"]

            self.__additional_results[w.symbol + "_height"] = weir_wl_up - (weir_q / w.c_weir)**(2.0/3.0)

        # NOTE: If we call super() first, adding output time series with
        # set_time series has no effect, as e.g. PIMIxin/CSVMixin have already
        # written their export file. That is why we do it at the end instead.
        super().post()

    def extract_results(self, *args, **kwargs):
        results = super().extract_results(*args, **kwargs)
        results.update(self.__additional_results)
        return results


def plot_operating_points(optimization_problem, output_folder, results):
    # This is for post processing, for plotting the flow and height of a weir
    import matplotlib.pyplot as plt
    import numpy as np

    for w in optimization_problem.weirs():
        weir_name = w.symbol
        weir_flow_results = results[weir_name + '.Q']
        water_level = results[weir_name + '.HQUp.H']

        # Calculating the working area of the weir
        hmargin = 0.25 * (w.hw_max - w.hw_min)
        hs = np.linspace(w.hw_min - hmargin, w.hw_max + w.q_max/(w.c_weir) + hmargin, 100)
        q_max_h_plot = w.slope * (hs - w.h_nom) + w.q_nom
        q_min_h_plot = w.c_weir * np.fmax(0, (hs - w.hw_max))**1.5
        slope_max = w.q_max / ((w.q_max / w.c_weir) ** (2.0/3.0))
        q_min_h_lin_plot = slope_max * hs - slope_max * w.hw_max
        q_max_th_plot = w.c_weir * np.fmax(0, (hs - w.hw_min))**1.5
        q_min_plot = np.full((100), w.q_min)
        q_max_plot = np.full((100), w.q_max)

        # Plotting the working area of the weir
        plt.clf()
        plt.plot(hs, q_min_plot, 'k-')
        plt.plot(hs, q_max_h_plot, 'k-')
        plt.plot(hs, q_max_plot, 'k-')
        plt.plot(hs, q_min_h_plot, 'g-')
        plt.plot(hs, q_min_h_lin_plot, 'k-')
        plt.plot(hs, q_max_th_plot, 'b-')
        plt.plot(water_level, weir_flow_results, 'r+', markeredgewidth=2)
        plt.ylim(w.q_min - 0.5, w.q_max * 1.2)
        plt.title('Working area of the weir')
        plt.xlabel('Water level [m]')
        plt.ylabel('Flow [$m^3\\,s^{-1}$]')
        plt.xlim(min(hs), max(hs))
        save_name = weir_name + "_working_area" + ".png"
        fname = os.path.join(output_folder, save_name)
        plt.savefig(fname, bbox_inches='tight', pad_inches=0.1)

        if w.q_max > 2 * max(weir_flow_results[1:]):
            logger.warning('The given maximum weir flow ({}) is much higher than the actual maximum flow ({}). '
                           'This might lead to an unncessarily big linearization error.'.format(
                                w.q_max, max(weir_flow_results)))
