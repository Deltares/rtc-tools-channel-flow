import os

from rtctools.optimization.collocated_integrated_optimization_problem import \
    CollocatedIntegratedOptimizationProblem
from rtctools.optimization.goal_programming_mixin import (GoalProgrammingMixin,
                                                          StateGoal)
from rtctools.optimization.modelica_mixin import ModelicaMixin
from rtctools.optimization.pi_mixin import PIMixin
from rtctools.util import run_optimization_problem

from rtctools_channel_flow.pumping_station_mixin import \
    PumpingStation, PumpingStationMixin, plot_operating_points


class WaterLevelRangeGoal(StateGoal):
    """
    Goal that tries to keep the water level minum and maximum water level,
    the values of which are read from the optimization problem.
    """

    state = 'storage.HQ.H'

    priority = 1

    def __init__(self, optimization_problem):
        self.target_min = optimization_problem.wl_min
        self.target_max = optimization_problem.wl_max

        _range = self.target_max - self.target_min
        self.function_range = (self.target_min - _range, self.target_max + _range)

        super(WaterLevelRangeGoal, self).__init__(optimization_problem)


class Example(PumpingStationMixin, GoalProgrammingMixin, PIMixin, ModelicaMixin,
              CollocatedIntegratedOptimizationProblem):
    """
    An example showing the basic usage of the PumpingStationMixin. It consists of
    two goals:
    1. Keep water level in the acceptable range.
    2. Minimize power usage for doing so.
    """

    # Set the target minimum and maximum water levels.
    wl_min, wl_max = (-0.5, 0)

    def __init__(self, *args, **kwargs):
        super(Example, self).__init__(*args, **kwargs)

        self.__output_folder = kwargs[
            'output_folder']  # So we can write our pictures to it

        # Here we define a list of pumping stations, each consisting of a list
        # of pumps. In our case, there is only one pumping station containing
        # a single pump.
        self.__pumping_stations = [PumpingStation(self, 'pumpingstation1',
                                                  pump_symbols=[
                                                      'pumpingstation1.pump1'])]

    def pumping_stations(self):
        # This is the method that we must implement. It has to return a list of
        # PumpingStation objects, which we already initialized in the __init__
        # function. So here we just return that list.
        return self.__pumping_stations

    def path_goals(self):
        goals = super(Example, self).path_goals()
        goals.append(WaterLevelRangeGoal(self))
        return goals

    def solver_options(self):
        options = super(Example, self).solver_options()
        return options

    def post(self):
        super(Example, self).post()

        results = self.extract_results()

        # TODO: Currently we use hardcoded references to pump1. It would be
        # prettier if we could generalize this so we can handle an arbitrary
        # number of pumps. It would also be prettier to replace hardcoded
        # references to e.g. pumpingstation1.pump1__power with something like
        # pumpingstation1.pump.power(), if at all possible.

        # Calculate the total amount of energy used. Note that QHP fit was
        # made to power in W, and that our timestep is 1 hour.
        powers = results['pumpingstation1.pump1__power'][1:]
        total_power = sum(powers) / 1000
        print("Total power = {} kWh".format(total_power))

        # Make plots
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
        import numpy as np

        plt.style.use('ggplot')
        '''
        def unite_legends(axes):
            h, l = [], []
            for ax in axes:
                tmp = ax.get_legend_handles_labels()
                h.extend(tmp[0])
                l.extend(tmp[1])
            return h, l
        '''
        # Plot #1: Data over time. X-axis is always time.
        f, axarr = plt.subplots(4, sharex=True)
        # TODO: Do not use private API of PIMixin
        times = self.times()

        axarr[0].set_ylabel('Water level\n[m]')
        axarr[0].plot(times, results['storage_level'], label='Polder',
                      linewidth=2, color='b')
        axarr[0].plot(times, self.wl_max * np.ones_like(times), label='Polder Max',
                      linewidth=2, color='r', linestyle='--')
        axarr[0].plot(times, self.wl_min * np.ones_like(times), label='Polder Min',
                      linewidth=2, color='g', linestyle='--')
        ymin, ymax = axarr[0].get_ylim()
        axarr[0].set_ylim(ymin - 0.1, ymax + 0.1)

        axarr[1].set_ylabel('Water level\n[m]')
        axarr[1].plot(times, self.get_timeseries('H_ext', 0).values, label='Sea',
                      linewidth=2, color='b')
        ymin, ymax = axarr[1].get_ylim()
        axarr[1].set_ylim(ymin - 0.5, ymax + 0.5)

        axarr[2].set_ylabel('Energy price\n[EUR/kWh]')
        axarr[2].step(times, self.get_timeseries('energy_price', 0).values,
                      label='Energy price',
                      linewidth=2, color='b')
        ymin, ymax = axarr[2].get_ylim()
        axarr[2].set_ylim(-0.1, ymax + 0.1)

        axarr[3].set_ylabel('Discharge [m^3/s]')
        axarr[3].step(times, results['pumpingstation1.pump1.Q'], label='Pump',
                      linewidth=2, color='b')
        axarr[3].plot(times, self.get_timeseries('Q_in', 0).values, label='Inflow',
                      linewidth=2, color='g')
        ymin, ymax = axarr[3].get_ylim()

        axarr[3].set_ylim(-0.05 * (ymax - ymin), ymax * 1.1)
        axarr[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        f.autofmt_xdate()

        # Shrink each axis by 20% and put a legend to the right of the axis
        for i in range(len(axarr)):
            box = axarr[i].get_position()
            axarr[i].set_position([box.x0, box.y0, box.width * 0.8, box.height])
            axarr[i].legend(loc='center left', bbox_to_anchor=(1, 0.5),
                            frameon=False)

        # Output Plot
        f.set_size_inches(8, 9)
        plt.savefig(os.path.join(self.__output_folder, 'overall_results.png'),
                    bbox_inches='tight', pad_inches=0.1)

        # Plot the working area with the operating points of the pump.
        plot_operating_points(self, self._output_folder)


# Run
run_optimization_problem(Example, base_folder='..')
