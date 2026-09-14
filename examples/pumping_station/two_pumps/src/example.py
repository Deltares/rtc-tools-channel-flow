import os
from datetime import timedelta

from rtctools.optimization.collocated_integrated_optimization_problem import CollocatedIntegratedOptimizationProblem
from rtctools.optimization.goal_programming_mixin import GoalProgrammingMixin, StateGoal
from rtctools.optimization.linearized_order_goal_programming_mixin import LinearizedOrderGoalProgrammingMixin
from rtctools.optimization.modelica_mixin import ModelicaMixin
from rtctools.optimization.pi_mixin import PIMixin
from rtctools.util import run_optimization_problem

from rtctools_hydraulic_structures.pumping_station_mixin import \
    MinimizePumpCostGoal, PumpStatusGoal, PumpingStation, PumpingStationMixin, plot_operating_points


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

        super().__init__(optimization_problem)


class Example(PumpingStationMixin,
              LinearizedOrderGoalProgrammingMixin, GoalProgrammingMixin,
              PIMixin, ModelicaMixin, CollocatedIntegratedOptimizationProblem):
    """
    An example showing the basic usage of the PumpingStationMixin. It consists of two goals:
    1. Keep water level in the acceptable range.
    2. Minimize power usage for doing so.
    """

    # Set the target minimum and maximum water levels.
    wl_min, wl_max = (-0.5, 0)

    pumpingstation_history_constraints = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.__output_folder = kwargs['output_folder']  # So we can write our pictures to it

        # Here we define a list of pumping stations, each consisting of a list
        # of pumps. In our case, there is only one pumping station containing
        # a single pump.
        self.__pumping_stations = [PumpingStation(self, 'pumpingstation1',
                                                  pump_symbols=['pumpingstation1.pump1',
                                                                'pumpingstation1.pump2'],
                                                  energy_price_symbols=['energy_price_P1',
                                                                        'energy_price_P2'],
                                                  status_history='{pump}_status_hist')]

    def pumping_stations(self):
        # This is the method that we must implement. It has to return a list of
        # PumpingStation objects, which we already initialized in the __init__
        # function. So here we just return that list.
        return self.__pumping_stations

    def goals(self):
        goals = super().goals()
        goals.append(MinimizePumpCostGoal(self))

        for ps in self.pumping_stations():
            for p in ps.pumps():
                g = PumpStatusGoal(self, p)
                g.priority = -1
                goals.append(g)

        return goals

    def path_goals(self):
        goals = super().path_goals()
        goals.append(WaterLevelRangeGoal(self))
        return goals

    def solver_options(self):
        options = super().solver_options()
        options['solver'] = 'cbc'
        options['casadi_solver'] = 'qpsol'
        return options

    def goal_programming_options(self):
        options = super().goal_programming_options()
        options['keep_soft_constraints'] = True
        options['constraint_relaxation'] = 1E-5
        return options

    def post(self):
        super().post()

        results = self.extract_results()

        # TODO: Currently we use hardcoded references to pump1. It would be
        # prettier if we could generalize this so we can handle an arbitrary
        # number of pumps. It would also be prettier to replace hardcoded
        # references to e.g. pumpingstation1.pump1__power with something like
        # pumpingstation1.pump.power(), if at all possible.

        # Calculate the total amount of energy used. Note that QHP fit was
        # made to power in W, and that our timestep is 1 hour.
        power_pump_1 = results['pumpingstation1.pump1__power'][1:]
        power_pump_2 = results['pumpingstation1.pump2__power'][1:]
        powers = power_pump_1 + power_pump_2

        total_power = sum(powers)/1000
        print("Total power = {} kWh".format(total_power))

        # Get index of first time step
        i_start = list(self.get_timeseries('energy_price_P1').times).index(self.initial_time)

        # Calculate total energy cost
        price_pump_1 = self.get_timeseries('energy_price_P1').values[i_start+1:]
        price_pump_2 = self.get_timeseries('energy_price_P2').values[i_start+1:]
        cost = sum(power_pump_1 * price_pump_1 + power_pump_2 * price_pump_2) / 1000
        print("Total cost = {} Euro".format(cost))

        # Make plots
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
        import numpy as np

        plt.style.use('ggplot')

        def unite_legends(axes):
            handles, labels = [], []
            for ax in axes:
                tmp = ax.get_legend_handles_labels()
                handles.extend(tmp[0])
                labels.extend(tmp[1])
            return handles, labels

        # Plot #1: Data over time. X-axis is always time.
        f, axarr = plt.subplots(5, sharex=True)
        times = [self.io.reference_datetime + timedelta(seconds=s) for s in self.times()]

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
        axarr[1].plot(times, self.get_timeseries('H_ext', 0).values[i_start:], label='Sea',
                      linewidth=2, color='b')
        ymin, ymax = axarr[1].get_ylim()
        axarr[1].set_ylim(ymin - 0.5, ymax + 0.5)

        axarr[2].set_ylabel('Energy price\n[EUR/kWh]')
        axarr[2].step(times, self.get_timeseries('energy_price_P1', 0).values[i_start:], label='Energy price',
                      linewidth=2, color='b')
        axarr[2].step(times, self.get_timeseries('energy_price_P2', 0).values[i_start:], label='Energy price',
                      linewidth=2, linestyle='--', color='r')
        ymin, ymax = axarr[2].get_ylim()
        axarr[2].set_ylim(-0.1, ymax + 0.1)

        axarr[3].set_ylabel('Discharge\n[$\\mathdefault{m^3\\!/s}$]')
        axarr[3].step(times, results['pumpingstation1.pump1.Q'], label='Pump 1',
                      linewidth=2, color='b')
        axarr[3].step(times, results['pumpingstation1.pump2.Q'], label='Pump 2',
                      linewidth=2, linestyle='--', color='r')
        axarr[3].plot(times, self.get_timeseries('Q_in', 0).values[i_start:], label='Inflow',
                      linewidth=2, color='g')
        ymin, ymax = axarr[3].get_ylim()

        axarr[3].set_ylim(-0.05 * (ymax - ymin), ymax * 1.1)
        axarr[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

        axarr[4].set_ylabel('Pump speed\n[$\\mathdefault{min^{-1}}$]')
        axarr[4].step(times, results['pumpingstation1.pump1_speed'], label='Pump 1',
                      linewidth=2, color='b')
        axarr[4].step(times, results['pumpingstation1.pump2_speed'], label='Pump 2',
                      linewidth=2, linestyle='--', color='r')
        ymin, ymax = axarr[4].get_ylim()
        axarr[4].set_ylim(-0.05 * (ymax - ymin), ymax * 1.1)
        axarr[4].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

        for ax in axarr:
            ax.set_xlim(min(times), max(times))

        f.autofmt_xdate()

        # Shrink each axis by 20% and put a legend to the right of the axis
        for i in range(len(axarr)):
            box = axarr[i].get_position()
            axarr[i].set_position([box.x0, box.y0, box.width * 0.8, box.height])
            axarr[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

        # Output Plot
        f.set_size_inches(8, 9)
        plt.savefig(os.path.join(self.__output_folder, 'overall_results.png'), bbox_inches='tight', pad_inches=0.1)

        # Plot the working area with the operating points of the pump.
        plot_operating_points(self, self.__output_folder)


# Run
run_optimization_problem(Example, base_folder='..')
