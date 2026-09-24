from datetime import timedelta

from rtctools.optimization.collocated_integrated_optimization_problem \
    import CollocatedIntegratedOptimizationProblem
from rtctools.optimization.csv_mixin import CSVMixin
from rtctools.optimization.goal_programming_mixin import GoalProgrammingMixin, StateGoal
from rtctools.optimization.modelica_mixin import ModelicaMixin
from rtctools.util import run_optimization_problem

from rtctools_channel_flow.weir_mixin import Weir, WeirMixin, plot_operating_points

# There are two water level targets, with different priority.
# The water level should stay in the required range during all the simulation


class WLRangeGoal_01(StateGoal):

    def __init__(self, optimization_problem):
        self.state = 'Branch1.H'
        self.priority = 1

        self.target_min = optimization_problem.get_timeseries('h_min_branch1')
        self.target_max = optimization_problem.get_timeseries('h_max_branch1')

        super().__init__(optimization_problem)


class WLRangeGoal_02(StateGoal):

    def __init__(self, optimization_problem):
        self.state = 'Branch2.H'
        self.priority = 2

        self.target_min = optimization_problem.get_timeseries('h_min_branch2')
        self.target_max = optimization_problem.get_timeseries('h_max_branch2')

        super().__init__(optimization_problem)


class WeirExample(WeirMixin, GoalProgrammingMixin, CSVMixin, ModelicaMixin, CollocatedIntegratedOptimizationProblem):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__weirs = [Weir(self, 'weir1')]
        self.__output_folder = kwargs['output_folder']  # So we can write our pictures to it

    def weirs(self):
        return self.__weirs

    def solver_options(self):
        options = super().solver_options()
        solver = options['solver']
        options[solver]['allowable_gap'] = 0.005
        options[solver]['print_level'] = 2
        return options

    def path_goals(self):
        goals = super().path_goals()
        goals.append(WLRangeGoal_01(self))
        goals.append(WLRangeGoal_02(self))
        return goals

    def post(self):
        super().post()
        results = self.extract_results()

        # Make plots
        import matplotlib.dates as mdates
        import matplotlib.pyplot as plt
        import os

        plt.style.use('ggplot')

        def unite_legends(axes):
            handles, labels = [], []
            for ax in axes:
                tmp = ax.get_legend_handles_labels()
                handles.extend(tmp[0])
                labels.extend(tmp[1])
            return handles, labels

        # Plot #1: Data over time. X-axis is always time.
        f, axarr = plt.subplots(4, sharex=True)

        times = [self.io.reference_datetime + timedelta(seconds=s) for s in self.times()]
        weir = self.weirs()[0]

        axarr[0].set_ylabel('Water level\n[m]')
        axarr[0].plot(times, results['branch_1_water_level'], label='Upstream',
                      linewidth=2, color='b')
        axarr[0].plot(times, self.get_timeseries('h_min_branch1').values, label='Upstream Min',
                      linewidth=2, color='r', linestyle='--')
        axarr[0].plot(times, self.get_timeseries('h_max_branch1').values, label='Upstream Max',
                      linewidth=2, color='g', linestyle='--')
        ymin, ymax = axarr[0].get_ylim()
        axarr[0].set_ylim(ymin - 0.1, ymax + 0.1)

        axarr[1].set_ylabel('Water level\n[m]')
        axarr[1].plot(times, results['branch_2_water_level'], label='Downstream',
                      linewidth=2, color='b')
        axarr[1].plot(times, self.get_timeseries('h_max_branch2').values, label='Downstream Max',
                      linewidth=2, color='r', linestyle='--')
        axarr[1].plot(times, self.get_timeseries('h_min_branch2').values, label='Downstream Min',
                      linewidth=2, color='g', linestyle='--')
        ymin, ymax = axarr[1].get_ylim()
        axarr[1].set_ylim(ymin - 0.1, ymax + 0.1)

        axarr[2].set_ylabel('Discharge\n[$\\mathdefault{m^3\\!/s}$]')
        # We need the first point for plotting, but its value does not
        # necessarily make sense as it is not included in the optimization.
        weir_flow = results['WeirFlow1']
        weir_height = results["weir1_height"]
        minimum_water_level_above_weir = (weir_flow-weir.q_nom)/weir.slope + weir.h_nom

        minimum_weir_height = minimum_water_level_above_weir - ((weir_flow/weir.c_weir)**(2.0/3.0))

        minimum_weir_height[0] = minimum_weir_height[1]

        weir_flow = results['WeirFlow1']
        weir_flow[0] = weir_flow[1]

        axarr[2].step(times, weir_flow, label='Weir',
                      linewidth=2, color='b')
        axarr[2].step(times, self.get_timeseries('upstream_q_ext').values, label='Inflow',
                      linewidth=2, color='r', linestyle='--')
        axarr[2].step(times, -1 * self.get_timeseries('downstream_q_ext').values, label='Outflow',
                      linewidth=2, color='g', linestyle='--')
        ymin, ymax = axarr[2].get_ylim()
        axarr[2].set_ylim(-0.1, ymax + 0.1)

        weir_height = results["weir1_height"]
        weir_height[0] = weir_height[1]

        axarr[3].set_ylabel('Weir height\n[m]')
        axarr[3].step(times, weir_height, label='Weir',
                      linewidth=2, color='b')
        ymin, ymax = axarr[3].get_ylim()
        ymargin = 0.1 * (ymax - ymin)
        axarr[3].set_ylim(ymin - ymargin, ymax + ymargin)
        axarr[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        f.autofmt_xdate()

        # Shrink each axis by 20% and put a legend to the right of the axis
        for ax in axarr:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.set_xlim(min(times), max(times))
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

        # Output Plot
        f.set_size_inches(8, 9)
        plt.savefig(os.path.join(self.__output_folder, 'overall_results.png'),
                    bbox_inches='tight', pad_inches=0.1)

        plot_operating_points(self, self.__output_folder, results)


if __name__ == "__main__":
    run_optimization_problem(WeirExample)
