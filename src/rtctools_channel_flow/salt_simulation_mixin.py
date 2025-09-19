import logging
logger = logging.getLogger("rtctools")
import pandas as pd
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import os

def area_coef(storage_areas, number):
    coef = 0.0
    for idx in range(len(storage_areas)):
        coef +=storage_areas[number] / storage_areas[idx]
    return coef

def set_mid_q_wl_closing(self, net_q_storages, net_q_system):
    storage_areas = [self.parameters()[x + '.A'] for x in self.active_storage_names]

    for idx, storage_name in enumerate(self.active_storage_names):
        if idx == 0:
            q_net_previous = 0
        if idx == len(self.active_storage_names) - 1:
            pass
        else:
            q_net_box_needed = net_q_system / area_coef(storage_areas, idx)

            q_mid_set = q_net_previous + net_q_storages[idx] - q_net_box_needed
            q_net_previous = q_mid_set
            if q_mid_set < -0.001:
                logger.debug(
                    "{} has negative downstream discharge, check the laterals.".format(
                        storage_name))
                #raise Exception
            if self.upstream_open_boundary:
                connector_name  = self.connector_names[idx + 1]
            else:
                connector_name = self.connector_names[idx]
            self.set_var(connector_name + '_middle_discharge', q_mid_set)
            print('Q mid set for {} at {}: {}'.format(storage_name, self.get_current_time(), q_mid_set))

def set_mid_q_q_closing(self, net_q_storages, net_q_system):

    for idx, storage_name in enumerate(self.active_storage_names):
        if idx == 0:
            q_net_previous = 0
        if idx == len(self.active_storage_names) - 1:
            q_mid_set = q_net_previous + net_q_storages[idx]
            self.set_var(storage_name + '_qforcing_advective', -q_mid_set)
            print('Q down set for {} at {}: {}'.format(storage_name, self.get_current_time(), q_mid_set))
        else:
            q_mid_set = q_net_previous + net_q_storages[idx]
            q_net_previous = q_mid_set
            if q_mid_set < -0.001:
                logger.debug(
                    "{} has negative downstream discharge, check the laterals.".format(
                        storage_name))
                raise Exception
            if self.upstream_open_boundary:
                connector_name  = self.connector_names[idx + 1]
            else:
                connector_name = self.connector_names[idx]
            self.set_var(connector_name + '_middle_discharge', q_mid_set)
            print('Q mid set for {} at {}: {}'.format(storage_name, self.get_current_time(), q_mid_set))


class SaltSimulationMixin():
    """
    Class for plotting results.
    """

    def initialize(self):

        if self.use_set_vales_file:
            #Read and overwrite timeseries_import.csv for ease of use
            logger.info("Timeseries import is overwritten with values from set_values")
            input_data = pd.read_csv('..\\input\\set_values.csv')
            concentrations_data = pd.read_csv('..\\input\\initial_state.csv')

            upstream_storage_name = self.active_storage_names[0]
            downsream_storage_name = self.active_storage_names[-1]

            upstream_q_forcing_in = input_data.iloc[0, input_data.columns.get_loc(upstream_storage_name + '_qforcing_in')]
            downstream_q_forcing_in = input_data.iloc[0, input_data.columns.get_loc(downsream_storage_name + '_qforcing_in')]
            downstream_q_forcing_out = input_data.iloc[0, input_data.columns.get_loc(downsream_storage_name + '_qforcing_out')]

            if not self.upstream_open_boundary:
                logger.info(
                    "No upstream open boundary, sea concentration is set to the value in initial state")
            sea_concnetration = concentrations_data.iloc[
                    0, concentrations_data.columns.get_loc(self.storage_names[0] + '.C')]
            if not self.downstream_open_boundary:
                logger.info("No downstream open boundary, sea concentration is set to the value in initial state")
            lake_concnetration = concentrations_data.iloc[0, concentrations_data.columns.get_loc(self.storage_names[-1] + '.C')]

            storage1_mforcing_in = sea_concnetration * upstream_q_forcing_in
            storage3_mforcing_in = lake_concnetration * downstream_q_forcing_in

            ts_import = pd.read_csv('..\\input\\timeseries_import_orig.csv')

            ts_import[upstream_storage_name + '_qforcing_in'] = ts_import[upstream_storage_name + '_qforcing_in'].map(lambda x: upstream_q_forcing_in)
            ts_import[downsream_storage_name + '_qforcing_in'] = ts_import[downsream_storage_name + '_qforcing_in'].map(lambda x: downstream_q_forcing_in)
            ts_import[upstream_storage_name + '_mforcing_in'] = ts_import[upstream_storage_name + '_mforcing_in'].map(lambda x: storage1_mforcing_in)
            ts_import[downsream_storage_name + '_mforcing_in'] = ts_import[downsream_storage_name + '_mforcing_in'].map(lambda x: storage3_mforcing_in)
            ts_import[downsream_storage_name + '_qforcing_out'] = ts_import[downsream_storage_name + '_qforcing_out'].map(lambda x: downstream_q_forcing_out)


            ts_import.to_csv('..\\input\\timeseries_import.csv', index=False)

        for name in self.connector_names:
            self.set_var(name + '_middle_discharge', 0.0)

        super().initialize()

    def read(self):
        input_database = pd.read_csv(self._input_folder + "/timeseries_import.csv", parse_dates=True, index_col=[0])
        input_database.to_csv(self._input_folder + "/timeseries_import.csv", date_format='%Y-%m-%d %H:%M:%S')
        print("Changed timeseries input format.")
        super().read()

    def update(self, dt):
        # Get the time step
        if dt < 0:
            dt = self.get_time_step()
        time_step = self.get_current_time() / dt

        QForcing = np.zeros(len(self.active_storage_names))
        for idx, name in enumerate(self.active_storage_names):
            QForcing[idx] = self.io.get_timeseries(name + '_qforcing_in')[1][int(time_step+1)] + self.io.get_timeseries(name + '_qforcing_out')[1][int(time_step+1)]

        net_q_system = np.sum(self.ZSF_Q) + np.sum(QForcing)
        net_q_storages = QForcing.copy()
        net_q_storages[0] += self.ZSF_Q[0]
        net_q_storages[-1] += self.ZSF_Q[1]

        if self.water_level_closing:
            set_mid_q_wl_closing(self, net_q_storages, net_q_system)
        else:
            set_mid_q_q_closing(self, net_q_storages, net_q_system)

        super().update(dt)


    def extra_equations(self):
        equations = super().extra_equations()
        equation_list = []
        variables = self.get_variables()

        for name in self.connector_names:
            var_in = name + '.HQUp.M[1]'
            var_out = name + '.HQDown.M[1]'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx + var_out_mx)

            var_in = name + '.HQDown.Q'
            var_out = name + '_middle_discharge'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx + var_out_mx)

            var_in = name + '.HQUp.Q'
            var_out = name + '_middle_discharge'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx - var_out_mx)

        for idx, name in enumerate(self.active_storage_names):
            var_in = name + '.QForcing[1]'
            var_out = name + '_qforcing_in'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx - var_out_mx)

            var_in = name + '.QForcing[2]'
            var_out = name + '_qforcing_out'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx - var_out_mx)

            var_in = name + '.MForcing[1]'
            var_out = name + '_mforcing_in'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            equation_list.append(var_in_mx - var_out_mx)

            var_in = name + '.MForcing[2]'
            var_out = name + '_qforcing_out'
            var_out_mul = name + '.HQDown.C[1]'
            var_in_mx = variables[var_in]
            var_out_mx = variables[var_out]
            var_out_mul_mx = variables[var_out_mul]
            equation_list.append(var_in_mx - var_out_mx * var_out_mul_mx)

            if idx == 0:
               var_in = name + '.QForcing[4]'
               var_out = name + '_qforcing_advective'
               var_in_mx = variables[var_in]
               var_out_mx = variables[var_out]
               equation_list.append(var_in_mx - var_out_mx)

               var_in = name + '.MForcing[4]'
               var_out = name + '_mforcing_advective'
               var_in_mx = variables[var_in]
               var_out_mx = variables[var_out]
               equation_list.append(var_in_mx - var_out_mx)

            if idx == len(self.active_storage_names)-1:
                var_in = name + '.QForcing[4]'
                var_out = name + '_qforcing_advective'
                var_in_mx = variables[var_in]
                var_out_mx = variables[var_out]
                equation_list.append(var_in_mx - var_out_mx)

                var_in = name + '.MForcing[4]'
                var_out = name + '_qforcing_advective'
                var_out_mul = name + '.HQDown.C[1]'
                var_in_mx = variables[var_in]
                var_out_mx = variables[var_out]
                var_out_mul_mx = variables[var_out_mul]
                equation_list.append(var_in_mx - var_out_mx * var_out_mul_mx)

            if idx==0 or idx == len(self.active_storage_names) - 1:
               #ZSF
               
               var_in = name + '.QForcing[3]'
               var_out = name + '_qforcing_ZSF'
               var_in_mx = variables[var_in]
               var_out_mx = variables[var_out]
               equation_list.append(var_in_mx - var_out_mx)
               logger.debug("Appending equation: in: {} and out: {}".format(var_in, var_out))
               
               var_in = name + '.MForcing[3]'
               var_out = name + '_mforcing_ZSF'
               var_in_mx = variables[var_in]
               var_out_mx = variables[var_out]
               equation_list.append(var_in_mx - var_out_mx)
               logger.debug("Appending equation: in: {} and out: {}".format(var_in, var_out))

        equations += equation_list
        return equations

    # Plotting
    def get_output_variables(self):

        variables = super().get_output_variables().copy()

        for idx, name in enumerate(self.active_storage_names):
           variables.extend([ name + ".MForcing[2]"])
           if idx == len(self.active_storage_names) - 1 or idx == 0:
               variables.extend([name + "_qforcing_ZSF"])
               variables.extend([name + "_mforcing_ZSF"])
               variables.extend([name + "_qforcing_advective"])
               variables.extend([name + ".MForcing[4]"])

        for name in self.connector_names:
            variables.extend([name + ".flux_q1_s1"])
            variables.extend([name + ".HQUp.Q"])

        for name in self.storage_names:
            variables.extend([name + ".V"])

        return variables



    def post(self):

        super().post()

        results = self.extract_results()
        np.set_printoptions(suppress=True)

        color_list=['darkviolet','orange','green','magenta', 'yellow','deepskyblue','black','forestgreen','brown','pink']
        f, axarr = plt.subplots(10, sharex=True)
        plt.subplots_adjust(left=0.1, bottom=0.1,
                    top=0.95, wspace=0.4, hspace=0.85)
        times= self.times()/3600

        #Plot 1
        if self.upstream_open_boundary:
           axarr[0].plot(times, results['concentration_' + self.storage_names[0]], label= self.storage_names[0] + 'C',
                      linewidth=2, color='blue')

        for idx, storage_name in enumerate(self.active_storage_names):
           axarr[0].plot(times, results['concentration_' + storage_name], label= storage_name + 'C',
                      linewidth=2, color=color_list[idx])

        if self.downstream_open_boundary:
           axarr[0].plot(times, results['concentration_' + self.storage_names[-1]], label= storage_name + 'C',
                      linewidth=2, color='red')

        axarr[0].set_ylabel('Concentration\n[kg/m3]')
        ymin, ymax = axarr[0].get_ylim()
        axarr[0].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[0].legend()
        axarr[0].set_title("Concentrations", fontsize = 10)

        #Plot 2
        for idx, storage_name in enumerate(self.active_storage_names):
           axarr[1].plot(times, results[storage_name + '.V'] / self.parameters()[storage_name + '.A'],
                      linewidth=2,  color=color_list[idx], linestyle='-')#, label='H_' +storage_name)
        if self.upstream_open_boundary:
           axarr[1].plot(times, results[self.storage_names[0] + '.V'] / self.parameters()[self.storage_names[0] + '.A'], label='meer',
                      linewidth=2, color='b')
        if self.downstream_open_boundary:
           axarr[1].plot(times, results[self.storage_names[-1] + '.V'] /  self.parameters()[self.storage_names[-1] + '.A'], label='zee',
                      linewidth=2, color='r', linestyle='--')
        axarr[1].set_ylabel('Water level\n[m]')
        ymin, ymax = axarr[1].get_ylim()
        axarr[1].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[1].legend()
        axarr[1].set_title("Water levels", fontsize = 10)

        # Plot 3
        for idx, connector_name in enumerate(self.connector_names):
            if self.upstream_open_boundary and idx == 0:
                axarr[2].plot(times, results[connector_name + '.flux_q1_s1'],
                              label='q_uit_' + connector_name, linewidth=2, color='b')
            #elif self.downstream_open_boundary and idx==len(self.connector_names)-1:
            #    axarr[2].plot(times, results[connector_name + '.flux_q1_s1'],
            #                  label='q_uit_' + connector_name, linewidth=2, color=color_list[idx])
            else:
                axarr[2].plot(times, results[connector_name + '.flux_q1_s1'],
                              label='q_uit_' + connector_name, linewidth=2,
                              color=color_list[idx])

        axarr[2].set_ylabel('Discharge\n[m3/s]')
        ymin, ymax = axarr[2].get_ylim()
        axarr[2].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[2].legend()
        axarr[2].set_title("Dispersive transport discharge", fontsize = 10)

        # Plot 4
        #Todo: the timeseires is referred as a given name, it does not find it with the . name
        for idx, connector_name in enumerate(self.connector_names):
            if self.upstream_open_boundary and idx == 0:
                axarr[3].plot(times, results[connector_name + '_M_Up'],
                              label='M_uit_' + connector_name, linewidth=2, color='b')
            else:
                axarr[3].plot(times, results[connector_name + '_M_Up'],
                              label='M_uit_' + connector_name, linewidth=2,
                              color=color_list[idx])

        axarr[3].set_ylabel('Mass flux\n[kg/s]')
        ymin, ymax = axarr[3].get_ylim()
        axarr[3].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[3].legend()
        axarr[3].set_title("Total (dispersive and advective) mass flux", fontsize = 10)

        #Plot 5
        for idx, connector_name in enumerate(self.connector_names):
            if self.upstream_open_boundary and idx == 0:
                pass
            elif self.downstream_open_boundary and idx==len(self.connector_names)-1:
                pass
            else:
                axarr[4].plot(times, results[connector_name + '.HQUp.Q'],
                              label='q_adv,mid_from_' + connector_name, linewidth=2,
                              color=color_list[idx])

        for idx, storage_name in enumerate(self.active_storage_names):
           if abs(sum(self.io.get_timeseries(storage_name+'_qforcing_in')[1]))>0.0001:
               axarr[4].plot(times, self.io.get_timeseries(storage_name+'_qforcing_in')[1], label='q_in_' + storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           if abs(sum(self.io.get_timeseries(storage_name + '_qforcing_out')[1])) > 0.0001:
               axarr[4].plot(times, -self.io.get_timeseries(storage_name+'_qforcing_out')[1], label='q_out_' + storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           plt.yticks(np.arange(-100, 100, 50))
        axarr[4].set_ylabel('Discharge\n[m3/s]')
        axarr[4].set_title("Advective transport discharge", fontsize = 10)


        # Plot 6
        # Todo: this can only take values outside the mixin...

        for idx, storage_name in enumerate(self.active_storage_names):
           if abs(sum(self.io.get_timeseries(storage_name+'_qforcing_in')[1]))>0.0001:
               axarr[5].plot(times, self.io.get_timeseries(storage_name+'_qforcing_in')[1], label='q_in_' + storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           if abs(sum(self.io.get_timeseries(storage_name + '_qforcing_out')[1])) > 0.0001:
               axarr[5].plot(times, -self.io.get_timeseries(storage_name + '_qforcing_out')[1], label='q_out_' +storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           plt.yticks(np.arange(-100, 100, 50))



        axarr[5].set_ylabel('Discharge\n[m3/s]')
        ymin, ymax = axarr[5].get_ylim()
        axarr[5].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[5].legend()
        axarr[5].set_title('Extra inflow', fontsize = 10)


        # Plot 7
        for idx, storage_name in enumerate(self.active_storage_names):
           if abs(sum(self.io.get_timeseries(storage_name+'_mforcing_in')[1]))>0.0001:
               axarr[6].plot(times, self.io.get_timeseries(storage_name+'_mforcing_in')[1], label='m_in_' + storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           if abs(sum(results[storage_name + '.MForcing[2]'])) > 0.0001:
               axarr[6].plot(times, -results[storage_name + '.MForcing[2]'], label='m_out_' +storage_name,
                      linewidth=2, linestyle='--', color=color_list[idx])
           #plt.yticks(np.arange(-100, 100, 50))

        axarr[6].set_ylabel('mass flux\n[kg/s]')
        ymin, ymax = axarr[6].get_ylim()
        axarr[6].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[6].legend()
        axarr[6].set_title('Extra flux', fontsize = 10)


        # Plot 7
        #Todo: maybe here the names are not too generic
        axarr[7].plot(times, results[self.active_storage_names[0] + '_qforcing_ZSF'], label='ZSF up',
                      linewidth=2, color='b')
        axarr[7].plot(times, results[self.active_storage_names[-1] + '_qforcing_ZSF'], label='ZSF down',
                      linewidth=2, color='r', linestyle='--')
        axarr[7].plot(times, -results[self.active_storage_names[-1] + '_qforcing_advective'], label='downstream_flushing',
                      linewidth=2, color='g', linestyle='--')
        ymin, ymax = axarr[7].get_ylim()
        axarr[7].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[7].legend()
        axarr[7].set_ylabel('Discharge\n[m3/s]')
        axarr[7].set_title('ZSF and flushing discharge', fontsize=10)

        # Plot 8
        # Todo: maybe here the names are not too generic
        axarr[8].plot(times, results[self.active_storage_names[0] + '_mforcing_ZSF'], label='ZSF up',
                      linewidth=2, color='b')
        axarr[8].plot(times, results[self.active_storage_names[-1] + '_mforcing_ZSF'], label='ZSF down',
                      linewidth=2, color='r', linestyle='--')
        axarr[8].plot(times, -results[self.active_storage_names[-1] + '.MForcing[4]'], label='downstream_flushing',
                      linewidth=2, color='g', linestyle='--')
        axarr[8].set_ylabel('Mass flux\n[kg/s]')
        ymin, ymax = axarr[8].get_ylim()
        axarr[8].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[8].legend()
        axarr[8].set_title('ZSF and flushing flux', fontsize=10)
        plt.sca(axarr[8])
        plt.yticks([0.0, np.round(ymax/2,-1)])

        # Plot 8
        axarr[9].plot(times, self.io.get_timeseries('upstream_discharge')[1], label='up',
                      linewidth=2, color='b')
        axarr[9].plot(times, self.io.get_timeseries('downstream_discharge')[1], label='down',
                      linewidth=2, color='r', linestyle='--')
        axarr[9].set_ylabel('Discharge\n[m3/s]')
        ymin, ymax = axarr[9].get_ylim()
        axarr[9].set_ylim(ymin - 0.1, ymax + 0.1)
        axarr[9].legend()
        axarr[9].set_title('Discharge boundaries', fontsize=10)
        plt.sca(axarr[9])
        plt.yticks([0.0, np.round(ymax/2,-1)])

        axarr[-1].set_xlabel('Time [h]')
        f.autofmt_xdate()
        for ax in axarr:
            ax.set_xlim(min(times), max(times))

        # Shrink each axis by 20% and put a legend to the right of the axis
        for i in range(len(axarr)):
            box = axarr[i].get_position()
            axarr[i].set_position([box.x0, box.y0, box.width * 0.8, box.height])
            axarr[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
        

        # Output Plot
        f.set_size_inches(8, 9)

        plt.savefig(os.path.join(self._output_folder, 'overall_results.png'), bbox_inches='tight', pad_inches=0.1, dpi=300)

        df = pd.read_csv('..\\output\\timeseries_export.csv')
        #small_df = df[['concentration_storage1', 'concentration_storage1', 'concentration_storage3', 'connector_1.HQUp.Q', 'connector_2.HQUp.Q','storage3_qforcing_advective']].copy()
        #small_df=small_df.rename(columns = {'storage3_qforcing_advective':'downstream_sluiting_q'})
        #small_df.to_csv('..\\output\\timeseries_export_short.csv')