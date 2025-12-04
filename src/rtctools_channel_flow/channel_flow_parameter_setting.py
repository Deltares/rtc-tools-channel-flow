import logging
from rtctools_channel_flow.calculate_parameters import GetLinearSVVariables

logger = logging.getLogger("rtctools")


class ChannelFlowParameterSettingOpimizationMixin:
    """
    Sets parameters for channel flow blocks using in optimization mode.

    Supported blocks:
    Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV

    :cvar LinearisedSV: True if LinearisedSV branch is used.
                        Default is ``None``.
    :cvar sv_branches:  List of LinearisedSV branches to set parameters for.
                        Default is ``None``.
    """

    LinearisedSV = None
    sv_branches = None

    def parameters(self, ensemble_member):
        """
        Set the parameters for the channel flow blocks.
        """
        p = super().parameters(ensemble_member)

        if self.LinearisedSV:
            p = self.set_linear_sv_parameters(p=p)
            if self.use_dynamic_nominals:
                p = self.set_linear_sv_dynamic_nominal(p=p)
        return p

    def set_linear_sv_parameters(self, p):
        """
        Set the parameters for the block:
        Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV.

        :param p: The parameters of the model.
        """
        p["step_size"] = 0.0

        use_semi_implicit = False
        use_convective_acceleration = True
        use_upwind = False

        for channel in self.sv_branches:
            p[channel + ".use_convective_acceleration"] = use_convective_acceleration
            p[channel + ".use_upwind"] = use_upwind
            variableGetter = GetLinearSVVariables(
                n_level_nodes=int(p[channel + ".n_level_nodes"]),
                length=float(p[channel + ".length"]),
                h_b_up=float(p[channel + ".H_b_up"]),
                h_b_down=float(p[channel + ".H_b_down"]),
                q_nominal=float(p[channel + ".Q_nominal"]),
                width=float(p[channel + ".width"]),
                y_nominal=float(p[channel + ".H_nominal"]),
                y_nominal_down=float(p[channel + ".H_nominal_down"]),
                friction_coefficient=float(p[channel + ".friction_coefficient"]),
            )
            variableGetter.getVariables()
            for i in range(len(variableGetter.t0)):
                p[channel + ".T0[" + str(i + 1) + "]"] = variableGetter.t0[i]
            for i in range(len(variableGetter.v0)):
                p[channel + ".V0[" + str(i + 1) + "]"] = variableGetter.v0[i]
            for i in range(len(variableGetter.delta)):
                p[channel + ".Delta[" + str(i + 1) + "]"] = variableGetter.delta[i]
            for i in range(len(variableGetter.gamma)):
                p[channel + ".Gamma[" + str(i + 1) + "]"] = variableGetter.gamma[i]
            for i in range(len(variableGetter.c0)):
                p[channel + ".C0[" + str(i + 1) + "]"] = variableGetter.c0[i]
            logger.debug(
                f"Set Linear SV parameters for channel {channel} for channel flow"
                " block Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV"
            )

        return p

    def set_linear_sv_dynamic_nominal(self, p):
        """
        Set the dynamic nominal water levels for the block:
        Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV.

        :param p: The parameters of the model.
        """
        for channel in self.sv_branches:
            reach_number = int(channel.split("reach_", 1)[1])
            nominal_level = self.get_timeseries("H_reach_" + str(reach_number)).values[
                self.timeseries_import.forecast_index
            ]

            value = nominal_level - p[channel + ".H_b_up"]
            p[channel + ".H_nominal"] = value

            if reach_number < 10:
                nominal_level = self.get_timeseries(
                    "H_reach_" + str(reach_number + 1)
                ).values[self.timeseries_import.forecast_index]
            else:
                nominal_level = self.get_timeseries(
                    "H_reach_" + str(reach_number)
                ).values[self.timeseries_import.forecast_index]
            value = nominal_level - p[channel + ".H_b_down"]
            p[channel + ".H_nominal_down"] = value
        return p
