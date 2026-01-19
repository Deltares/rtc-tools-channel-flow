import logging
from rtctools_channel_flow.calculate_parameters import (
    GetLinearSVVariables,
    GetIDZVariables,
)

from functools import lru_cache

logger = logging.getLogger("rtctools")


@lru_cache(None)
def inform_once(logger, msg):
    logger.info(msg)


class ChannelFlowParameterSettingOpimizationMixin:
    """
    Sets parameters for channel flow blocks using in optimization mode.

    Supported blocks:
    Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV
    Deltares.ChannelFlow.Hydraulic.Branches.IDZ

    :cvar linearised_sv: True if LinearisedSV branch is used.
                        Default is ``None``.
    :cvar linearised_sv_branches:  List of LinearisedSV branches to set parameters for.
                        Default is ``None``.
    :cvar linearised_sv_use_dynamic_nominals: True if dynamic nominal water depths
                        should be used for LinearisedSV branches.
                        Default is ``False``.
    :cvar linearised_sv_nominal_levels: dict
        A dictionary with nominal water depths for each branch.
        For each branch, nominals can be provided for up and down.
        values can be a fixed value or a timeseries.

        key: branch name,
            key: "H_b_up" or "H_b_down":
            value: nominal water depths. Can be a fixed value or a timeseries.
                These are the values around which the linearization is done.
                Possible types:

                - fixed value: float
                - timeseries name: str, name of the timeseries to use for nominal depths.
                    Note that the value at the start of the optimization run will be
                    used as the nominal depth.
            key: "Q_nominal":
                These are the values around which the linearization is done.
                Possible types:

                - fixed value: float
                - timeseries name: str, name of the timeseries to use for nominal flows.
                    Note that the value at the start of the optimization run will be
                    used as the nominal flow.


        example:
            {my_branch_name: {
                "H_b_up": 2.5,
                "H_b_down": "H_nominal_down_timeseries",
                "Q_nominal": "Q_nominal_timeseries_name"
            }}
        :cvar idz: True if IDZ branch is used.
        Default is ``None``.

    :cvar idz_branches: List of IDZ branches to set parameters for.
        Each entry should be the full Modelica path to an
        ``Deltares.ChannelFlow.Hydraulic.Branches.IDZ`` block.
        Default is ``None``.

    :cvar idz_use_dynamic_nominals: True if dynamic nominal water depths
        should be used for IDZ branches.
        Default is ``False``.

    :cvar idz_nominal_levels: dict
        A dictionary with nominal water depths for each IDZ branch.

        For each branch, a nominal water depth can be provided that is
        used as the linearization point for the IDZ formulation.

        key: branch name
            key: "H_nominal":
                Nominal water depth for the IDZ branch.

                Possible types:
                - fixed value: float
                - timeseries name: str, name of the timeseries to use for
                  nominal depth.

                Note that if a timeseries is provided, the value at the
                start of the optimization run will be used as the nominal
                water depth.

        example:
            {
                my_idz_branch_name: {
                    "H_nominal": 2.0
                }
            }
    """

    linearised_sv = None
    linearised_sv_branches = None
    linearised_sv_use_dynamic_nominals = False
    linearised_sv_nominal_levels = None
    idz = None
    idz_branches = None
    idz_use_dynamic_nominals = False

    def parameters(self, ensemble_member):
        """
        Set the parameters for the channel flow blocks.
        """
        p = super().parameters(ensemble_member)

        if self.linearised_sv:
            # check if linearised_sv_branches is set
            if self.linearised_sv_branches is None:
                raise ValueError(
                    "List of linearised_sv_branches is not provided while linearised_sv "
                    "is True. Cannot set parameters for LinearisedSV block."
                )
            p = self.set_linear_sv_parameters(p=p)
            if self.linearised_sv_use_dynamic_nominals:
                # check if linearised_sv_nominal_levels is set
                if self.linearised_sv_nominal_levels is None:
                    raise ValueError(
                        "Dictionary of linearised_sv_nominal_levels is not provided "
                        "while linearised_sv_use_dynamic_nominals is True. Cannot set "
                        "dynamic nominal water depths for LinearisedSV block."
                    )
                p = self.set_linear_sv_dynamic_nominal(p=p)

        if self.idz:
            if self.idz_branches is None:
                raise ValueError(
                    "List of idz is not provided while idz "
                    "is True. Cannot set parameters for IDZ block."
                )
            p = self.set_idz_parameters(p=p)
            if self.idz_use_dynamic_nominals:
                # check if linearised_sv_nominal_levels is set
                if self.idz_nominal_levels is None:
                    raise ValueError(
                        "Dictionary of idz_nominal_levels is not provided "
                        "while idz_use_dynamic_nominals is True. Cannot set "
                        "dynamic nominal water depths for IDZ block."
                    )
                p = self.set_idz_dynamic_nominal(p=p)
        return p

    def set_linear_sv_parameters(self, p):
        """
        Set the parameters for the block:
        Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV.

        :param p: The parameters of the model.
        :return: Updated parameters with LinearisedSV parameters.
        """

        # TODO: expand this method to accept the following parameters:
        # p["step_size"] = 0.0
        # use_semi_implicit = False
        # use_convective_acceleration = True
        # use_upwind = False

        required_params = [
            ".n_level_nodes",
            ".length",
            ".H_b_up",
            ".H_b_down",
            ".Q_nominal",
            ".width",
            ".H_nominal",
            ".H_nominal_down",
            ".friction_coefficient",
        ]

        for channel in self.linearised_sv_branches:
            # check if all required parameters are present
            for param in required_params:
                if channel + param not in p:
                    raise ValueError(
                        f"Parameter {channel + param} is required to set "
                        "LinearisedSV parameters but is missing. "
                        "This parameter can be set via the declaration of "
                        f"the LinearisedSV branch, {channel}, in the model "
                        ".mo file."
                    )
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

    def set_idz_parameters(self, p):
        """
        Set the parameters for the block:
        ``Deltares.ChannelFlow.Hydraulic.Branches.IDZ``.

        This method computes and assigns the internal IDZ parameters
        based on geometric properties, hydraulic characteristics, and
        nominal operating conditions of the channel.

        The following parameters must be defined in the model prior
        to calling this method (typically via the Modelica ``.mo`` file):

        Required parameters per IDZ branch:
            - ``.length`` :
                Channel length [m]
            - ``.H_b`` :
                Bed level [m]
            - ``.Q_nominal`` :
                Nominal discharge [m³/s]
            - ``.width`` :
                Bottom width of the channel [m]
            - ``.H_nominal`` :
                Nominal water depth [m]
            - ``.friction_coefficient`` :
                Hydraulic friction coefficient [-]
            - ``.side_slope`` :
                Side slope of the channel cross section [-]

        Based on these inputs, the following IDZ parameters are computed
        and written to the parameter dictionary:

            - ``p11, p12, p21, p22`` :
                Linearized state-space coefficients of the IDZ model
            - ``Au`` :
                Upstream wetted cross-sectional area [m²]
            - ``Ad`` :
                Downstream wetted cross-sectional area [m²]
            - ``Delay_in_hour`` :
                Hydraulic delay through the channel [h]

        :param p: dict
            Dictionary containing the current model parameters.
        :return: dict
            Updated parameter dictionary including computed IDZ parameters.
        """

        required_params = [
            ".length",
            ".H_b_up",
            ".H_b_down",
            ".Q_nominal",
            ".width",
            ".H_nominal",
            ".friction_coefficient",
            ".side_slope",
        ]

        for channel in self.idz_branches:
            # check if all required parameters are present
            for param in required_params:
                if channel + param not in p:
                    raise ValueError(
                        f"Parameter {channel + param} is required to set "
                        "IDZ parameters but is missing. "
                        "This parameter can be set via the declaration of "
                        f"the IDZ branch, {channel}, in the model "
                        ".mo file."
                    )
            variableGetter = GetIDZVariables(
                length=float(p[channel + ".length"]),
                h_b_up=float(p[channel + ".H_b_up"]),
                h_b_down=float(p[channel + ".H_b_down"]),
                q_nominal=float(p[channel + ".Q_nominal"]),
                width=float(p[channel + ".width"]),
                y_nominal=float(p[channel + ".H_nominal"]),
                friction_coefficient=float(p[channel + ".friction_coefficient"]),
                side_slope=float(p[channel + ".side_slope"]),
            )
            variableGetter.getVariables()
            p[channel + ".p11"] = variableGetter.p11
            p[channel + ".p12"] = variableGetter.p12
            p[channel + ".p21"] = variableGetter.p21
            p[channel + ".p22"] = variableGetter.p22
            p[channel + ".Au"] = variableGetter.Au
            p[channel + ".Ad"] = variableGetter.Ad
            p[channel + ".Delay_in_hour"] = variableGetter.Delay_in_hour
            logger.debug(
                f"Set IDZ parameters for channel {channel} for channel flow"
                " block Deltares.ChannelFlow.Hydraulic.Branches.IDZ"
            )

        return p

    def set_idz_dynamic_nominal(self, p):
        """
        Set the dynamic nominal water depths for the block:
        Deltares.ChannelFlow.Hydraulic.Branches.IDZ.
        If the required nominal data is not provided, the
        default nominal depths are used.

        :param p: The parameters of the model.

        :return: Updated parameters with dynamic nominal water depths.
        """
        raise NotImplementedError("Dynamic nominal idz is not yet implemented.")

        return p
