import logging
import numpy as np
from rtctools_channel_flow.calculate_parameters import GetLinearSVVariables

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

    :cvar linearised_sv: True if LinearisedSV branch is used.
                        Default is ``None``.
    :cvar linearised_sv_branches:  List of LinearisedSV branches to set parameters for.
                        Default is ``None``.
    :cvar linearised_sv_use_dynamic_nominals: True if dynamic nominal water levels
                        should be used for LinearisedSV branches.
                        Default is ``False``.
    :cvar linearised_sv_nominal_levels: dict
        A dictionary with nominal water levels for each branch.
        For each branch, nominals can be provided for up and down.
        values can be a fixed value or a timeseries.

        key: branch name,
            key: "H_b_up" or "H_b_down":
            value: nominal water levels. Can be a fixed value or a timeseries.
                These are the values around which the linearization is done.
                Possible types:

                - fixed value: float
                - timeseries name: str, name of the timeseries to use for nominal levels.
                    Note that the value at the start of the optimization run will be
                    used as the nominal level.
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
    """

    linearised_sv = None
    linearised_sv_branches = None
    linearised_sv_use_dynamic_nominals = False
    linearised_sv_nominal_levels = None

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
                        "dynamic nominal water levels for LinearisedSV block."
                    )
                p = self.set_linear_sv_dynamic_nominal(p=p)
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

    def set_linear_sv_dynamic_nominal(self, p):
        """
        Set the dynamic nominal water levels for the block:
        Deltares.ChannelFlow.Hydraulic.Branches.LinearisedSV.
        If the required nominal data is not provided, the
        default nominal levels are used.

        :param p: The parameters of the model.

        :return: Updated parameters with dynamic nominal water levels.
        """
        for channel in self.linearised_sv_branches:
            if channel not in self.linearised_sv_nominal_levels:
                logger.warning(
                    f"Nominal water levels for channel {channel} are not provided "
                    "while linearised_sv_use_dynamic_nominals is True. Cannot set "
                    "dynamic nominal water levels for LinearisedSV block. Using "
                    "default nominal levels."
                )
                continue
            for key in ["H_b_up", "H_b_down", "Q_nominal"]:
                # check if key is present for channel
                if key not in self.linearised_sv_nominal_levels[channel]:
                    continue
                # check if nominal is a float
                if isinstance(self.linearised_sv_nominal_levels[channel][key], float):
                    nominal_level = self.linearised_sv_nominal_levels[channel][key]
                # check if nominal is a timeseries name
                elif isinstance(self.linearised_sv_nominal_levels[channel][key], str):
                    # check if timeseries exists
                    ts_name = self.linearised_sv_nominal_levels[channel][key]
                    if ts_name not in self.io.get_timeseries_names():
                        logger.warning(
                            f"Nominal water level timeseries {ts_name} for channel "
                            f"{channel} and key {key} does not exist. Cannot set "
                            "dynamic nominal water levels for LinearisedSV block."
                        )
                        continue
                    nominal_level = self.get_timeseries(
                        ts_name
                    ).values[self.timeseries_import.forecast_index]
                    # check if timeseries has a value at forecast_index
                    if np.isnan(nominal_level):
                        logger.warning(
                            f"Nominal water level timeseries {ts_name} "
                            f"does not have a value at forecast index "
                            f"{self.timeseries_import.forecast_index}. Cannot set "
                            "dynamic nominal water levels for LinearisedSV block."
                        )
                        continue
                else:
                    logger.warning(
                        f"Nominal water level for channel {channel} and key {key} is "
                        "not a float or a timeseries name. Cannot set dynamic nominal "
                        "water levels for LinearisedSV block."
                    )
                    continue

                if key == "H_b_up":
                    depth = nominal_level - p[channel + ".H_b_up"]
                    # nominal should be positive and non-zero
                    # if value is zero net nominal to 1
                    if depth <= 0:
                        raise ValueError(
                            f"Calculated dynamic nominal water level for channel {channel} is non-positive. "
                            f"Calculated value: {depth}. Provided nominal water level is {nominal_level} "
                            f"and upstream bed level is {p[channel + '.H_b_up']}."
                        )
                    p[channel + ".H_nominal"] = depth
                    msg = (
                        f"Set dynamic nominal {channel + '.H_nominal'} water level for channel {channel} "
                        f"to {depth}."
                    )
                    inform_once(logger, msg)
                elif key == "H_b_down":
                    depth = nominal_level - p[channel + ".H_b_down"]
                    # nominal should be positive and non-zero
                    # if value is zero net nominal to 1
                    if depth <= 0:
                        raise ValueError(
                            f"Calculated dynamic nominal water level for channel {channel} is non-positive. "
                            f"Calculated value: {depth}. Provided nominal water level is {nominal_level} "
                            f"and downstream bed level is {p[channel + '.H_b_down']}."
                        )
                    p[channel + ".H_nominal_down"] = depth
                    msg = (
                        f"Set dynamic nominal {channel + '.H_nominal_down'} water level for channel {channel} "
                        f"to {depth}."
                    )
                    inform_once(logger, msg)
                elif key == "Q_nominal":
                    p[channel + ".Q_nominal"] = nominal_level
                    logger.debug(
                        f"Set dynamic nominal {channel + '.Q_nominal'} flow for channel {channel} "
                        f"to {nominal_level}."
                    )

        return p
