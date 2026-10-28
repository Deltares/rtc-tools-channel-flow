import casadi as ca
import numpy as np
import pytest

from rtctools_channel_flow.pumping_station_mixin import PumpingStationMixin

_solve_power_subproblem = PumpingStationMixin._PumpingStationMixin__solve_power_subproblem


class _StubMixin:
    """Just enough of PumpingStationMixin to run its power-subproblem check."""

    _pumpingstation_ineq_relative_error = (
        PumpingStationMixin._pumpingstation_ineq_relative_error
    )
    _pumpingstation_absolute_error = PumpingStationMixin._pumpingstation_absolute_error

    def __init__(self, obj_power, h, q):
        self._result = (obj_power, np.array([h, q]))

    def _psmixin_working_area_constraints(self, *args):
        return []

    def _PumpingStationMixin__solve_hq_subproblem(self, H, Q, f, constraints):
        return self._result


@pytest.mark.parametrize("excess", [0.0, 1e-6], ids=["equal", "power_eq_smaller"])
def test_power_subproblem_check_with_casadi_dm(excess):
    # power_eq is evaluated through a casadi Function, so it starts out as a DM.
    # When it is <= obj_power it ends up in np.allclose's atol, which recursed
    # without bound on casadi >=3.8 (abs() of a DM no longer returns a float).
    H = ca.SX.sym("H")
    Q = ca.SX.sym("Q")
    h, q = 2.0, 3.0
    obj_power = h * q + excess

    stub = _StubMixin(obj_power, h, q)
    result_power, (result_h, result_q) = _solve_power_subproblem(
        stub, None, None, None, H, Q, [H * Q]
    )

    assert result_power == obj_power
    assert (result_h, result_q) == (h, q)


def test_power_subproblem_check_still_rejects_mismatch():
    H = ca.SX.sym("H")
    Q = ca.SX.sym("Q")
    stub = _StubMixin(7.0, 2.0, 3.0)

    with pytest.raises(Exception, match="failed to solve power subproblem"):
        _solve_power_subproblem(stub, None, None, None, H, Q, [H * Q])
