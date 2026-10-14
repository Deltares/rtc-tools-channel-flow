"""Compile and solve every example, guarding against unit-resolution failures.

The library's Modelica source is otherwise untested: nothing else here parses a
`.mo` file. These tests exist to catch a regression in how `Modelica.Units.SI.*`
types resolve, which is why a missing-class error always fails, even where Bug 3
is expected.

Two known upstream defects stop the examples short of a solution, so only unit
resolution is asserted unconditionally:

- pymoca drops `min`/`max` modifiers on MSL-typed connector variables, so
  bounds arrive as `(-inf, inf)` and rtc-tools rejects the goals
  (rtc-tools#1798 "Bug 3"). Confirmed on 0.11.2; 0.9.3 is unaffected; appears
  fixed on 0.12.0, though that rests on one example (`weir/basic` XPASSes while
  the other five stop at an earlier 0.12.0 failure). The condition below is an
  open-ended `>= 0.10` because the first affected release was not bisected.
  Marked expected, strictly: once CI's pymoca pin is bumped past a release that
  fixes this, the XPASS fails the run and the marker has to go.
- casadi >=3.8 recurses until the stack is exhausted on `np.allclose()` applied
  to a `DM`, which `pumping_station_mixin` does. Unrelated to Modelica; it stops
  three examples before any optimisation runs.

With pymoca 0.9.3 and casadi 3.7.2 all six examples solve, `rc=0`.
"""

import re
import subprocess
import sys
from pathlib import Path

import casadi
import pymoca
import pytest

EXAMPLES = Path(__file__).parent.parent / "examples"

# Bug 3 reaches us as a bare Exception from rtc-tools, so it cannot be told
# apart from any other failure by type. Raise our own to keep the xfail narrow.
class BoundsLost(Exception):
    """pymoca discarded min/max on an MSL-typed variable (rtc-tools#1798)."""


def _pymoca_drops_modifiers():
    major, minor, *_ = pymoca.__version__.split(".")
    return (int(major), int(minor)) >= (0, 10)


PYMOCA_DROPS_MODIFIERS = _pymoca_drops_modifiers()

BOUNDS_LOST = re.compile(
    r"\(-inf, inf\)|No function range specified for goal",
)
MISSING_CLASS = re.compile(
    r"Could not find class|ClassNotFoundError",
)
CASADI_RECURSION_MARKERS = ("RecursionError", "_numpy_array_function_dispatch")


def _example_dirs():
    return sorted(p.parent.parent for p in EXAMPLES.glob("*/*/src/example.py"))


def _example_id(path):
    return f"{path.parent.name}/{path.name}"


@pytest.mark.parametrize("example", _example_dirs(), ids=_example_id)
@pytest.mark.xfail(
    PYMOCA_DROPS_MODIFIERS,
    raises=BoundsLost,
    strict=True,
    reason="pymoca >=0.10 drops min/max on MSL-typed variables (rtc-tools#1798)",
)
def test_example_solves(example):
    # example.py calls run_optimization_problem() at import time and resolves
    # its input/output folders relative to src/, so it has to run as a script.
    result = subprocess.run(
        [sys.executable, "example.py"],
        cwd=example / "src",
        capture_output=True,
        text=True,
        timeout=1800,
    )
    output = result.stdout + result.stderr

    # Never expected, on any pymoca version: this is what the suite is for.
    assert not MISSING_CLASS.search(output), (
        f"Modelica class failed to resolve in {_example_id(example)}:\n"
        f"{output[-2000:]}"
    )

    if result.returncode != 0 and BOUNDS_LOST.search(output):
        raise BoundsLost(output[-1000:])

    # casadi >=3.8 recursion (see the module docstring): not this library's
    # defect, and not what these tests are guarding.
    if result.returncode != 0 and all(m in output for m in CASADI_RECURSION_MARKERS):
        pytest.skip(f"casadi/numpy recursion in casadi {casadi.__version__}")

    assert result.returncode == 0, output[-2000:]
