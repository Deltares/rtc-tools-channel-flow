"""Compile and solve every example, guarding against unit-resolution failures.

The library's Modelica source is otherwise untested: nothing else here parses a
`.mo` file. These tests exist to catch a regression in how `Modelica.Units.SI.*`
types resolve, which is why a missing-class error always fails, even where Bug 3
is expected.

One known upstream defect stops the examples short of a solution, so only unit
resolution is asserted unconditionally: pymoca drops `min`/`max` modifiers on
MSL-typed connector variables, so bounds arrive as `(-inf, inf)` and rtc-tools
rejects the goals (rtc-tools#1798 "Bug 3"). Confirmed on 0.11.2; 0.9.3 is
unaffected; appears fixed on 0.12.0, though that rests on one example
(`weir/basic` XPASSes while the other five stop at an earlier 0.12.0 failure).
The condition below is an open-ended `>= 0.10` because the first affected
release was not bisected. Marked expected, strictly: once CI's pymoca pin is
bumped past a release that fixes this, the XPASS fails the run and the marker
has to go.

With pymoca 0.9.3 all six examples solve, `rc=0`, on casadi 3.7 and 3.8.
"""

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

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
INVALID_SOLUTION = re.compile(
    r"too (?:few|many) degrees of freedom|"
    r"relative/absolute (?:power|head loss) error exceedence",
    re.IGNORECASE,
)


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
def test_example_solves(example, tmp_path):
    working_copy = tmp_path / "example"
    working_copy.mkdir()
    for folder in ("input", "model", "src"):
        shutil.copytree(
            example / folder,
            working_copy / folder,
            ignore=shutil.ignore_patterns("*.pymoca_cache", "*.pickle", "__pycache__"),
        )
    (working_copy / "output").mkdir()

    # example.py calls run_optimization_problem() at import time and resolves
    # its input/output folders relative to src/, so it has to run as a script.
    result = subprocess.run(
        [sys.executable, "example.py"],
        cwd=working_copy / "src",
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

    assert result.returncode == 0, output[-2000:]
    assert "Solver succeeded" in output, output[-2000:]
    assert not INVALID_SOLUTION.search(output), output[-4000:]

    if _example_id(example) == "weir/basic":
        with (working_copy / "output" / "timeseries_export.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))

        assert len(rows) > 1
        # Parameters of the example's free-flow weir: width=10 m, coefficient=0.61.
        c_weir = (2 * 9.81) ** 0.5 * (2 / 3) * 0.61 * 10
        for row in rows[1:]:
            flow = float(row["WeirFlow1"])
            assert -1e-6 <= flow <= 1.0 + 1e-6
            if flow > 1e-6:
                crest = float(row["branch_1_water_level"]) - (flow / c_weir) ** (2 / 3)
                assert 1.7 - 1e-5 <= crest <= 3.0 + 1e-5
                assert crest >= float(row["branch_2_water_level"]) - 1e-5
