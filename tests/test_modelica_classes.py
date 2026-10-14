"""Resolve every class in the Modelica library.

The examples instantiate only PumpingStation, Orifice, Weir and Storage.Linear,
which leaves most of the library unexercised -- including the files that carry
the most unit references, such as Hydraulic/Branches/LinearisedSV.mo and the
Salt elements. Flattening every class covers them.

These are hard assertions on all pymoca versions: a class that cannot resolve
means a type name is wrong or an import is missing, which is a defect in this
library regardless of which pymoca is installed.
"""

import copy
import importlib.resources
import inspect
import os
from pathlib import Path

import pytest
from pymoca import ast, parser
from pymoca.tree import flatten

LIBRARY = Path(__file__).parent.parent / "src" / "rtctools_channel_flow" / "modelica"


def _parse(text):
    # bypass_cache arrived with pymoca 0.11; 0.9 has no text cache at all.
    if "bypass_cache" in inspect.signature(parser.parse).parameters:
        return parser.parse(text, bypass_cache=True)
    return parser.parse(text)


@pytest.fixture(scope="module")
def tree():
    """The MSL plus this library, parsed once."""
    msl = Path(str(importlib.resources.files("rtctools_standard_library").joinpath("modelica")))
    root = None
    for folder in (msl, LIBRARY):
        for dirpath, _dirs, files in os.walk(folder):
            for name in sorted(f for f in files if f.endswith(".mo")):
                parsed = _parse(Path(dirpath, name).read_text(encoding="utf-8"))
                if root is None:
                    root = parsed
                else:
                    root.extend(parsed)
    return root


def _classes():
    # package.mo files declare a package rather than an instantiable class.
    return sorted(p for p in LIBRARY.rglob("*.mo") if p.name != "package.mo")


def _class_name(path):
    return str(path.relative_to(LIBRARY).with_suffix("")).replace(os.sep, ".")


@pytest.mark.parametrize("model", _classes(), ids=_class_name)
def test_class_resolves(tree, model):
    # flatten() resolves with find_class(copy=False) and mutates the tree, so
    # each class needs its own copy or the result depends on test order.
    flatten(copy.deepcopy(tree), ast.ComponentRef.from_string(_class_name(model)))
