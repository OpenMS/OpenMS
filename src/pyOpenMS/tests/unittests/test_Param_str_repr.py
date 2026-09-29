## ----------------------------------------------------------------------------
## $Maintainer: $
## $Authors: Timo Sachsenberg $
## ----------------------------------------------------------------------------
"""
Tests for the human-readable string forms of Param and ParamEntry
(OpenMS issue #10097: ``print(param)`` used to require ``print(param.asDict())``).

* ``str(param)`` / ``print(param)`` lists one entry per line as
  ``key = value (restrictions) [tags]  # description``
* ``repr(param)`` gives ``Param({...})`` mirroring ``Param.asDict()``
* ``ParamEntry`` has matching ``__str__`` / ``__repr__`` (using its leaf name)
"""

import contextlib
import io
import math

import pyopenms


def _make_param():
    """Build a small Param with sections, restrictions, tags and descriptions."""
    p = pyopenms.Param()
    p.setValue("mode", "fast", "Processing mode")
    p.setValidStrings("mode", ["fast", "exact"])
    p.setValue("ms_levels", [1, 2], "")
    p.setValue("name", "unnamed")
    p.setValue("algorithm:threshold", 0.5, "Intensity threshold", ["advanced"])
    p.setMinFloat("algorithm:threshold", 0.0)
    p.setMaxFloat("algorithm:threshold", 1.0)
    p.setValue("algorithm:iterations", 10, "Number of\niterations")
    p.setMinInt("algorithm:iterations", 1)
    return p


def test_empty_param():
    """An empty Param has the same str() and repr()."""
    p = pyopenms.Param()
    assert repr(p) == "Param({})"
    assert str(p) == "Param({})"


def test_repr_mirrors_asDict():
    """repr() contains exactly the dict that asDict() returns."""
    p = _make_param()
    r = repr(p)
    assert r.startswith("Param({") and r.endswith("})")
    assert eval(r[len("Param("):-1]) == p.asDict()


def test_repr_is_evaluable():
    """eval(repr(p)) reconstructs all keys and values via Param(dict)."""
    # repr() -> Param(dict) reconstructs all keys and values
    p = _make_param()
    p2 = eval(repr(p), {"Param": pyopenms.Param})
    assert isinstance(p2, pyopenms.Param)
    assert p2.asDict() == p.asDict()
    assert repr(p2) == repr(p)


def test_dict_constructor():
    """Param(dict) matches Param.from_dict() and leaves other constructors intact."""
    d = {"a": 1, "algorithm:threshold": 0.5, "name": "x", "levels": [1, 2]}
    p = pyopenms.Param(d)
    assert p.asDict() == d
    assert p == pyopenms.Param.from_dict(d)
    assert pyopenms.Param({}) == pyopenms.Param()
    # the existing constructors are unaffected
    assert pyopenms.Param(p) == p


def test_str_one_line_per_entry():
    """str() emits one line per entry in Param iteration order."""
    p = _make_param()
    lines = str(p).split("\n")
    assert len(lines) == p.size()
    # Param iteration order: entries of a section come before its subsections
    assert lines == [
        "mode = 'fast' (valid: 'fast', 'exact')  # Processing mode",
        "ms_levels = [1, 2]",
        "name = 'unnamed'",
        "algorithm:threshold = 0.5 (min=0.0, max=1.0) [advanced]  # Intensity threshold",
        # multi-line descriptions are flattened so each entry stays on one line
        "algorithm:iterations = 10 (min=1)  # Number of iterations",
    ]


def test_print_uses_str():
    """print(param) writes str(param)."""
    p = _make_param()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print(p)
    assert buf.getvalue() == str(p) + "\n"


def test_str_of_algorithm_defaults_does_not_raise():
    """A real algorithm parameter set formats without error."""
    # a real, large parameter set with sections, tags, restrictions and long descriptions
    p = pyopenms.PeakPickerHiRes().getDefaults()
    s = str(p)
    assert len(s.split("\n")) == p.size()
    assert "signal_to_noise = " in s
    assert "Param({" in repr(p)


def test_repr_handles_non_finite_floats():
    """nan/inf values are spelled float('nan') etc. so repr() stays evaluable."""
    p = pyopenms.Param({"a": float("nan"), "b": float("inf"), "c": float("-inf"),
                        "l": [1.0, float("nan"), float("inf")]})
    r = repr(p)
    assert "float('nan')" in r and "float('inf')" in r and "float('-inf')" in r
    p2 = eval(r, {"Param": pyopenms.Param})
    assert math.isnan(p2["a"]) and p2["b"] == math.inf and p2["c"] == -math.inf
    assert p2["l"][0] == 1.0 and math.isnan(p2["l"][1]) and p2["l"][2] == math.inf
    assert "float('nan')" in repr(p.getEntry("a"))


def test_str_flattens_tags_and_descriptions():
    """Line breaks in tags or descriptions never break the one-line-per-entry format."""
    p = pyopenms.Param()
    p.setValue("k", "v", "first\nsecond", ["tag\nnext"])
    assert str(p) == "k = 'v' [tag next]  # first second"


def test_param_entry():
    """ParamEntry has matching __repr__/__str__ using its leaf name."""
    p = _make_param()
    e = p.getEntry("algorithm:threshold")
    assert repr(e) == (
        "ParamEntry(name='threshold', value=0.5, "
        "description='Intensity threshold', tags=['advanced'])"
    )
    assert str(e) == "threshold = 0.5 (min=0.0, max=1.0) [advanced]  # Intensity threshold"

    plain = p.getEntry("name")
    assert repr(plain) == "ParamEntry(name='name', value='unnamed', description='', tags=[])"
    assert str(plain) == "name = 'unnamed'"
