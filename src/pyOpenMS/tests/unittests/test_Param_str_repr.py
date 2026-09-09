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

import pyopenms


def _make_param():
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
    p = pyopenms.Param()
    assert repr(p) == "Param({})"
    assert str(p) == "Param({})"


def test_repr_mirrors_asDict():
    p = _make_param()
    r = repr(p)
    assert r.startswith("Param({") and r.endswith("})")
    assert eval(r[len("Param("):-1]) == p.asDict()


def test_str_one_line_per_entry():
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
    p = _make_param()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print(p)
    assert buf.getvalue() == str(p) + "\n"


def test_str_of_algorithm_defaults_does_not_raise():
    # a real, large parameter set with sections, tags, restrictions and long descriptions
    p = pyopenms.PeakPickerHiRes().getDefaults()
    s = str(p)
    assert len(s.split("\n")) == p.size()
    assert "signal_to_noise = " in s
    assert "Param({" in repr(p)


def test_param_entry():
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
