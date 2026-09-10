# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# $Maintainer: Tom David Müller $
# $Authors: Tom David Müller $
# --------------------------------------------------------------------------

"""
Tests for Python bool support in Param.

OpenMS has no boolean parameter type: a flag is a string parameter holding
'true'/'false' restricted to valid strings ['true', 'false'] (the convention
TOPPBase, ParamXMLFile and the INI reader use). pyOpenMS translates between
that convention and Python bool:

* assigning True/False stores 'true'/'false' and, for a key without string
  restrictions, marks the entry as a boolean parameter
* boolean parameters read back as True/False from getValue(), [], get(),
  items(), values(), asDict(), to_dict() and ParamEntry.value
"""

import os
import tempfile

import pytest

import pyopenms

BOOL_VALID = ["true", "false"]


def _assert_flag(p, key, expected):
    """Assert that key is a boolean parameter of the given value in p."""
    assert p.getValueType(key) == pyopenms.ValueType.STRING_VALUE
    assert p.getValidStrings(key) == BOOL_VALID
    assert p.getEntry(key).valid_strings == BOOL_VALID
    assert p.getValue(key) is expected
    assert p[key] is expected


class TestParamBoolWrite:
    def test_setValue_bool(self):
        """setValue(key, True/False) stores a flag string and marks the entry as boolean."""
        p = pyopenms.Param()
        p.setValue("flag", True, "a flag")
        _assert_flag(p, "flag", True)
        assert p.getDescription("flag") == "a flag"
        p.setValue("flag", False)
        _assert_flag(p, "flag", False)

    def test_setValue_bool_with_tags(self):
        """setValue with description and tags keeps them alongside the boolean marker."""
        p = pyopenms.Param()
        p.setValue("flag", False, "a flag", ["advanced"])
        _assert_flag(p, "flag", False)
        assert p.hasTag("flag", "advanced")

    def test_setitem_bool(self):
        """p[key] = bool creates a boolean parameter."""
        p = pyopenms.Param()
        p["a"] = True
        p["b"] = False
        _assert_flag(p, "a", True)
        _assert_flag(p, "b", False)

    def test_from_dict_bool(self):
        """Param.from_dict accepts bool values next to other types."""
        p = pyopenms.Param.from_dict({"flag": False, "n": 1, "s": "x"})
        _assert_flag(p, "flag", False)
        assert p["n"] == 1
        assert p["s"] == "x"

    def test_update_dict_bool(self):
        """Param.update with a dict accepts bool values."""
        p = pyopenms.Param()
        p.update({"flag": True})
        _assert_flag(p, "flag", True)

    def test_update_param_copies_bool_marker(self):
        """Param.update from another Param keeps flags readable as bool."""
        src = pyopenms.Param()
        src["flag"] = True
        dst = pyopenms.Param()
        dst.update(src)
        _assert_flag(dst, "flag", True)

    def test_param_entry_ctor_bool(self):
        """ParamEntry(name, True, ...) creates a boolean entry; a string does not."""
        e = pyopenms.ParamEntry("flag", True, "a flag")
        assert e.value is True
        assert e.valid_strings == BOOL_VALID
        e2 = pyopenms.ParamEntry("s", "true", "a string")
        assert e2.value == "true"
        assert e2.valid_strings == []

    def test_param_entry_value_setter(self):
        """Assigning a bool to ParamEntry.value marks an unrestricted entry as boolean."""
        e = pyopenms.ParamEntry("x", "abc", "")
        e.value = False
        assert e.value is False
        assert e.valid_strings == BOOL_VALID
        e.value = True
        assert e.value is True
        f = pyopenms.ParamEntry("y", "abc", "")
        f.value = "def"
        assert f.value == "def"
        assert f.valid_strings == []

    def test_existing_restrictions_are_kept(self):
        """A bool assigned to an entry with other valid strings leaves the restriction alone."""
        # a bool assigned to an entry with other valid strings does not touch them
        p = pyopenms.Param()
        p.setValue("mode", "none")
        p.setValidStrings("mode", ["none", "fixed"])
        p["mode"] = True
        assert p.getValidStrings("mode") == ["none", "fixed"]
        assert p["mode"] == "true"  # not a boolean parameter -> stays a str
        # checkDefaults reports it like any other invalid string
        defaults = pyopenms.Param()
        defaults.setValue("mode", "none")
        defaults.setValidStrings("mode", ["none", "fixed"])
        with pytest.raises(Exception):
            p.checkDefaults("test", defaults)

    def test_bool_in_list_rejected(self):
        """Lists containing bools raise TypeError on every write path."""
        p = pyopenms.Param()
        with pytest.raises(TypeError):
            p.setValue("l", [True, False])
        with pytest.raises(TypeError):
            p["l"] = [True, False]
        with pytest.raises(TypeError):
            pyopenms.Param.from_dict({"l": [True, False]})
        assert "l" not in p

    def test_unsupported_type_message(self):
        """Unsupported value types raise a TypeError naming the accepted types."""
        p = pyopenms.Param()
        with pytest.raises(TypeError, match="bool"):
            p["x"] = object()


class TestParamBoolRead:
    def _make(self):
        """Build a Param with two flags, an int and a string."""
        p = pyopenms.Param()
        p["sec:on"] = True
        p["sec:off"] = False
        p["n"] = 3
        p["s"] = "text"
        return p

    def test_read_paths(self):
        """Every read path returns boolean parameters as Python bool."""
        p = self._make()
        assert p.getValue("sec:on") is True
        assert p.getValue("sec:off") is False
        assert p["sec:on"] is True
        assert p.get("sec:off") is False
        assert p.get("missing", "dflt") == "dflt"
        d = p.asDict()
        assert d["sec:on"] is True and d["sec:off"] is False
        assert d["n"] == 3 and d["s"] == "text"
        assert p.to_dict() == d
        assert dict(p.items()) == d
        assert p.values() == [d[k] for k in p.keys()]
        assert p.getEntry("sec:on").value is True
        assert {e.name: e.value for e in p}["on"] is True

    def test_string_true_on_flag_reads_as_bool(self):
        """Assigning 'true'/'false' strings to a flag keeps it a boolean parameter."""
        p = pyopenms.Param()
        p["flag"] = False
        p["flag"] = "true"  # plain string keeps the restriction -> still a flag
        _assert_flag(p, "flag", True)
        p.setValue("flag", "false")
        _assert_flag(p, "flag", False)

    def test_invalid_string_on_flag_is_returned_unchanged(self):
        """An invalid string on a flag reads back unchanged (str)."""
        p = pyopenms.Param()
        p["flag"] = True
        p["flag"] = "maybe"
        assert p["flag"] == "maybe"
        assert p.getValidStrings("flag") == BOOL_VALID

    def test_unrestricted_string_is_not_converted(self):
        """'true' on a key without the flag restriction stays a str."""
        p = pyopenms.Param()
        p["s"] = "true"
        assert p["s"] == "true"
        assert p.getValidStrings("s") == []
        p.setValidStrings("s", ["false", "true"])  # reversed order is not the convention
        assert p["s"] == "true"

    def test_other_valid_strings_stay_str(self):
        """String parameters with other valid strings are untouched."""
        p = pyopenms.Param()
        p["mode"] = "a"
        p.setValidStrings("mode", ["a", "b"])
        assert p["mode"] == "a"

    def test_round_trip_dict(self):
        """from_dict(to_dict()) preserves boolean parameters."""
        p = self._make()
        q = pyopenms.Param.from_dict(p.to_dict())
        assert q == p
        assert q.to_dict() == p.to_dict()
        _assert_flag(q, "sec:on", True)
        _assert_flag(q, "sec:off", False)

    def test_round_trip_ini(self):
        """ParamXMLFile store/load preserves boolean parameters."""
        p = self._make()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "p.ini")
            pyopenms.ParamXMLFile().store(path, p)
            q = pyopenms.Param()
            pyopenms.ParamXMLFile().load(path, q)
        _assert_flag(q, "sec:on", True)
        _assert_flag(q, "sec:off", False)
        assert q["n"] == 3 and q["s"] == "text"


class TestParamBoolAlgorithm:
    def test_algorithm_flag_reads_as_bool(self):
        """Algorithm default flags read back as bool."""
        gf = pyopenms.GaussFilter()
        p = gf.getDefaults()
        assert p["use_ppm_tolerance"] is False
        assert p.getValue("write_log_messages") is False
        assert p.asDict()["use_ppm_tolerance"] is False

    def test_algorithm_flag_set_with_bool(self):
        """Setting an algorithm flag with a bool passes checkDefaults/setParameters."""
        gf = pyopenms.GaussFilter()
        p = gf.getParameters()
        p["use_ppm_tolerance"] = True
        assert p.getValidStrings("use_ppm_tolerance") == BOOL_VALID
        assert p.getValueType("use_ppm_tolerance") == pyopenms.ValueType.STRING_VALUE
        p.checkDefaults("GaussFilter", gf.getDefaults())
        gf.setParameters(p)
        assert gf.getParameters()["use_ppm_tolerance"] is True
        p.setValue("use_ppm_tolerance", False)
        gf.setParameters(p)
        assert gf.getParameters()["use_ppm_tolerance"] is False

    def test_algorithm_invalid_flag_value_still_rejected(self):
        """An invalid flag string is still rejected by setParameters."""
        gf = pyopenms.GaussFilter()
        p = gf.getParameters()
        p.setValue("use_ppm_tolerance", "maybe")
        assert p["use_ppm_tolerance"] == "maybe"
        with pytest.raises(RuntimeError):
            gf.setParameters(p)
