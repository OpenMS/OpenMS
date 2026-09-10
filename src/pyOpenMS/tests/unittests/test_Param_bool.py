# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# $Maintainer: Tom David Müller $
# $Authors: Tom David Müller $
# --------------------------------------------------------------------------

"""
Tests for Python bool support in Param and meta values.

OpenMS has no boolean parameter type: a flag is a string parameter holding
'true'/'false' restricted to valid strings ['true', 'false'] (the convention
TOPPBase, ParamXMLFile and the INI reader use). pyOpenMS accepts a Python
bool on assignment and stores it by that convention:

* Param.setValue/[]/update/from_dict and ParamEntry accept True/False and
  store 'true'/'false'; a key without string restrictions is marked as a
  boolean parameter (valid strings ['true', 'false'])
* setMetaValue() accepts True/False and stores 'true'/'false' as well
* reading keeps returning the string, so existing ``== "true"`` checks and
  ``value in valid_strings`` keep working; Param.getBool() is the explicit
  conversion (like the C++ ParamValue::toBool())
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
    assert p.getValue(key) == ("true" if expected else "false")
    assert p[key] == ("true" if expected else "false")
    assert p.getBool(key) is expected


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

    def test_param_entry_ctor_bool(self):
        """ParamEntry(name, True, ...) creates a boolean entry; a string does not."""
        e = pyopenms.ParamEntry("flag", True, "a flag")
        assert e.value == "true"
        assert e.valid_strings == BOOL_VALID
        e2 = pyopenms.ParamEntry("s", "true", "a string")
        assert e2.value == "true"
        assert e2.valid_strings == []

    def test_param_entry_value_setter(self):
        """Assigning a bool to ParamEntry.value marks an unrestricted entry as boolean."""
        e = pyopenms.ParamEntry("x", "abc", "")
        e.value = False
        assert e.value == "false"
        assert e.valid_strings == BOOL_VALID
        e.value = True
        assert e.value == "true"
        f = pyopenms.ParamEntry("y", "abc", "")
        f.value = "def"
        assert f.value == "def"
        assert f.valid_strings == []

    def test_existing_restrictions_are_kept(self):
        """A bool assigned to an entry with other valid strings leaves the restriction alone."""
        p = pyopenms.Param()
        p.setValue("mode", "none")
        p.setValidStrings("mode", ["none", "fixed"])
        p["mode"] = True
        assert p.getValidStrings("mode") == ["none", "fixed"]
        assert p["mode"] == "true"
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

    def test_reads_stay_strings(self):
        """Every read path keeps returning the flag strings (no silent bool conversion)."""
        p = self._make()
        assert p.getValue("sec:on") == "true"
        assert p["sec:off"] == "false"
        assert p.get("sec:on") == "true"
        assert p.get("missing", "dflt") == "dflt"
        d = p.asDict()
        assert d == {"sec:on": "true", "sec:off": "false", "n": 3, "s": "text"}
        assert p.to_dict() == d
        assert dict(p.items()) == d
        assert p.values() == [d[k] for k in p.keys()]
        assert p.getEntry("sec:on").value == "true"
        assert {e.name: e.value for e in p}["on"] == "true"
        # values remain members of their restrictions
        assert p["sec:on"] in p.getValidStrings("sec:on")

    def test_getBool(self):
        """getBool converts 'true'/'false' to bool, independent of restrictions."""
        p = self._make()
        assert p.getBool("sec:on") is True
        assert p.getBool("sec:off") is False
        p["plain"] = "true"  # no restrictions
        assert p.getBool("plain") is True
        p["plain"] = "false"
        assert p.getBool("plain") is False

    def test_getBool_rejects_non_flags(self):
        """getBool raises for missing keys, non-string values and other strings."""
        p = self._make()
        with pytest.raises(Exception):
            p.getBool("missing")
        with pytest.raises(Exception):
            p.getBool("n")
        with pytest.raises(Exception):
            p.getBool("s")
        p["sec:on"] = "maybe"
        with pytest.raises(Exception):
            p.getBool("sec:on")

    def test_string_assignment_keeps_flag(self):
        """Assigning 'true'/'false' strings to a flag keeps its restriction."""
        p = pyopenms.Param()
        p["flag"] = False
        p["flag"] = "true"
        _assert_flag(p, "flag", True)
        p.setValue("flag", "false")
        _assert_flag(p, "flag", False)

    def test_round_trip_dict(self):
        """from_dict(to_dict()) preserves flag values."""
        p = self._make()
        q = pyopenms.Param.from_dict(p.to_dict())
        assert q == p
        assert q.to_dict() == p.to_dict()
        assert q.getBool("sec:on") is True and q.getBool("sec:off") is False

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
    def test_algorithm_flag_reads(self):
        """Algorithm default flags read back as strings and convert with getBool."""
        gf = pyopenms.GaussFilter()
        p = gf.getDefaults()
        assert p["use_ppm_tolerance"] == "false"
        assert p.getBool("use_ppm_tolerance") is False
        assert p.getBool("write_log_messages") is False

    def test_algorithm_flag_set_with_bool(self):
        """Setting an algorithm flag with a bool passes checkDefaults/setParameters."""
        gf = pyopenms.GaussFilter()
        p = gf.getParameters()
        p["use_ppm_tolerance"] = True
        assert p["use_ppm_tolerance"] == "true"
        assert p.getValidStrings("use_ppm_tolerance") == ["true", "false"]
        assert p.getValueType("use_ppm_tolerance") == pyopenms.ValueType.STRING_VALUE
        p.checkDefaults("GaussFilter", gf.getDefaults())
        gf.setParameters(p)
        assert gf.getParameters()["use_ppm_tolerance"] == "true"
        assert gf.getParameters().getBool("use_ppm_tolerance") is True
        p.setValue("use_ppm_tolerance", False)
        gf.setParameters(p)
        assert gf.getParameters().getBool("use_ppm_tolerance") is False

    def test_partial_param_with_string_flags(self):
        """A partial Param with string flags behaves the same as one built from defaults."""
        gf = pyopenms.GaussFilter()
        p = pyopenms.Param()
        p["use_ppm_tolerance"] = "false"
        gf.setParameters(p)
        q = gf.getParameters()
        assert q["use_ppm_tolerance"] == "false"
        assert q.getBool("use_ppm_tolerance") is False
        p["use_ppm_tolerance"] = True
        gf.setParameters(p)
        assert gf.getParameters()["use_ppm_tolerance"] == "true"
        assert gf.getParameters().getBool("use_ppm_tolerance") is True

    def test_algorithm_invalid_flag_value_still_rejected(self):
        """An invalid flag string is still rejected by setParameters."""
        gf = pyopenms.GaussFilter()
        p = gf.getParameters()
        p.setValue("use_ppm_tolerance", "maybe")
        assert p["use_ppm_tolerance"] == "maybe"
        with pytest.raises(RuntimeError):
            gf.setParameters(p)


class TestMetaValueBool:
    def test_setMetaValue_bool(self):
        """setMetaValue accepts a bool and stores the flag string."""
        dp = pyopenms.DataProcessing()
        dp.setMetaValue("flag", True)
        assert dp.getMetaValue("flag") == "true"
        dp.setMetaValue("flag", False)
        assert dp.getMetaValue("flag") == "false"

    def test_copy_parameters_into_metadata(self):
        """Copying algorithm parameters (including flags set as bool) into metadata works."""
        gf = pyopenms.GaussFilter()
        p = gf.getParameters()
        p["use_ppm_tolerance"] = True
        dp = pyopenms.DataProcessing()
        for key, value in p.items():
            dp.setMetaValue(key, value)
        assert dp.getMetaValue("use_ppm_tolerance") == "true"
        assert dp.getMetaValue("gaussian_width") == p["gaussian_width"]

    def test_bool_in_meta_list_rejected(self):
        """Lists containing bools are rejected by setMetaValue."""
        dp = pyopenms.DataProcessing()
        with pytest.raises(TypeError):
            dp.setMetaValue("l", [True, False])
