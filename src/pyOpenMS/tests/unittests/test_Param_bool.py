# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# $Maintainer: Tom David Müller $
# $Authors: Tom David Müller $
# --------------------------------------------------------------------------

"""
Tests for Python bool support in Param.

OpenMS has no boolean parameter type: a flag is a string parameter holding
'true'/'false', usually restricted to exactly these two values. The C++ side
reads flags by content (ParamValue::toBool) and pyOpenMS does the same, so a
boolean parameter is a Python bool on every path, whatever its origin:

* a string parameter holding 'true'/'false' reads back as True/False from
  getValue(), [], get(), items(), values(), asDict(), to_dict() and
  ParamEntry.value, unless its restrictions allow other values as well
* its restrictions are [True, False] and getValueType() is BOOL_VALUE
* assigning True/False stores 'true'/'false' and gives a key without
  restrictions the restriction [True, False] (for INI/CTD rendering)
"""

import os
import tempfile

import pytest

import pyopenms

BOOL_VALID = [True, False]
BOOL_TYPE = pyopenms.ValueType.BOOL_VALUE
STR_TYPE = pyopenms.ValueType.STRING_VALUE


def _assert_flag(p, key, expected, valid=BOOL_VALID):
    """Assert that key reads as the given bool on every path of p."""
    assert p.getValueType(key) == BOOL_TYPE
    assert p.getValidStrings(key) == valid
    assert p.getEntry(key).valid_strings == valid
    assert p.getValue(key) is expected
    assert p[key] is expected
    assert p.get(key) is expected
    assert p.getEntry(key).value is expected
    assert p.asDict()[key] is expected


def _ini_round_trip(p):
    """Store p as INI and load it back."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "p.ini")
        pyopenms.ParamXMLFile().store(path, p)
        q = pyopenms.Param()
        pyopenms.ParamXMLFile().load(path, q)
        return q, open(path).read()


class TestParamBoolWrite:
    def test_setValue_bool(self):
        """setValue(key, True/False) stores a flag and marks a new key as boolean."""
        p = pyopenms.Param()
        p.setValue("flag", True, "a flag")
        _assert_flag(p, "flag", True)
        assert p.getDescription("flag") == "a flag"
        p.setValue("flag", False)
        _assert_flag(p, "flag", False)

    def test_setValue_bool_with_tags(self):
        """setValue with description and tags keeps them alongside the restriction."""
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

    def test_update_param_copies_restriction(self):
        """Param.update from another Param keeps the flag restriction."""
        src = pyopenms.Param()
        src["flag"] = True
        dst = pyopenms.Param()
        dst.update(src)
        _assert_flag(dst, "flag", True)

    def test_param_entry_ctor_bool(self):
        """ParamEntry(name, True, ...) and ParamEntry(name, 'true', ...) both read as bool."""
        e = pyopenms.ParamEntry("flag", True, "a flag")
        assert e.value is True
        assert e.valid_strings == BOOL_VALID
        e2 = pyopenms.ParamEntry("s", "true", "a string")
        assert e2.value is True
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

    def test_bool_on_entry_with_other_restrictions(self):
        """A bool assigned to a tri-state parameter is stored but stays a str (restrictions win)."""
        p = pyopenms.Param()
        p.setValue("mode", "none")
        p.setValidStrings("mode", ["none", "fixed"])
        p["mode"] = True
        assert p.getValidStrings("mode") == ["none", "fixed"]
        assert p["mode"] == "true"
        assert p.getValueType("mode") == STR_TYPE
        # checkDefaults reports it like any other invalid string
        defaults = pyopenms.Param()
        defaults.setValue("mode", "none")
        defaults.setValidStrings("mode", ["none", "fixed"])
        with pytest.raises(Exception):
            p.checkDefaults("test", defaults)

    def test_bool_in_list_rejected(self):
        """Lists containing bools raise TypeError on every write path and at every position."""
        p = pyopenms.Param()
        for bad in ([True, False], [1, True], [1.0, False], [2, 3, True]):
            with pytest.raises(TypeError):
                p.setValue("l", bad)
            with pytest.raises(TypeError):
                p["l"] = bad
            with pytest.raises(TypeError):
                pyopenms.Param.from_dict({"l": bad})
        assert "l" not in p
        # meta values are out of scope: they take no bool, but the list guard applies there too
        s = pyopenms.MSSpectrum()
        with pytest.raises(TypeError):
            s.setMetaValue("b", True)
        with pytest.raises(TypeError):
            s.setMetaValue("l", [1, True])

    def test_unsupported_type_message(self):
        """Unsupported value types raise a TypeError naming the accepted types."""
        p = pyopenms.Param()
        with pytest.raises(TypeError, match="bool"):
            p["x"] = object()


class TestParamBoolRestrictions:
    def test_valid_strings_of_flag_are_bools(self):
        """Restrictions of a boolean parameter are [True, False] in either declared order."""
        p = pyopenms.Param()
        p.setValue("a", "true")
        p.setValidStrings("a", ["true", "false"])
        p.setValue("b", "true")
        p.setValidStrings("b", ["false", "true"])
        for key in ("a", "b"):
            _assert_flag(p, key, True)
            assert True in p.getValidStrings(key)
            assert p[key] in p.getValidStrings(key)

    def test_setValidStrings_accepts_bools(self):
        """setValidStrings(key, [True, False]) (any order) stores the canonical OpenMS flag restriction."""
        for order in ([True, False], [False, True]):
            p = pyopenms.Param()
            p["s"] = "false"
            p.setValidStrings("s", order)
            _assert_flag(p, "s", False)
            q, ini = _ini_round_trip(p)
            assert 'type="bool"' in ini
            _assert_flag(q, "s", False)

    def test_param_entry_valid_strings_setter(self):
        """ParamEntry.valid_strings accepts bools and reads back as [True, False]."""
        e = pyopenms.ParamEntry("s", "true", "")
        e.valid_strings = [False, True]
        assert e.valid_strings == BOOL_VALID
        assert e.value is True
        e.valid_strings = ["x", "y"]
        assert e.valid_strings == ["x", "y"]
        assert e.value == "true"

    def test_invalid_restriction_elements(self):
        """Restriction elements other than str/bytes/bool raise TypeError."""
        p = pyopenms.Param()
        p["s"] = "x"
        with pytest.raises(TypeError):
            p.setValidStrings("s", [1, 2])
        e = pyopenms.ParamEntry("s", "x", "")
        with pytest.raises(TypeError):
            e.valid_strings = [1.5]
        # a single string is not a list of valid strings
        with pytest.raises(TypeError):
            p.setValidStrings("s", "abc")
        with pytest.raises(TypeError):
            e.valid_strings = "abc"

    def test_other_valid_strings_stay_str(self):
        """String parameters with other valid strings are untouched."""
        p = pyopenms.Param()
        p["mode"] = "a"
        p.setValidStrings("mode", ["a", "b"])
        assert p["mode"] == "a"
        assert p.getValidStrings("mode") == ["a", "b"]
        assert p.getValueType("mode") == STR_TYPE

    def test_tri_state_stays_str_everywhere(self):
        """'true' on a parameter that also allows 'auto' stays a str on every path."""
        p = pyopenms.Param()
        p.setValue("tri", "auto")
        p.setValidStrings("tri", ["auto", "true", "false"])
        p["tri"] = "true"
        assert p["tri"] == "true"
        assert p.getValue("tri") == "true"
        assert p.getEntry("tri").value == "true"
        assert p.asDict()["tri"] == "true"
        assert p.getValueType("tri") == STR_TYPE
        assert p.getValidStrings("tri") == ["auto", "true", "false"]
        p["tri"] = "auto"
        assert p["tri"] == "auto"

    def test_type_follows_current_restrictions_in_any_order(self):
        """The Python type is derived from the entry's current value and restrictions, whatever the call order."""
        # value first, restriction later: reads as bool until 'auto' is allowed as well
        p = pyopenms.Param()
        p["tri"] = "true"
        assert p["tri"] is True
        assert p.getValueType("tri") == BOOL_TYPE
        p.setValidStrings("tri", ["auto", "true", "false"])
        assert p["tri"] == "true"
        assert p.getValueType("tri") == STR_TYPE
        # restriction first, value later: str from the start
        q = pyopenms.Param()
        q.setValue("tri", "auto")
        q.setValidStrings("tri", ["auto", "true", "false"])
        q["tri"] = "true"
        assert q["tri"] == "true"
        # narrowing the restriction back to the two bool strings makes it a bool again
        q.setValidStrings("tri", ["true", "false"])
        assert q["tri"] is True
        assert q.getValidStrings("tri") == BOOL_VALID
        # at every point all read paths agree with each other
        for param in (p, q):
            v = param["tri"]
            for other in (param.getValue("tri"), param.getEntry("tri").value, param.asDict()["tri"]):
                assert other == v and type(other) is type(v)

    def test_getValueType(self):
        """getValueType() is BOOL_VALUE exactly when the value reads back as bool."""
        p = pyopenms.Param()
        p["a"] = True
        p["b"] = "false"
        p["c"] = "hello"
        p["d"] = 1
        assert p.getValueType("a") == BOOL_TYPE
        assert p.getValueType("b") == BOOL_TYPE
        assert p.getValueType("c") == STR_TYPE
        assert p.getValueType("d") == pyopenms.ValueType.INT_VALUE
        p["a"] = "maybe"
        assert p.getValueType("a") == STR_TYPE
        assert p["a"] == "maybe"
        assert int(BOOL_TYPE) == 7
        assert BOOL_TYPE not in (
            pyopenms.ValueType.STRING_VALUE, pyopenms.ValueType.INT_VALUE,
            pyopenms.ValueType.DOUBLE_VALUE, pyopenms.ValueType.STRING_LIST,
            pyopenms.ValueType.INT_LIST, pyopenms.ValueType.DOUBLE_LIST,
            pyopenms.ValueType.EMPTY_VALUE)


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

    def test_string_true_reads_as_bool(self):
        """'true'/'false' strings read back as bool, with and without the restriction."""
        p = pyopenms.Param()
        p["flag"] = False
        p["flag"] = "true"
        _assert_flag(p, "flag", True)
        p.setValue("flag", "false")
        _assert_flag(p, "flag", False)
        p["s"] = "true"
        _assert_flag(p, "s", True, valid=[])
        p["s"] = "false"
        _assert_flag(p, "s", False, valid=[])
        assert pyopenms.Param.from_dict({"f": True}) == pyopenms.Param.from_dict({"f": "true"})

    def test_invalid_string_on_flag_is_returned_unchanged(self):
        """An invalid string on a flag reads back unchanged (str)."""
        p = pyopenms.Param()
        p["flag"] = True
        p["flag"] = "maybe"
        assert p["flag"] == "maybe"
        assert p.getValueType("flag") == STR_TYPE
        assert p.getValidStrings("flag") == BOOL_VALID

    def test_round_trip_dict(self):
        """from_dict(to_dict()) preserves boolean parameters."""
        p = self._make()
        q = pyopenms.Param.from_dict(p.to_dict())
        assert q == p
        assert q.to_dict() == p.to_dict()
        _assert_flag(q, "sec:on", True)
        _assert_flag(q, "sec:off", False)

    def test_round_trip_ini(self):
        """ParamXMLFile store/load preserves boolean parameters, restricted or not."""
        p = self._make()
        p["plain"] = "true"  # no restriction: stored as type="string"
        q, ini = _ini_round_trip(p)
        _assert_flag(q, "sec:on", True)
        _assert_flag(q, "sec:off", False)
        _assert_flag(q, "plain", True, valid=[])
        assert q["n"] == 3 and q["s"] == "text"


class TestParamBoolConsistency:
    """The same flag must read as bool no matter how the Param was built."""

    def test_reversed_restriction_order_in_cpp_defaults(self):
        """C++ classes declaring ['false', 'true'] read like those declaring ['true', 'false']."""
        gf = pyopenms.GaussFilter().getDefaults()
        _assert_flag(gf, "use_ppm_tolerance", False)
        ffm = pyopenms.FeatureFindingMetabo().getDefaults()
        _assert_flag(ffm, "report_convex_hulls", False)
        _assert_flag(ffm, "enable_RT_filtering", True)
        _assert_flag(pyopenms.PeakIntegrator().getDefaults(), "fit_EMG", False)

    def test_unrestricted_cpp_default(self):
        """A C++ default declared as 'true' without restriction reads as bool."""
        p = pyopenms.MRMFeatureFinderScoring().getDefaults()
        _assert_flag(p, "use_ms1_ion_mobility", True, valid=[])

    def test_partial_param_setParameters(self):
        """A partial Param with string- or bool-valued flags reads as bool before and after setParameters."""
        for value in ("false", False):
            p = pyopenms.Param()
            p["use_ppm_tolerance"] = value
            assert p["use_ppm_tolerance"] is False
            gf = pyopenms.GaussFilter()
            gf.setParameters(p)
            assert gf.getParameters()["use_ppm_tolerance"] is False
            assert gf.getParameters()["write_log_messages"] is False
        for value in ("true", True):
            p = pyopenms.Param()
            p["use_ppm_tolerance"] = value
            gf = pyopenms.GaussFilter()
            gf.setParameters(p)
            assert gf.getParameters()["use_ppm_tolerance"] is True

    def test_setDefaults(self):
        """Param.setDefaults on a partial Param leaves every flag readable as bool."""
        p = pyopenms.Param()
        p["use_ppm_tolerance"] = "false"
        p.setDefaults(pyopenms.GaussFilter().getDefaults())
        assert p["use_ppm_tolerance"] is False
        assert p["write_log_messages"] is False
        assert p.getValueType("use_ppm_tolerance") == BOOL_TYPE

    def test_merge_and_insert_both_directions(self):
        """merge()/insert() give bool whichever side carried the restriction."""
        a = pyopenms.Param()
        a["flag"] = "true"
        b = pyopenms.Param()
        b["flag"] = True
        ab = pyopenms.Param(a)
        ab.merge(b)
        ba = pyopenms.Param(b)
        ba.merge(a)
        assert ab["flag"] is True and ba["flag"] is True
        c = pyopenms.Param()
        c.insert("", a)
        d = pyopenms.Param()
        d.insert("", b)
        assert c["flag"] is True and d["flag"] is True
        assert c.getValueType("flag") == BOOL_TYPE

    def test_ini_round_trip_of_cpp_defaults(self):
        """INI store/load of defaults gives bool for both restriction orders."""
        q, _ = _ini_round_trip(pyopenms.FeatureFindingMetabo().getDefaults())
        _assert_flag(q, "report_convex_hulls", False)
        q, _ = _ini_round_trip(pyopenms.GaussFilter().getDefaults())
        _assert_flag(q, "use_ppm_tolerance", False)

    def test_idfilter_like_enum_with_false_member(self):
        """A parameter whose restrictions are ['false', 'sequence', ...] is a str on every path."""
        p = pyopenms.Param()
        p.setValue("best", "false")
        p.setValidStrings("best", ["false", "sequence", "sequence+charge"])
        assert p["best"] == "false"
        assert p.getValueType("best") == STR_TYPE
        assert p.getValidStrings("best") == ["false", "sequence", "sequence+charge"]
        p["best"] = "sequence"
        assert p["best"] == "sequence"


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
        assert p.getValueType("use_ppm_tolerance") == BOOL_TYPE
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
