## ----------------------------------------------------------------------------
## $Maintainer: $
## $Authors: $
## ----------------------------------------------------------------------------
"""
Boolean parameters in pyOpenMS (OpenMS issue #10116, step 6).

OpenMS has no boolean value type. A boolean parameter is a scalar string restricted to
exactly 'true' and 'false'; ``Param.ParamEntry.isBool()`` is the project-wide definition.
Such a parameter is a real Python ``bool`` on every read path and takes a real ``bool`` on
every write path, and nothing is inferred from a value's content -- so a value's Python
type cannot change just by copying it, which is what these tests pin down.
"""
import pytest

import pyopenms


def _tristate():
    """A Param with a parameter whose restrictions allow more than true/false."""
    p = pyopenms.Param()
    p.setValue("mode", "auto")
    p.setValidStrings("mode", ["auto", "true", "false"])
    return p


# --------------------------------------------------------------------------- reading


def test_flag_reads_as_bool():
    p = pyopenms.GaussFilter().getDefaults()
    assert p["use_ppm_tolerance"] is False
    assert p.getValue("use_ppm_tolerance") is False
    assert p.asDict()["use_ppm_tolerance"] is False
    assert dict(p.items())["use_ppm_tolerance"] is False
    assert p.get("use_ppm_tolerance") is False
    assert p.getEntry("use_ppm_tolerance").value is False
    assert p.isBool("use_ppm_tolerance")


def test_reversed_restriction_order_also_reads_as_bool():
    # FeatureFindingMetabo declares this one as {"false","true"} -- isBool() is
    # order-independent, so Python must not care either.
    p = pyopenms.FeatureFindingMetabo().getDefaults()
    assert p["report_convex_hulls"] is False
    assert p.isBool("report_convex_hulls")


def test_values_and_repr_show_bools():
    p = pyopenms.Param()
    p["flag"] = True
    assert p.values() == [True]
    assert repr(p) == "Param({'flag': True})"
    assert "flag = True" in str(p)


def test_unrestricted_true_string_stays_a_string():
    # No content sniffing: only the restrictions decide, never the value.
    p = pyopenms.Param()
    p.setValue("s", "true")
    assert p["s"] == "true"
    assert not p.getEntry("s").isBool()
    assert not p.isBool("s")


def test_tristate_stays_a_string():
    p = _tristate()
    assert p["mode"] == "auto"
    p.setValue("mode", "false")
    assert p["mode"] == "false"
    assert not p.isBool("mode")


def test_string_list_of_true_false_is_not_boolean():
    # isBool() requires a scalar; a list keeps both its value and its restrictions as str.
    p = pyopenms.Param()
    p.setValue("choices", ["true", "false"])
    p.setValidStrings("choices", ["true", "false"])
    assert p["choices"] == ["true", "false"]
    assert p.getValidStrings("choices") == ["true", "false"]
    assert not p.isBool("choices")


# --------------------------------------------------------------------------- writing


def test_bool_on_a_new_key_defines_a_boolean_parameter():
    p = pyopenms.Param()
    p["flag"] = True
    assert p["flag"] is True
    assert p.isBool("flag")
    # the stamped restrictions are what make boolean-ness survive copies, dict round trips
    # and INI round trips; getValidStrings() refuses to spell them out, isBool() reports it
    p2 = pyopenms.Param()
    p2.setValue("flag", False)
    assert p2["flag"] is False
    assert p2.isBool("flag")


def test_bool_on_an_existing_boolean_sets_only_the_value():
    p = pyopenms.GaussFilter().getDefaults()
    tags = p.getTags("use_ppm_tolerance")
    p["use_ppm_tolerance"] = True
    assert p["use_ppm_tolerance"] is True
    assert p.isBool("use_ppm_tolerance")
    assert p.getTags("use_ppm_tolerance") == tags


def test_setValue_with_bool_keeps_description_and_tags():
    p = pyopenms.Param()
    p.setValue("flag", True, "a flag", ["advanced"])
    assert p["flag"] is True
    assert p.isBool("flag")
    assert p.getDescription("flag") == "a flag"
    assert p.getTags("flag") == ["advanced"]


@pytest.mark.parametrize("setup", [
    lambda p: (p.setValue("mode", "auto"), p.setValidStrings("mode", ["auto", "true", "false"])),
    lambda p: p.setValue("mode", "hello"),          # unrestricted string
    lambda p: p.setValue("mode", 1),                # int
    lambda p: p.setValue("mode", 1.5),              # float
    lambda p: p.setValue("mode", ["a", "b"]),       # string list
])
def test_bool_on_a_non_boolean_entry_raises(setup):
    p = pyopenms.Param()
    setup(p)
    with pytest.raises(TypeError) as exc:
        p["mode"] = True
    assert "setValidStrings" in str(exc.value)
    # and nothing was changed
    assert not p.isBool("mode")


def test_string_assignment_to_a_boolean_is_still_unvalidated():
    # isBool() ignores the current value on purpose (test_ParamEntry_isBool.py), and a
    # hand-edited INI can produce this, so writing a string must keep working.
    p = pyopenms.Param()
    p["flag"] = True
    p.setValue("flag", "false")
    assert p["flag"] is False
    p.setValue("flag", "auto")
    assert p.isBool("flag")


def test_unsupported_value_type_raises_TypeError_naming_the_key():
    p = pyopenms.Param()
    with pytest.raises(TypeError) as exc:
        p["x"] = {1, 2}
    assert "'x'" in str(exc.value)


@pytest.mark.parametrize("value", [[1, True], [True, 1], [1.0, True], [True]])
def test_bools_inside_lists_are_rejected_at_every_position(value):
    p = pyopenms.Param()
    with pytest.raises(TypeError):
        p["x"] = value


# ------------------------------------------------------------------ restrictions API


def test_getValidStrings_raises_for_a_boolean():
    p = pyopenms.Param()
    p["flag"] = True
    with pytest.raises(TypeError) as exc:
        p.getValidStrings("flag")
    assert "isBool" in str(exc.value)
    with pytest.raises(TypeError):
        p.getEntry("flag").valid_strings


def test_getValidStrings_still_works_for_a_tristate():
    p = _tristate()
    assert p.getValidStrings("mode") == ["auto", "true", "false"]
    assert p.getEntry("mode").valid_strings == ["auto", "true", "false"]


def test_setValidStrings_turns_a_string_into_a_boolean():
    p = pyopenms.Param()
    p.setValue("s", "false")
    assert p["s"] == "false"
    p.setValidStrings("s", ["true", "false"])
    assert p["s"] is False
    assert p.isBool("s")


def test_getValueType_still_reports_the_stored_type():
    # No invented BOOL_VALUE: the value really is stored as a string.
    p = pyopenms.Param()
    p["flag"] = True
    assert p.getValueType("flag") == pyopenms.ValueType.STRING_VALUE


# --------------------------------------------------------------------------- round trips


def _roundtrip_param():
    p = pyopenms.Param()
    p["flag"] = True
    p["off"] = False
    p.setValue("mode", "auto")
    p.setValidStrings("mode", ["auto", "true", "false"])
    p.setValue("s", "true")
    p.setValue("n", 3)
    return p


def _assert_shapes(p):
    assert p["flag"] is True
    assert p["off"] is False
    assert p["mode"] == "auto"
    assert p["s"] == "true"
    assert p["n"] == 3
    assert p.isBool("flag") and p.isBool("off")
    assert not p.isBool("mode") and not p.isBool("s")


def test_dict_roundtrip_preserves_types():
    p = _roundtrip_param()
    _assert_shapes(pyopenms.Param(p.to_dict()))
    _assert_shapes(pyopenms.Param.from_dict(p.asDict()))


def test_repr_roundtrip_preserves_types():
    p = _roundtrip_param()
    _assert_shapes(eval(repr(p), {"Param": pyopenms.Param}))


def test_self_update_is_idempotent():
    p = _roundtrip_param()
    p.update(p)
    p.update(p.asDict())
    _assert_shapes(p)


def test_assigning_every_value_back_is_idempotent():
    # the shape of test_000.py's _testParam gauntlet
    p = _roundtrip_param()
    for k in p.keys():
        value = p[k]
        p[k] = value
        assert p[k] == value
    _assert_shapes(p)


def test_update_into_an_empty_param_keeps_boolean_ness():
    p = _roundtrip_param()
    q = pyopenms.Param()
    q.update(p)
    _assert_shapes(q)


def test_update_does_not_overwrite_the_targets_restrictions():
    p = _tristate()
    partial = pyopenms.Param()
    partial.setValue("mode", "false")
    p.update(partial)
    assert p["mode"] == "false"      # still a tri-state here, so still a str
    assert not p.isBool("mode")


def test_filtered_update_from_a_dict_adds_nothing():
    p = pyopenms.Param()
    p["flag"] = True
    before = p.asDict()
    p.update({"unknown_flag": True}, 1)
    assert p.asDict() == before
    assert not p.exists("unknown_flag")


def test_insert_and_copy_preserve_boolean_ness():
    p = _roundtrip_param()
    nested = pyopenms.Param()
    nested.insert("master:", p)
    assert nested["master:flag"] is True
    _assert_shapes(nested.copy("master:", True))


def test_ini_roundtrip_preserves_boolean_ness(tmp_path):
    p = _roundtrip_param()
    path = str(tmp_path / "roundtrip.ini")
    handler = pyopenms.ParamXMLFile()
    handler.store(path, p)
    loaded = pyopenms.Param()
    handler.load(path, loaded)
    _assert_shapes(loaded)


# ------------------------------------------------------------- algorithms and metadata


def test_setParameters_roundtrip_for_a_flag():
    gf = pyopenms.GaussFilter()
    p = gf.getParameters()
    p["use_ppm_tolerance"] = True
    gf.setParameters(p)
    assert gf.getParameters()["use_ppm_tolerance"] is True


def test_partial_param_with_a_bool_applies():
    gf = pyopenms.GaussFilter()
    partial = pyopenms.Param()
    partial["use_ppm_tolerance"] = True
    gf.setParameters(partial)
    assert gf.getParameters()["use_ppm_tolerance"] is True


def test_partial_param_cannot_turn_a_tristate_into_a_boolean():
    # Assigning True to a new key stamps {"true","false"} on it; the algorithm's own
    # restrictions must win, or a tri-state parameter would report as boolean afterwards.
    algo = pyopenms.OpenPepXLAlgorithm()
    partial = pyopenms.Param()
    partial["algorithm:deisotope"] = True
    algo.setParameters(partial)
    applied = algo.getParameters()
    assert applied["algorithm:deisotope"] == "true"
    assert not applied.isBool("algorithm:deisotope")
    assert applied.getValidStrings("algorithm:deisotope") == ["true", "false", "auto"]


def test_params_can_be_copied_into_meta_values():
    dp = pyopenms.DataProcessing()
    for key, value in pyopenms.GaussFilter().getDefaults().items():
        dp.setMetaValue(key, value)
    # DataValue carries no restrictions, so a flag comes back as the canonical string --
    # the same thing C++ writes via DefaultParamHandler::writeParametersToMetaValues.
    assert dp.getMetaValue("use_ppm_tolerance") == "false"


def test_setMetaValue_rejects_bools_inside_lists():
    dp = pyopenms.DataProcessing()
    with pytest.raises(TypeError):
        dp.setMetaValue("x", [1, True])


# ----------------------------------------------------- a boolean holding a bad value


def _broken_boolean():
    p = pyopenms.Param()
    p["flag"] = True
    p.setValue("flag", "auto")   # legal: isBool() ignores the current value
    return p


def test_reading_a_boolean_with_an_illegal_value_raises_ValueError():
    p = _broken_boolean()
    assert p.isBool("flag")
    for read in (lambda: p["flag"],
                 lambda: p.getValue("flag"),
                 lambda: p.get("flag"),
                 lambda: p.items(),
                 lambda: p.values(),
                 lambda: p.asDict(),
                 lambda: p.getEntry("flag").value):
        with pytest.raises(ValueError) as exc:
            read()
        assert "flag" in str(exc.value)


def test_str_and_repr_of_a_broken_boolean_do_not_raise():
    p = _broken_boolean()
    text = str(p)
    assert "'auto'" in text
    # the restrictions are the only hint why every accessor raises, so they stay visible
    assert "valid:" in text
    assert repr(p) == "Param({'flag': 'auto'})"
