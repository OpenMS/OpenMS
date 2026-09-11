"""Tests for Param.ParamEntry.isBool()."""
import pyopenms


def test_param_entry_isBool():
    p = pyopenms.Param()

    p.setValue("flag", "false")
    p.setValidStrings("flag", ["true", "false"])
    assert p.getEntry("flag").isBool()

    # reversed order and current value do not matter
    p.setValue("flag", "true")
    p.setValidStrings("flag", ["false", "true"])
    assert p.getEntry("flag").isBool()

    # unrestricted true/false strings are not boolean
    p.setValue("unrestricted", "true")
    assert not p.getEntry("unrestricted").isBool()

    # additional allowed values are not boolean
    p.setValue("tristate", "auto")
    p.setValidStrings("tristate", ["auto", "true", "false"])
    assert not p.getEntry("tristate").isBool()

    # other types are not boolean
    p.setValue("int", 1)
    assert not p.getEntry("int").isBool()
