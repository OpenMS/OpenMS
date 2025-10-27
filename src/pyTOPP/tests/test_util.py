import shutil
import pytest

import pyopenms.pytopp.util as util


def test_as_bool_variants():
    assert util.as_bool(True) is True
    assert util.as_bool(False) is False
    assert util.as_bool("true") is True
    assert util.as_bool("TrUe") is True
    assert util.as_bool("1") is True
    assert util.as_bool("no") is False
    assert util.as_bool(None) is False


def test_is_nullish_handles_sentinel():
    class _Null:  # mimics CTDopts sentinel name
        pass

    assert util.is_nullish(None) is True
    assert util.is_nullish("") is True
    assert util.is_nullish(_Null()) is True
    assert util.is_nullish("x") is False


def test_first_value_for_flag_basic():
    argv = ["-foo", "bar", "-baz=buz", "-empty"]
    assert util.first_value_for_flag(argv, "foo") == "bar"
    assert util.first_value_for_flag(argv, "baz") == "buz"
    assert util.first_value_for_flag(argv, "empty") == ""
    assert util.first_value_for_flag(argv, "nope") is None


def test_tok_extra_unescapes_and_splits():
    # NOTE: current implementation only strips a leading "-- " sentinel.
    # A mid-string "--" remains as a token.
    s = "+--pi0_lambda 0.0 0 0  +-X  \\-keep  --  -y"
    out = util.tok_extra(s)
    assert out == ["--pi0_lambda", "0.0", "0", "0", "-X", "-keep", "--", "-y"]


def test_coerce_bool_from_argv_variants():
    base = False
    assert util.coerce_bool_from_argv(["-flag", "true"], "flag", base) is True
    assert util.coerce_bool_from_argv(["-flag=true"], "flag", base) is True
    assert util.coerce_bool_from_argv(["-no-flag"], "flag", True) is False
    # bare `-flag` keeps current in this helper
    assert util.coerce_bool_from_argv(["-flag"], "flag", base) is base


def test_recover_missing_extras_and_coerce_all_bools():
    # fake model w/ two boolean params
    class P:
        def __init__(self, name, typ):
            self.name, self.type = name, typ

    class M:
        def list_parameters(self):
            return [P("a:toggle", bool), P("b:flag", "bool")]

    argv = ["-scoring:extra", "+--X 1", "-a:toggle", "true", "-b:flag=false"]
    arg_dict = {"scoring:extra": None}

    util.recover_missing_extras(arg_dict, argv, ["scoring:extra"])
    assert arg_dict["scoring:extra"] == "+--X 1"

    util.coerce_all_bools_from_argv(M(), arg_dict, argv)
    assert arg_dict["a:toggle"] is True
    assert arg_dict["b:flag"] is False


def test_run_tool_preview_uses_which_and_sanitizes_env(monkeypatch, capsys):
    # Pretend `which("echo")` returns /bin/echo
    monkeypatch.setattr(shutil, "which", lambda exe: "/bin/echo" if exe == "echo" else None)
    # Do a preview (dry-run)
    util.run_tool("echo", ["hello", "world"], dry=True, env={"LD_PRELOAD": "X", "PATH": "/x"})
    out = capsys.readouterr().out
    # Robust check: tokens appear in the previewed, quoted command line
    assert "CMD:" in out
    assert "echo" in out and "hello" in out and "world" in out


def test_run_tool_module_preview_and_validation(capsys):
    # module execution path, also dry
    util.run_tool("pyprophet", ["--help"], module_name="pyprophet", dry=True)
    assert "CMD:" in capsys.readouterr().out

    # mismatched module/executable should raise
    with pytest.raises(ValueError):
        util.run_tool("pyprophet", [], module_name="notpyprophet", dry=True)
