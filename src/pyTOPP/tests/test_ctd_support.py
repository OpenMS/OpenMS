import os
import tempfile
import textwrap

import pyopenms.pytopp.ctdsupport as cs


# ---------- tiny fakes to stand in for a CTDopts model/params ----------

class _P:
    """Param stub that mimics the attributes used by CTDHelpPrinter."""
    def __init__(self, name, typ="string", **kw):
        self.name = name
        self.type = typ
        self.required = kw.get("required", False)
        self.is_list = kw.get("is_list", False)
        self.description = kw.get("description", "")
        self.choices = kw.get("choices", None)
        self.restrictions = kw.get("restrictions", None)
        self.file_formats = kw.get("file_formats", None)
        self.tags = kw.get("tags", [])


class _Model:
    """
    Minimal stand-in for a CTDopts model. It exposes exactly what
    CTDHelpPrinter and PyTOPPArgParser use.
    """
    def __init__(self, params, defaults=None, name="Tool", desc=None, docurl=None):
        self._params = list(params)
        self.name = name
        self._defaults = dict(defaults or {})
        self.opt_attribs = {}
        if desc is not None:
            self.opt_attribs["description"] = desc
        if docurl is not None:
            self.opt_attribs["docurl"] = docurl

    def list_parameters(self):
        return list(self._params)

    def get_defaults(self):
        # May be a flat dict with colon keys; CTDHelpPrinter flattens anyway.
        return dict(self._defaults)

    def parse_cl_args(self, tokens, prefix="-"):
        """
        Simplified parser:
        - Accepts '-key value' and '-key=value'
        - Returns a dict containing ONLY flags seen (real CTDopts returns defaults too,
          but the parser filters using 'explicit' so this is fine).
        """
        out = {}
        for i, tok in enumerate(tokens):
            if not tok.startswith(prefix) or tok.startswith("+--"):
                continue
            if "=" in tok:
                key, val = tok[len(prefix):].split("=", 1)
                out[key] = val
            else:
                key = tok[len(prefix):]
                # If next token is a value-like (not another flag), capture it
                if i + 1 < len(tokens) and not tokens[i + 1].startswith("-"):
                    out[key] = tokens[i + 1]
                else:
                    out[key] = True  # bare flag
        return out

    def validate_args(self, mapping):
        # Identity validation for tests
        return dict(mapping)

    def write_ctd(self, path, values=None):
        # Write something small so code paths that call this won't fail
        with open(path, "w", encoding="utf-8") as fh:
            if values is None:
                fh.write("<ctd defaults/shape stub>\n")
            else:
                fh.write(str(values) + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# CTDHelpPrinter helper coverage
# ─────────────────────────────────────────────────────────────────────────────

def test_typename_choices_and_formats_helpers():
    assert cs.CTDHelpPrinter._typename(_P("b", bool)) == "bool"
    assert cs.CTDHelpPrinter._typename(_P("f", "input-file")) == "input-file"

    assert cs.CTDHelpPrinter._choices_text(_P("c", choices=["a", "b"])) == "a, b"
    assert cs.CTDHelpPrinter._choices_text(_P("r", restrictions=("x", "y"))) == "x, y"
    assert cs.CTDHelpPrinter._choices_text(_P("n")) is None

    assert cs.CTDHelpPrinter._formats_text(_P("fmt", file_formats=["mzML", "mzXML"])) == "mzML, mzXML"
    assert cs.CTDHelpPrinter._formats_text(_P("nofmt")) is None


def test_section_mapping_and_advanced_tagging():
    cfg = cs.CTDHelpConfig()

    assert cs.CTDHelpPrinter._section_for("in", cfg) == "I/O"
    assert cs.CTDHelpPrinter._section_for("infer:ipf_max_transition_pep", cfg) == "IPF"
    assert cs.CTDHelpPrinter._section_for("scoring:level", cfg) == "Scoring"
    assert cs.CTDHelpPrinter._section_for("weird:thing", cfg) == "General"

    assert cs.CTDHelpPrinter._is_advanced(_P("normal"), cfg) is False
    assert cs.CTDHelpPrinter._is_advanced(_P("adv", tags=["advanced"]), cfg) is True
    cfg.advanced_by_name.add("hard")
    assert cs.CTDHelpPrinter._is_advanced(_P("hard"), cfg) is True


def test_guess_usage_detects_in_and_out():
    params = [
        _P("in", "input-file", required=True),
        _P("out", "output-file"),
        _P("threads", "int"),
    ]
    model = _Model(params, name="Demo")
    printer = cs.CTDHelpPrinter(model)
    line = printer._guess_usage(cs.CTDHelpConfig(binary_name="demo"))
    assert "demo -in <input-file>" in line
    assert "[-out <output-file>]" in line


def test_guess_usage_when_no_io_params():
    params = [_P("threads", "int")]
    model = _Model(params, name="Demo")
    printer = cs.CTDHelpPrinter(model)
    line = printer._guess_usage(cs.CTDHelpConfig(binary_name="demo"))
    assert line.strip().endswith("[options]")


def test_print_help_smoke(capsys):
    params = [
        _P("in", "input-file", required=True, description="Input file"),
        _P("out", "output-file", description="Output file"),
        _P("scoring:level", "string", description="Level", choices=["ms1", "ms2"]),
        _P("flag", bool, description="Enable something"),
    ]
    defaults = {"out": "", "flag": False, "scoring:level": "ms2"}
    model = _Model(params, defaults, name="Printer", desc="A test tool", docurl="http://example")
    cs.CTDHelpPrinter(model, cfg=cs.CTDHelpConfig(binary_name="printer")).print(advanced=False)

    out = capsys.readouterr().out
    assert "Printer" in out
    assert "Usage:" in out
    assert "I/O:" in out
    assert "Scoring:" in out
    assert "Docs: http://example" in out
    # defaults & choices visible
    assert "default: ms2" in out
    assert "choices: ms1, ms2" in out
    # tip is printed in basic mode
    assert "help-advanced" in out


# ─────────────────────────────────────────────────────────────────────────────
# PyTOPPArgParser behavior
# ─────────────────────────────────────────────────────────────────────────────

def _build_parser_for_defaults():
    params = [
        _P("out", "output-file"),
        _P("threads", "int"),
        _P("scoring:level", "string"),
        _P("infer:run_protein", bool),
    ]
    defaults = {
        "out": "",
        "threads": "1",
        "scoring:level": "ms1",
        "infer:run_protein": "true",
    }
    model = _Model(params, defaults, name="ParserTool")
    return cs.PyTOPPArgParser(model)


def test_pure_parse_merges_defaults_ini_cli(monkeypatch):
    parser = _build_parser_for_defaults()

    # simulate no -ini (so defaults in model + CLI explicit only)
    argv = [
        "-scoring:level", "ms2",
        "-infer:run_protein", "false",
        "+--pi0_lambda", "0", "0", "0",   # forwarded; must be ignored by unknown-flag check
    ]
    result = parser.parse(argv, mode="pure")
    # CLI overrides defaults (deep-merge) and forwarded tokens do not interfere
    assert result["scoring:level"] == "ms2"
    assert result["infer:run_protein"] == "false"
    # unchanged default preserved
    assert result["threads"] == "1"


def test_find_unknown_flags_ignores_forwarded_and_numbers():
    parser = _build_parser_for_defaults()
    tokens = ["+--pi0_lambda", "0", "0", "0", "-unknown", "42", "-0.1"]
    unknowns = parser._find_unknown_flags(tokens)
    # only '-unknown' should be flagged
    assert "-unknown" in unknowns
    # forwarded and negative number literals are ignored
    assert all(u != "+--pi0_lambda" for u in unknowns)
    assert all(u != "-0.1" for u in unknowns)


def test_explicit_and_filter_cli_defaults_drop_non_explicit_and_nulls():
    parser = _build_parser_for_defaults()

    # Fake a CTDopts-style parse_cl_args result where one key is _Null
    null = cs._Null()  # sentinel imported by module
    cli_args = {
        "scoring": {"level": "ms1"},
        "unused": null,
    }
    explicit = {"scoring:level", "unused"}
    filtered = parser._filter_cli_defaults(cli_args, explicit)
    assert filtered == {"scoring:level": "ms1"}


def test_split_directives_extracts_and_preserves_rest():
    parser = _build_parser_for_defaults()
    argv = ["-ini", "file.ini", "-write_ini", "-write_param_ctd", "out.ctd", "-x", "1"]
    directives, rest = parser._split_directives(argv)

    assert directives["input_ctd"] == "file.ini"
    assert directives["write_tool_ctd"] is True
    assert directives["write_param_ctd"] == "out.ctd"
    assert rest == ["-x", "1"]


def test_merge_flat_and_helpers_roundtrip():
    parser = _build_parser_for_defaults()
    merged = parser._merge_flat({"a:b": 1}, {"a": {"c": 2}}, {"a:b": 3})
    assert merged == {"a:b": 3, "a:c": 2}

    nested = parser._as_nested({"x:y": 1, "x:z": 2, "w": 3})
    assert nested == {"x": {"y": 1, "z": 2}, "w": 3}

    flat_tuple = parser._flatten({"x": {"y": 1, "z": 2}, "w": 3})
    # keys are tuples when as_string=False
    assert flat_tuple == {("x", "y"): 1, ("x", "z"): 2, ("w",): 3}


def test_pyopenms_mode_returns_values_and_param_object(monkeypatch):
    """
    We don't require a real pyopenms Param; the code guards the sync and returns
    (validated, openms_param). Verify the tuple shape without raising.
    """
    parser = _build_parser_for_defaults()

    class FakeParam:
        def update(self, other, overwrite):
            pass

    fake_param = FakeParam()

    # Ensure the temporary file path is valid; the code is exception-safe anyway.
    out = parser._parse_pyopenms(["-threads", "8"], openms_param=fake_param)
    # _parse_pyopenms returns a tuple (validated, openms_param)
    assert isinstance(out, tuple) and len(out) == 2
    validated, returned_param = out
    assert returned_param is fake_param
    assert validated["threads"] == "8"
