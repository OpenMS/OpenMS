from __future__ import annotations

import os
import shutil
import sys
import tempfile
import textwrap
from collections.abc import Mapping
from dataclasses import dataclass, field
from difflib import get_close_matches
from typing import Iterable, Mapping, Sequence

import CTDopts
import pyopenms as pms

# add near the top with the other CTDopts imports
from CTDopts.CTDopts import (
    CTDModel,
    _Null,
    args_from_file,
    flatten_dict,
    override_args,
    parse_cl_directives,
    set_nested_key,
)


@dataclass
class CTDHelpConfig:
    """Configuration for pretty-printing CTDopts-based CLI help.

    This config controls how `print_model_help` renders the two-column help output
    for a CTDopts `CTDModel`. You can customize cosmetics, column widths, which
    parameters count as "advanced", and how options are grouped into sections.

    Attributes:
        tool_name:
            Human-friendly title printed at the top. If omitted, falls back to
            `model.name` or `"Tool"`.
        binary_name:
            Program name shown in the `Usage:` line (e.g., `"PyProphet.py"`).
            If omitted, uses the basename of `sys.argv[0]`.
        show_directives:
            If `True`, append the OpenMS CTD/INI directives (-ini, -write_ini,
            -write_param_ctd) at the end of the help.
        show_tip_when_basic:
            If `True` and `advanced=False`, print a tip about
            `--help-advanced`/`--helphelp`.
        show_description:
            If `True`, print the model's description (if any) below the tool name.
        left_col_width:
            Width (in characters) of the left column containing option names.
            Long option names wrap onto the next line.
        min_total_width:
            Minimum terminal width assumed when wrapping. If the detected TTY
            width is narrower, this value is used to keep output readable.
        col_padding:
            Number of spaces between the left and right columns.
        advanced_tag:
            Parameter tag that marks an option as *advanced*. Any parameter
            whose CTD model includes this tag will be hidden in basic help.
        advanced_by_name:
            Explicit set of parameter names to treat as advanced regardless of
            tags (useful if you cannot or do not want to modify the CTD model).
        section_order:
            Order in which logical sections are rendered.
        section_map:
            Maps the top-level prefix of a parameter name (before the first
            colon) to a section name. For example, `"scoring:level"` -> prefix
            `"scoring"` -> section `"Scoring"`.
        io_param_names:
            Parameter names that are always placed in the *I/O* section (e.g.,
            `"in"`, `"out"`, `"input"`, `"output"`).
    """

    # Cosmetic
    tool_name: str | None = None
    binary_name: str | None = None  # e.g. "pyprophet" (default: sys.argv[0])
    show_directives: bool = True
    show_tip_when_basic: bool = True
    show_description: bool = True

    # Columns / layout
    left_col_width: int = 34  # width of the left key column
    min_total_width: int = 80  # minimum terminal width to assume
    col_padding: int = 2  # spaces between columns

    # Advanced parameter selection
    advanced_tag: str = "advanced"  # any param with this tag is "advanced"
    advanced_by_name: set[str] = field(default_factory=set)

    # Section grouping
    section_order: Sequence[str] = (
        "I/O",
        "General",
        "Scoring",
        "Inference",
        "Export",
    )
    # Map param-name prefixes (before first ':') to section names
    section_map: Mapping[str, str] = field(
        default_factory=lambda: {
            "scoring": "Scoring",
            "infer": "Inference",
            "export": "Export",
        }
    )
    # Names that always go to I/O
    io_param_names: set[str] = field(
        default_factory=lambda: {"in", "out", "input", "output"}
    )

class CTDHelpPrinter:
    """
    Pretty, two-column help renderer for CTDopts `CTDModel`s.

    Responsibilities:
      - Group parameters into sections (I/O, Scoring, Inference, Export, etc.)
      - Show a friendly `Usage:` line using simple heuristics.
      - Render defaults, choices, file formats, and "required" badges.
      - Hide "advanced" parameters unless requested.
      - Optionally print OpenMS CTD/INI directive flags.

    Usage:
        printer = CTDHelpPrinter(model, cfg=CTDHelpConfig(...))
        printer.print(advanced=False)  # basic
        printer.print(advanced=True)   # show all
    """

    def __init__(self, model, cfg: "CTDHelpConfig | None" = None) -> None:
        """
        Args:
            model: A CTDopts `CTDModel` instance.
            cfg:   Optional `CTDHelpConfig` to tweak layout/sections/labels.
        """
        self.model = model
        self.cfg = cfg

    # --------------------------- public API ---------------------------

    def print(self, *, advanced: bool) -> None:
        """Render help to stdout.

        Args:
            advanced: If True, include advanced parameters; else hide them.
        """
        cfg = self.cfg or CTDHelpConfig()
        tool_name = cfg.tool_name or getattr(self.model, "name", None) or "Tool"
        usage = self._guess_usage(cfg)

        # terminal width
        width = shutil.get_terminal_size(fallback=(cfg.min_total_width, 24)).columns
        width = max(cfg.min_total_width, width)

        defaults = self._get_defaults_map()
        params = [p for p in self.model.list_parameters() if advanced or not self._is_advanced(p, cfg)]

        # group params by section
        sections: dict[str, list] = {k: [] for k in cfg.section_order}
        for p in params:
            sections.setdefault(self._section_for(p.name, cfg), []).append(p)

        # header
        print(f"{tool_name}")
        if cfg.show_description:
            desc = (getattr(self.model, "opt_attribs", {}) or {}).get("description", "")
            if isinstance(desc, str) and desc.strip():
                print(textwrap.fill(desc.strip(), width=width))
                print()
            docurl = (getattr(self.model, "opt_attribs", {}) or {}).get("docurl", None)
            if docurl:
                print(f"Docs: {docurl}\n")

        print("Usage:")
        print(usage + "\n")

        left_max = cfg.left_col_width
        pad = cfg.col_padding

        for sec in cfg.section_order:
            group = sections.get(sec, [])
            if not group:
                continue
            print(f"{sec}:")
            for p in sorted(group, key=lambda x: x.name):
                name = p.name
                typ = self._typename(p)
                default_val = self._param_default(p, defaults)
                req = getattr(p, "required", False)
                is_list = getattr(p, "is_list", False)

                # left column text
                if self._flag_hint(p, default_val):
                    left = f"-{name}"
                else:
                    many = " ..." if is_list else ""
                    left = f"-{name} <{typ}>{many}"

                # right column details
                bits: list[str] = []
                desc = getattr(p, "description", "") or ""
                if desc:
                    bits.append(desc)
                if default_val not in (None, ""):
                    bits.append(f"default: {default_val}")
                if req:
                    bits.append("required")
                ch = self._choices_text(p)
                if ch:
                    bits.append(f"choices: {ch}")
                fmts = self._formats_text(p)
                if fmts:
                    bits.append(f"formats: {fmts}")

                right = " — ".join(bits) if bits else ""
                left_render = f"  {left}"

                if len(left_render) >= left_max:
                    print(left_render)
                    if right:
                        wrapped = textwrap.fill(
                            right,
                            width=width,
                            initial_indent=" " * (left_max + pad),
                            subsequent_indent=" " * (left_max + pad),
                        )
                        print(wrapped)
                    else:
                        print()
                else:
                    right_width = max(10, width - (left_max + pad))
                    right_wrapped = textwrap.fill(
                        right,
                        width=right_width,
                        initial_indent="",
                        subsequent_indent="",
                    )
                    first_line, *rest = right_wrapped.splitlines() or [""]
                    print(left_render.ljust(left_max + pad) + first_line)
                    for line in rest:
                        print(" " * (left_max + pad) + line)
            print()

        if cfg.show_directives:
            self._print_ctd_directives()
        if cfg.show_tip_when_basic and not advanced:
            print("\nTip: use --help-advanced (or --helphelp) for all options, including advanced.")

    # ------------------------ internal helpers ------------------------

    @staticmethod
    def _typename(p) -> str:
        """Return a human-readable type label for a CTDopts parameter."""
        t = getattr(p, "type", None)
        if isinstance(t, type):  # int/float/bool
            return t.__name__
        if isinstance(t, str):   # "input-file", "output-file", "string", ...
            return t
        return "value"

    @staticmethod
    def _choices_text(p) -> str | None:
        """Format a parameter's allowed values (choices or restrictions)."""
        if getattr(p, "choices", None):
            return ", ".join(map(str, p.choices))
        restr = getattr(p, "restrictions", None)
        if isinstance(restr, (list, tuple)) and restr:
            return ", ".join(map(str, restr))
        return None

    @staticmethod
    def _formats_text(p) -> str | None:
        """Format allowed file extensions for file parameters."""
        fmts = getattr(p, "file_formats", None)
        if fmts:
            return ", ".join(fmts)
        return None

    @staticmethod
    def _flag_hint(p, default_val) -> bool:
        """
        True if the parameter should render like a boolean flag.

        Booleans that default to False (and are not required) are rendered as
        `-flag` with no value placeholder, mirroring typical CLI behavior.
        """
        t = getattr(p, "type", None)
        is_bool = (t is bool) or (isinstance(t, str) and t.lower() == "bool")
        if not is_bool or getattr(p, "required", False):
            return False
        if default_val is None:  # unknown -> treat as False (flag)
            return True
        s = str(default_val).strip().lower()
        return s in ("false", "0", "none", "off", "")

    def _param_default(self, p, defaults_map: Mapping[str, Any]) -> Any | None:
        """
        Resolve the default value for a parameter.

        Resolution order:
          1) CTD defaults (from `model.get_defaults()`)
          2) `p.default` (CTDopts runtime default)
          3) `p.default_value` (legacy naming)
        """
        if p.name in defaults_map:
            return defaults_map[p.name]
        if hasattr(p, "default"):
            return p.default
        if hasattr(p, "default_value"):
            return p.default_value
        return None

    def _get_defaults_map(self) -> dict[str, Any]:
        """Build a flat name→default map from the CTD model."""
        try:
            return self._flatten(self.model.get_defaults(), as_string=True)
        except Exception:
            return {}

    @staticmethod
    def _section_for(name: str, cfg: "CTDHelpConfig") -> str:
        """Determine which help section a parameter belongs to."""
        if name in cfg.io_param_names:
            return "I/O"
        if name.startswith("infer:ipf"):
            return "IPF"
        head = name.split(":", 1)[0]
        return cfg.section_map.get(head, "General")

    @staticmethod
    def _is_advanced(p, cfg: "CTDHelpConfig") -> bool:
        """Check if a parameter should be hidden in basic help."""
        if p.name in cfg.advanced_by_name:
            return True
        tags = getattr(p, "tags", None) or []
        return cfg.advanced_tag in set(tags)

    def _guess_usage(self, cfg: "CTDHelpConfig") -> str:
        """Render a friendly single-line Usage: string."""
        binname = cfg.binary_name or (sys.argv[0].split("/")[-1] or "tool.py")
        params = self.model.list_parameters()

        def is_input(p):
            t = str(getattr(p, "type", "")).lower()
            return (
                p.name in {"in", "input", "input_file", "input_files"}
                or t.startswith("input-")
                or t == "input-file"
            )

        def is_output(p):
            t = str(getattr(p, "type", "")).lower()
            return (
                p.name in {"out", "output", "output_file"}
                or t.startswith("output-")
                or t == "output-file"
            )

        ins = [p for p in params if is_input(p)]
        outs = [p for p in params if is_output(p)]

        if ins:
            p_in = ins[0]
            in_typ = f"<{self._typename(p_in)}{'...' if getattr(p_in, 'is_list', False) else ''}>"
            out_part = ""
            if outs:
                p_out = outs[0]
                out_typ = f"<{self._typename(p_out)}{'...' if getattr(p_out, 'is_list', False) else ''}>"
                out_opt = "" if getattr(p_out, "required", False) else "["
                out_end = "" if getattr(p_out, "required", False) else "]"
                out_part = f" {out_opt}-out {out_typ}{out_end}"
            return f"  {binname} -in {in_typ}{out_part} [options]"
        return f"  {binname} [options]"

    @staticmethod
    def _print_ctd_directives() -> None:
        """Print OpenMS CTD/INI directive flags."""
        print("CTD / INI directives (from OpenMS):")
        print("  -ini <FILE>                Load parameters from CTD")
        print("  -write_ini <FILE>          Write tool CTD (defaults)")
        print("  -write_param_ctd <FILE>    Write CTD with current parameters")

    # --------------------- tiny local flatten helper ---------------------

    @staticmethod
    def _flatten(arg_dict: Mapping[str, Any], *, as_string: bool = False) -> dict:
        """Flatten nested dict to {':'.join(path): value} (or tuple path)."""
        result: dict = {}

        def rec(sub, level: list[str]):
            for k, v in (sub or {}).items():
                if isinstance(v, Mapping):
                    rec(v, level + [k])
                else:
                    result[tuple(level + [k])] = v

        rec(arg_dict or {}, [])
        if as_string:
            return {":".join(k): v for k, v in result.items()}
        return result


# code related to CTD support
def addParamToCTDopts(defaults, model):
    keys = defaults.keys()
    for key in keys:
        ctd_tags = defaults.getTags(key)
        value = defaults.getValue(key)
        desc = defaults.getDescription(key)

        ctd_required = False
        if "required" in ctd_tags:
            ctd_required = True

        ctd_type = None
        ctd_type_str = ""
        ctd_list = False

        if isinstance(value, int):
            ctd_type = int
            ctd_type_str = "int"
        elif isinstance(value, float):
            ctd_type = float
            ctd_type_str = "float"
        elif isinstance(value, str):
            ctd_type = str
            ctd_type_str = "str"
        elif isinstance(value, list):
            ctd_list = True
            # we can't determine the type based on an element so we need this helper function
            value_type = defaults.getValueType(key)
            if value_type == pms.ValueType.STRING_LIST:
                ctd_type = str
                ctd_type_str = "list of str"
            elif value_type == pms.ValueType.DOUBLE_LIST:
                ctd_type = float
                ctd_type_str = "list of float"
            elif value_type == pms.ValueType.INT_LIST:
                ctd_type = int
                ctd_type_str = "list of int"

        print(
            'Adding parameter: {0} with value: {1} and description: "{2}".'.format(
                key, value, desc
            )
        )
        print(
            "        required: {0} \t tags: {1} \t type: {2}.".format(
                ctd_required, ctd_tags, ctd_type_str
            )
        )

        model.add(
            key.decode(),
            required=ctd_required,
            type=ctd_type,
            default=value,
            is_list=ctd_list,
            description=desc,
        )

class PyTOPPArgParser:
    """
    OpenMS-style CLI/INI parser with CTDopts validation.

    Key semantics (identical in both modes):
        defaults  <-  INI values  <-  **explicit** CLI flags

    Modes:
        - "pure"      : Does not use pyopenms tool algo's, i.e. no OpenMS Param object. Returns a validated flat dict.
        - "pyopenms"  : Same merge+validate, and (optionally) synchronize an
                        OpenMS `Param` by writing/reading a CTD so your
                        algorithm layer sees the same values.

    Unknown flags:
        Fails early with a helpful message (close-match hint) instead of
        silently falling back to defaults.

    Forwarding:
        Tokens beginning with `+--` are ignored by validation—they're meant to
        be passed directly to the tool (e.g., pyprophet score --pi0_lambda 0 0 0). This is mainly for tools that have a lot of additional/advanced options that we do not necessarily want to add to the main CTD model.

    Usage:
        parser = PyTOPPArgParser(model)
        values = parser.parse(sys.argv[1:], mode="pure")

        # if you have an OpenMS Param you want kept in sync:
        import pyopenms as poms
        om_param = poms.Param()
        values = parser.parse(sys.argv[1:], mode="pyopenms", openms_param=om_param)
    """

    def __init__(self, model: CTDModel):
        self.model = model

    # ---------- public API -------------------------------------------------

    def parse(
        self,
        argv: list[str],
        *,
        mode: str = "pure",
        openms_param: Optional[Any] = None,
    ) -> dict[str, Any]:
        """
        Dispatch to the chosen parser.

        Args:
            argv: Raw tokens (usually `sys.argv[1:]`).
            mode: "pure" (no OpenMS Param store).
            openms_param: Optional pyOpenMS `Param` to synchronize in "pyopenms" mode.

        Returns:
            Validated mapping with colon-joined keys (e.g., "scoring:level").
        """
        if mode == "pure":
            return self._parse_pure(argv)
        if mode == "pyopenms":
            return self._parse_pyopenms(argv, openms_param=openms_param)
        raise ValueError("mode must be 'pure' or 'pyopenms'")

    # ---------- pure (non-pyOpenMS) path ----------------------------------

    def _parse_pure(self, argv: list[str]) -> dict[str, Any]:
        argv = self._strip_prog(argv)

        directives, remaining = self._split_directives(argv)
        self._fail_on_unknowns(remaining)

        # handle write-tool-ctd early
        if directives["write_tool_ctd"] is not None:
            out = self._default_or_path(directives["write_tool_ctd"], f"{self._tool_name()}.ctd")
            self.model.write_ctd(out)
            print(f"Wrote tool CTD: {out}")
            sys.exit(0)

        # base from -ini (if any)
        base_args: dict[str, Any] = {}
        if directives["input_ctd"] not in (None, True):
            base_args = args_from_file(directives["input_ctd"])

        # parse CLI and keep only flags the user explicitly typed
        cli_full = self.model.parse_cl_args(remaining, prefix="-")
        explicit = self._explicit_flags(remaining)
        cli_explicit = self._filter_cli_defaults(cli_full, explicit)

        # deep merge: defaults <- ini <- explicit CLI
        merged = self._merge_flat(self.model.get_defaults(), base_args, cli_explicit)

        # validate and maybe write param CTD
        validated = self.model.validate_args(merged)

        if directives["write_param_ctd"] is not None:
            out = self._default_or_path(
                directives["write_param_ctd"], f"{self._tool_name()}_params.ctd"
            )
            self.model.write_ctd(out, validated)
            print(f"Wrote parameter CTD: {out}")
            sys.exit(0)

        return validated

    # ---------- pyOpenMS path (optional) -----------------------------------

    def _parse_pyopenms(self, argv: list[str], *, openms_param: Optional[Any]) -> dict[str, Any]:
        """
        Same semantics as `_parse_pure`, but if `openms_param` is supplied, it is
        synchronized with the final *validated* values via a round-trip CTD.

        This keeps your pyOpenMS algorithm layer aligned with what the CLI/INI
        produced after merging.
        """
        argv = self._strip_prog(argv)

        # We still use our exact directive split for consistency, but accept
        # -ini/-write_ini/-write_param_ctd as in OpenMS tools.
        directives, remaining = self._split_directives(argv)
        self._fail_on_unknowns(remaining)

        if directives["write_tool_ctd"] is not None:
            out = self._default_or_path(directives["write_tool_ctd"], f"{self._tool_name()}.ctd")
            self.model.write_ctd(out)
            print(f"Wrote tool CTD: {out}")
            sys.exit(0)

        base_args: dict[str, Any] = {}
        if directives["input_ctd"] not in (None, True):
            base_args = args_from_file(directives["input_ctd"])

        cli_full = self.model.parse_cl_args(remaining, prefix="-")
        explicit = self._explicit_flags(remaining)
        cli_explicit = self._filter_cli_defaults(cli_full, explicit)

        merged = self._merge_flat(self.model.get_defaults(), base_args, cli_explicit)
        validated = self.model.validate_args(merged)

        # Write param CTD if requested
        if directives["write_param_ctd"] is not None:
            out = self._default_or_path(
                directives["write_param_ctd"], f"{self._tool_name()}_params.ctd"
            )
            self.model.write_ctd(out, validated)
            print(f"Wrote parameter CTD: {out}")
            sys.exit(0)

        # Synchronize a provided pyOpenMS Param with the final values
        if openms_param is not None:
            try:
                import pyopenms as pms  # local import to keep "pure" mode light
                tmp = tempfile.NamedTemporaryFile(suffix=".ini", delete=False)
                tmp.close()
                try:
                    self.model.write_ctd(tmp.name, validated)
                    param = pms.Param()
                    fh = pms.ParamXMLFile()
                    fh.load(tmp.name, param)
                    # overwrite all (True) so it mirrors the final values
                    openms_param.update(param, True)
                finally:
                    os.unlink(tmp.name)
            except Exception as e:
                # Don't crash the parse if pyOpenMS sync fails; surface a clear note.
                print(
                    textwrap.dedent(
                        f"""\
                        [warning] Unable to synchronize openms_param from validated values: {e}
                        Your parsed values are still returned; pyOpenMS Param may be out of sync."""
                    ),
                    file=sys.stderr,
                )

        return validated

    # ---------- helpers: directives, flags, merging ------------------------

    def _tool_name(self) -> str:
        return getattr(self.model, "name", None) or "tool"

    @staticmethod
    def _strip_prog(argv: list[str]) -> list[str]:
        if not argv:
            return []
        # if someone passed full sys.argv, drop the program name
        return argv[1:] if argv and argv[0] == sys.argv[0] else argv

    @staticmethod
    def _split_directives(argv: list[str]) -> tuple[dict, list[str]]:
        """
        Exact-match extraction of OpenMS directives:
          -ini [FILE]
          -write_ini [FILE]
          -write_param_ctd [FILE]
        """
        directives = {"input_ctd": None, "write_tool_ctd": None, "write_param_ctd": None}
        rest: list[str] = []
        i = 0
        while i < len(argv):
            t = argv[i]
            if t in ("-ini", "--ini"):
                val: Any = True
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    val = argv[i + 1]
                    i += 1
                directives["input_ctd"] = val
            elif t in ("-write_ini", "--write_ini"):
                val = True
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    val = argv[i + 1]
                    i += 1
                directives["write_tool_ctd"] = val
            elif t in ("-write_param_ctd", "--write_param_ctd"):
                val = True
                if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                    val = argv[i + 1]
                    i += 1
                directives["write_param_ctd"] = val
            else:
                rest.append(t)
            i += 1
        return directives, rest

    def _fail_on_unknowns(self, tokens: list[str]) -> None:
        unknowns = self._find_unknown_flags(tokens)
        if not unknowns:
            return
        names = [p.name for p in self.model.list_parameters()]
        lines = []
        for u in unknowns:
            cand = get_close_matches(u.lstrip("-"), names, n=1)
            hint = f"  (did you mean '-{cand[0]}'?)" if cand else ""
            lines.append(f"  {u}{hint}")
        print("Unrecognized option(s):\n" + "\n".join(lines), file=sys.stderr)
        sys.exit(2)

    def _find_unknown_flags(self, tokens: list[str]) -> list[str]:
        known = {f"-{p.name}" for p in self.model.list_parameters()}
        known.update({"-ini", "-write_ini", "-write_param_ctd"})
        seen, out = set(), []
        for tok in tokens:
            if not tok.startswith("-"):
                continue
            if tok.startswith("+--"):               # forwarded to sub-tools
                continue
            if len(tok) > 1 and tok[1].isdigit():  # numeric literals like -0.1
                continue
            flag = tok.split("=", 1)[0]
            if flag in known or flag in seen:
                continue
            seen.add(flag)
            out.append(flag)
        return out

    @staticmethod
    def _explicit_flags(tokens: list[str]) -> set[str]:
        """Return set of colon-keys (without leading '-') explicitly present in tokens."""
        exp: set[str] = set()
        for tok in tokens:
            if not tok.startswith("-"):
                continue
            if tok in ("-ini", "-write_ini", "-write_param_ctd"):
                continue
            if tok.startswith("+--"):                # forwarded, not a model flag
                continue
            if len(tok) > 1 and tok[1].isdigit():   # negative numbers like -0.1
                continue
            flag = tok.split("=", 1)[0].lstrip("-")
            exp.add(flag)
        return exp

    def _filter_cli_defaults(self, cli_args: dict, explicit: set[str]) -> dict[str, Any]:
        """
        CTDopts returns a dict pre-populated with defaults. Keep only the keys
        the user actually typed, and drop `_Null` sentinels.
        """
        flat = self._flatten(cli_args, as_string=True)
        out: dict[str, Any] = {}
        for k, v in flat.items():
            if k in explicit and not isinstance(v, _Null):
                out[k] = v
        return out

    @staticmethod
    def _default_or_path(val: Any, fallback: str) -> str:
        return fallback if val is True else str(val)

    # ---------- tiny flatten/deep-merge helpers (self-contained) -----------

    @staticmethod
    def _flatten(arg_dict: Mapping[str, Any], as_string: bool = False) -> dict:
        """Flatten nested dict to {tuple(path): value} or {':'.join(path): value}."""
        result: dict = {}

        def rec(sub, level: list[str]):
            for k, v in (sub or {}).items():
                if isinstance(v, Mapping):
                    rec(v, level + [k])
                else:
                    result[tuple(level + [k])] = v

        rec(arg_dict or {}, [])
        if as_string:
            return {":".join(k): v for k, v in result.items()}
        return result

    @staticmethod
    def _set_nested(d: dict, key_path: Iterable[str], value: Any) -> None:
        key_path = [key_path] if isinstance(key_path, str) else list(key_path)
        cur = d
        for k in key_path[:-1]:
            if k not in cur or not isinstance(cur[k], dict):
                cur[k] = {}
            cur = cur[k]
        cur[key_path[-1]] = value

    @classmethod
    def _as_nested(cls, d: Mapping[str, Any]) -> dict:
        """Convert flat colon-keys to nested dicts (recursively)."""
        if not isinstance(d, Mapping):
            return {}
        out: dict = {}
        for k, v in d.items():
            if isinstance(v, Mapping):
                v = cls._as_nested(v)
            if isinstance(k, str) and ":" in k:
                cls._set_nested(out, k.split(":"), v)
            else:
                if k in out and isinstance(out[k], Mapping) and isinstance(v, Mapping):
                    out[k] = cls._deep_merge(out[k], v)
                else:
                    out[k] = v
        return out

    @classmethod
    def _deep_merge(cls, a: Mapping[str, Any], b: Mapping[str, Any]) -> dict:
        """Deep merge where values in b override a."""
        res = dict(a)
        for k, v in (b or {}).items():
            if k in res and isinstance(res[k], Mapping) and isinstance(v, Mapping):
                res[k] = cls._deep_merge(res[k], v)
            else:
                res[k] = v
        return res

    @classmethod
    def _merge_flat(cls, *dicts_like: Mapping[str, Any]) -> dict[str, Any]:
        """
        Normalize inputs to nested, deep-merge (later wins), then flatten to
        colon-joined keys for CTDopts validation.
        """
        nested: dict = {}
        for d in dicts_like:
            nested = cls._deep_merge(nested, cls._as_nested(d or {}))
        return cls._flatten(nested, as_string=True)