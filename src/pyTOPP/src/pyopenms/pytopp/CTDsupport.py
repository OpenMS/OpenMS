from __future__ import annotations

import os
import shutil
import sys
import tempfile
import textwrap
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import CTDopts
import pyopenms as pms

# add near the top with the other CTDopts imports
from CTDopts.CTDopts import (
    CTDModel,
    args_from_file,
    flatten_dict,
    override_args,
    parse_cl_directives,
)


@dataclass
class HelpConfig:
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
    binary_name: str | None = None  # e.g. "PyProphet.py" (default: sys.argv[0])
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


def _typename(p) -> str:
    """Return a human-readable type label for a CTDopts parameter.

    Args:
        p: A CTDopts Parameter-like object with a `type` attribute.

    Returns:
        A string such as `"int"`, `"float"`, `"bool"`, or a CTD string type
        like `"input-file"`, `"output-file"`, `"string"`. Falls back to
        `"value"` if not recognizable.
    """
    t = getattr(p, "type", None)
    if isinstance(t, type):  # int/float/bool
        return t.__name__
    if isinstance(t, str):  # "input-file", "output-file", "string", ...
        return t
    return "value"


def _choices_text(p) -> str | None:
    """Format a parameter's allowed values (choices or restrictions).

    Args:
        p: A CTDopts Parameter-like object that may specify `choices` or
           `restrictions`.

    Returns:
        A comma-separated string of allowed values, or `None` if unrestricted.
    """
    if getattr(p, "choices", None):
        return ", ".join(map(str, p.choices))
    restr = getattr(p, "restrictions", None)
    if isinstance(restr, (list, tuple)) and restr:
        return ", ".join(map(str, restr))
    return None


def _formats_text(p) -> str | None:
    """Format allowed file extensions for file parameters.

    Args:
        p: A CTDopts Parameter-like object with optional `file_formats`.

    Returns:
        A comma-separated string of file formats (extensions) or `None` if not
        a file parameter or unspecified.
    """
    fmts = getattr(p, "file_formats", None)
    if fmts:
        return ", ".join(fmts)
    return None


def _param_default(p, defaults_map):
    """Resolve the default value for a parameter.

    Resolution order:
      1) CTD defaults (from `model.get_defaults()`)
      2) `p.default` (CTDopts runtime default)
      3) `p.default_value` (legacy naming)

    Args:
        p: A CTDopts Parameter-like object.
        defaults_map: A flat `dict[str, Any]` mapping parameter names to default
            values, typically produced by `flatten_dict(model.get_defaults(), as_string=True)`.

    Returns:
        The default value, or `None` if no default is known.
    """
    if p.name in defaults_map:
        return defaults_map[p.name]
    if hasattr(p, "default"):
        return p.default
    if hasattr(p, "default_value"):
        return p.default_value
    return None


def _flag_hint(p, default_val) -> bool:
    """Return True if the parameter should render like a boolean flag.

    Booleans that default to *False* (and are not required) are rendered as
    `-flag` with no value placeholder, mirroring typical CLI behavior.

    Args:
        p: CTDopts Parameter-like object.
        default_val: The resolved default value from `_param_default`.

    Returns:
        True if the parameter should be shown as a standalone flag, else False.
    """
    t = getattr(p, "type", None)
    is_bool = (t is bool) or (isinstance(t, str) and t.lower() == "bool")
    if not is_bool or getattr(p, "required", False):
        return False
    if default_val is None:  # unknown -> treat as False (flag)
        return True
    s = str(default_val).strip().lower()
    return s in ("false", "0", "none", "off", "")


def _get_defaults_map(model) -> dict[str, str]:
    """Build a flat name→default map from a CTD model.

    Args:
        model: A `CTDModel` instance.

    Returns:
        A `dict` mapping parameter names to their default values as strings.
        Returns an empty dict if defaults cannot be obtained.
    """
    try:
        return flatten_dict(model.get_defaults(), as_string=True)
    except Exception:
        return {}


def _section_for(name: str, cfg: HelpConfig) -> str:
    """Determine which help section a parameter belongs to.

    Args:
        name: Full parameter name as defined in the CTD model.
        cfg: HelpConfig controlling section mapping.

    Returns:
        The section title (e.g., `"I/O"`, `"Scoring"`, `"General"`, etc.).
    """
    if name in cfg.io_param_names:
        return "I/O"
    if name.startswith("infer:ipf"):
        return "IPF"
    head = name.split(":", 1)[0]
    return cfg.section_map.get(head, "General")


def _is_advanced(p, cfg: HelpConfig) -> bool:
    """Check if a parameter should be hidden in basic help.

    A parameter is considered advanced if:
      * Its name is present in `cfg.advanced_by_name`, or
      * It has a tag matching `cfg.advanced_tag`.

    Args:
        p: CTDopts Parameter-like object.
        cfg: HelpConfig with advanced-selection settings.

    Returns:
        True if advanced; False otherwise.
    """
    if p.name in cfg.advanced_by_name:
        return True
    tags = getattr(p, "tags", None) or []
    return cfg.advanced_tag in set(tags)


def _guess_usage(model, cfg: HelpConfig) -> str:
    """Render a friendly single-line Usage: string.

    Heuristic:
      * If we can identify a primary input parameter and (optionally) an output
        parameter, show something like:
          `<bin> -in <_InFile...> [-out <_OutFile>] [options]`
      * Otherwise, fall back to:
          `<bin> [options]`

    Args:
        model: CTD model to inspect.
        cfg: HelpConfig providing `binary_name`.

    Returns:
        A single Usage line (without the `'Usage:'` label).
    """
    binname = cfg.binary_name or (sys.argv[0].split("/")[-1] or "tool.py")
    params = model.list_parameters()

    # heuristic: pick first required input and (optional) output
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
        in_typ = (
            f"<{_typename(p_in)}{'...' if getattr(p_in, 'is_list', False) else ''}>"
        )
        out_part = ""
        if outs:
            p_out = outs[0]
            out_typ = f"<{_typename(p_out)}{'...' if getattr(p_out, 'is_list', False) else ''}>"
            out_opt = "" if getattr(p_out, "required", False) else "["
            out_end = "" if getattr(p_out, "required", False) else "]"
            out_part = f" {out_opt}-out {out_typ}{out_end}"
        return f"  {binname} -in {in_typ}{out_part} [options]"
    return f"  {binname} [options]"


def _print_ctd_directives():
    """Print OpenMS CTD/INI directive flags.

    This mirrors the flags commonly exposed by pyTOPP wrappers:
      - `-ini` to load parameters from a CTD,
      - `-write_ini` to write the tool's default CTD,
      - `-write_param_ctd` to write a CTD with current parameter values.
    """
    print("CTD / INI directives (from OpenMS):")
    print("  -ini <FILE>                Load parameters from CTD")
    print("  -write_ini <FILE>          Write tool CTD (defaults)")
    print("  -write_param_ctd <FILE>    Write CTD with current parameters")


def print_model_help(model, *, advanced: bool, cfg: HelpConfig | None = None) -> None:
    """Pretty, two-column help printer for CTDopts `CTDModel`.

    The output is grouped into sections (for example: I/O, General, Scoring, Inference, Export), wraps descriptions to the terminal width, and shows defaults, choices, formats, and required flags. When `advanced=False`, parameters tagged as advanced (or explicitly listed in the config) are hidden.

    Args:
        model:
            The CTDopts `CTDModel` to render.
        advanced:
            If `True`, show all parameters. If `False`, hide advanced ones.
        cfg:
            Optional `HelpConfig` to control formatting and selection. If not
            provided, a default `HelpConfig()` is used.

    Returns:
        None. Prints directly to stdout.
    """
    cfg = cfg or HelpConfig()
    tool_name = cfg.tool_name or getattr(model, "name", None) or "Tool"
    usage = _guess_usage(model, cfg)

    # terminal width
    width = shutil.get_terminal_size(fallback=(cfg.min_total_width, 24)).columns
    width = max(cfg.min_total_width, width)

    defaults = _get_defaults_map(model)
    params = [
        p for p in model.list_parameters() if advanced or not _is_advanced(p, cfg)
    ]

    # group params by section
    sections: dict[str, list] = {k: [] for k in cfg.section_order}
    for p in params:
        sections.setdefault(_section_for(p.name, cfg), []).append(p)

    # header
    print(f"{tool_name}")
    if cfg.show_description:
        desc = (model.opt_attribs.get("description", None)).strip()
        if desc and isinstance(desc, str):
            print(textwrap.fill(desc, width=width))
            print()

        docurl = model.opt_attribs.get("docurl", None)
        if cfg.show_description and docurl:
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
            typ = _typename(p)
            default_val = _param_default(p, defaults)
            req = getattr(p, "required", False)
            is_list = getattr(p, "is_list", False)

            # left column text
            if _flag_hint(p, default_val):
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
            ch = _choices_text(p)
            if ch:
                bits.append(f"choices: {ch}")
            fmts = _formats_text(p)
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
        _print_ctd_directives()
    if cfg.show_tip_when_basic and not advanced:
        print(
            "\nTip: use --help-advanced (or --helphelp) for all options, including advanced."
        )

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


def parseCTDCommandLine(argv, model, openms_param):
    # Configure CTDOpt to use OpenMS style on the command line.
    directives = parse_cl_directives(
        argv, input_ctd="ini", write_tool_ctd="write_ini", prefix="-"
    )

    if (
        directives["write_tool_ctd"] is not None
    ):  # triggered if -write_ini was provided on CML
        # if called with -write_ini write CTD
        model.write_ctd(directives["write_tool_ctd"])
        exit(0)
    elif directives["input_ctd"] is not None:  # read ctd/ini file
        model = CTDModel(from_file=directives["input_ctd"])
        #        print(model.get_defaults())

        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(directives["input_ctd"], param)
        openms_param.update(param, True)
        return model.get_defaults(), openms_param

    else:  # only command line options provided
        temp = tempfile.NamedTemporaryFile(
            suffix="ini"
        )  # makes sure we get a writable file
        tmp_name = temp.name
        temp.close()  # removes the file

        model.write_ctd(tmp_name)
        param = pms.Param()
        fh = pms.ParamXMLFile()
        fh.load(tmp_name, param)
        openms_param.update(param)
        os.remove(tmp_name)
        return model.parse_cl_args(argv), openms_param

def _split_openms_directives(argv):
    """
    Exact-match extraction of OpenMS-style directives from argv (no script name).
    Recognizes:
      -ini [FILE]
      -write_ini [FILE]
      -write_param_ctd [FILE]
    Returns: (directives_dict, remaining_argv)
    """
    directives = {"input_ctd": None, "write_tool_ctd": None, "write_param_ctd": None}
    rest = []
    i = 0
    while i < len(argv):
        t = argv[i]
        # match EXACTLY these flags; don't confuse -in with -ini
        if t in ("-ini", "--ini"):
            val = True
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

def parseCTDCommandLinePure(argv, model):
    """
    Parse OpenMS-style CLI for a CTDopts model for non-pyOpenMS tools (i.e. the tool does not internally use the OpenMS Param object).

    This implements the canonical merge semantics used by pyTOPP tools:

        defaults  <-  values from -ini CTD  <-  explicit CLI flags

    It returns a single **validated** dictionary of effective parameters.

    Supported directives (all single-dash, OpenMS style):
      - `-ini <FILE>`:
          Read a CTD/INI file and treat its values as the base.
      - `-write_ini [<FILE>]`:
          Write the tool’s schema/default CTD and exit.
          If `<FILE>` is omitted, `<model.name>.ctd` is used.
      - `-write_param_ctd [<FILE>]`:
          Write a CTD containing the **current effective values** (after merging
          CLI and `-ini`) and exit. If `<FILE>` is omitted,
          `<model.name>_params.ctd` is used.

    Notes:
      * If `argv` includes the script name as `argv[0]`, it is removed.
      * Only OpenMS-style options (single-dash, e.g. `-in`, `-out`, `-threads`)
        are parsed; GNU-style long options are not recognized here.
      * This function may call `sys.exit(0)` when a write directive is used.

    Args:
        argv (list[str]): Raw command-line tokens (typically `sys.argv` or
            `sys.argv[1:]`).
        model (CTDModel): A CTDopts model that provides:
            - `get_defaults() -> dict`
            - `parse_cl_args(tokens, prefix='-') -> dict`
            - `validate_args(mapping) -> dict`
            - `write_ctd(path, values: dict | None = None) -> None`

    Returns:
        dict[str, Any]: The **validated** mapping of parameter names to values.

    Side effects:
        May write CTD files and terminate the process with `sys.exit(0)` when
        `-write_ini` or `-write_param_ctd` is supplied.
    """
    # ensure we don't pass the script name
    if argv and argv[0] == sys.argv[0]:
        argv = argv[1:]

    directives, remaining = _split_openms_directives(argv)

    # Write tool CTD (schema/defaults) and exit
    if directives["write_tool_ctd"] is not None:
        out = directives["write_tool_ctd"]
        if out is True:
            out = f"{getattr(model, 'name', 'tool')}.ctd"
        model.write_ctd(out)
        print(f"Wrote tool CTD: {out}")
        sys.exit(0)

    # Load args from -ini (CTD) if provided
    base_args = {}
    if directives["input_ctd"] not in (None, True):
        base_args = args_from_file(directives["input_ctd"])

    # Parse remaining CLI with CTDopts (prefix '-')
    cli_args = model.parse_cl_args(remaining, prefix="-")

    merged = override_args(model.get_defaults(), base_args, cli_args)

    validated = model.validate_args(merged)

    # Write CTD with current values if requested
    if directives["write_param_ctd"] is not None:
        out = directives["write_param_ctd"]
        if out is True:
            out = f"{getattr(model, 'name', 'tool')}_params.ctd"
        model.write_ctd(out, validated)
        print(f"Wrote parameter CTD: {out}")
        sys.exit(0)

    return validated

