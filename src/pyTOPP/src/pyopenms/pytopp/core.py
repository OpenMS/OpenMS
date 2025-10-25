from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List

_TRUE_TOKENS  = {"true", "1", "yes", "on"}
_FALSE_TOKENS = {"false", "0", "no", "off"}


# ------------------- generic value helpers -------------------
def as_bool(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    return str(v).strip().lower() in _TRUE_TOKENS


def is_nullish(v) -> bool:
    # CTDopts uses a private sentinel class named "_Null" for “no value”
    return v is None or type(v).__name__ in {"_Null", "Null"} or v == ""


def first_value_for_flag(argv: list[str], name: str) -> str | None:
    """
    Return raw token following -<name>, or value from -<name>=VALUE if used.
    Returns "" if -<name> has no following token; returns None if not present.
    """
    flag = f"-{name}"
    for i, t in enumerate(argv):
        if t == flag:
            if i + 1 < len(argv):
                return argv[i + 1]
            return ""
        if t.startswith(flag + "="):
            return t.split("=", 1)[1]
    return None


def tok_extra(value) -> list[str]:
    """
    Split shell-like extra args; undo '+--', '+-' and '\\-' escapes.
    Accepts None/_Null/etc. and returns [].
    """
    if is_nullish(value) or not isinstance(value, str) or not value.strip():
        return []
    s = value.strip()
    if s.startswith("-- "):  # optional nicety to allow '-- ' sentinel
        s = s[3:]
    toks = shlex.split(s)
    out: list[str] = []
    for t in toks:
        if t.startswith("+--") or t.startswith("+-") or t.startswith("\\-"):
            out.append(t[1:])
        else:
            out.append(t)
    return out


# ------------------- argv reconciliation helpers -------------------
def is_bool_param(p) -> bool:
    t = getattr(p, "type", None)
    return (t is bool) or (isinstance(t, str) and t.lower() == "bool")


def coerce_bool_from_argv(argv: list[str], name: str, current: bool) -> bool:
    """
    If user gave -<name> true/false (or =true/=false), or -no-<name>, override.
    Otherwise keep current.
    """
    flag   = f"-{name}"
    noflag = f"-no-{name}"

    for i, tok in enumerate(argv):
        if tok == noflag:
            return False
        if tok == flag:
            if i + 1 < len(argv) and not argv[i + 1].startswith("-"):
                v = argv[i + 1].strip().lower()
                if v in _TRUE_TOKENS:  return True
                if v in _FALSE_TOKENS: return False
            # bare -flag -> treat as True
            continue
        if tok.startswith(flag + "="):
            v = tok.split("=", 1)[1].strip().lower()
            if v in _TRUE_TOKENS:  return True
            if v in _FALSE_TOKENS: return False
    return current


def recover_missing_extras(arg_dict: dict, argv: list[str], keys: list[str]) -> None:
    """
    Some parsers can drop empty strings for :extra params. If a key is missing/nullish
    but present in argv, restore its raw text so tok_extra() can split it.
    """
    for k in keys:
        if is_nullish(arg_dict.get(k)):
            raw = first_value_for_flag(argv, k)
            if raw is not None:
                arg_dict[k] = raw


def coerce_all_bools_from_argv(model, arg_dict: dict, argv: list[str]) -> None:
    for p in model.list_parameters():
        if is_bool_param(p):
            cur = as_bool(arg_dict.get(p.name, getattr(p, "default", False)))
            arg_dict[p.name] = coerce_bool_from_argv(argv, p.name, cur)


# ------------------- subprocess runner -------------------
def run_tool(executable: str,
             args: Iterable[str],
             *,
             env=None,
             dry: bool = False,
             module_name: str | None = None) -> None:
    """
    Run an external tool and echo a friendly CMD line.

    If module_name is provided and executable equals that module name,
    run it as "python -m <module_name>" to guarantee the current venv.
    """
    if module_name and executable == module_name:
        display = f"{sys.executable} -m {module_name}"
        cmd = [sys.executable, "-m", module_name, *args]
    else:
        display = executable
        cmd = [executable, *args]

    print("CMD:", display, *args)
    if dry:
        return
    p = subprocess.run(cmd, env=env)
    if p.returncode != 0:
        raise RuntimeError(f"{executable} failed with exit code {p.returncode}")
