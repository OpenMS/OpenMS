from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Mapping, Optional

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


# Variables that are risky to inherit when spawning child processes
_BAD_ENV_VARS = {
    "PYTHONPATH", "PYTHONHOME", "PYTHONWARNINGS", "PYTHONUSERBASE",
    "LD_PRELOAD", "LD_AUDIT", "LD_LIBRARY_PATH",
    "DYLD_INSERT_LIBRARIES", "DYLD_LIBRARY_PATH",
}

# conservative allowlist we keep if the caller asks for a “minimal” env
_MINIMAL_KEEP = {
    "PATH", "HOME", "TMPDIR", "TEMP", "TMP",
    "SystemRoot", "WINDIR", "COMSPEC", "USERPROFILE", "USERNAME", "LANG", "LC_ALL",
}

_MODNAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*$")

def _sanitize_env(env: Optional[Mapping[str, str]]) -> Mapping[str, str]:
    """
    Create a safe environment for the subprocess:
    - If env is None, start from os.environ and strip risky keys.
    - If env is provided, copy it and strip risky keys.
    """
    base = dict(os.environ) if env is None else dict(env)
    for bad in _BAD_ENV_VARS:
        base.pop(bad, None)
    return base

def _quote_cmd(parts: Iterable[str]) -> str:
    """For logging only; shows a safely-quoted command line."""
    return shlex.join(map(str, parts))

def run_tool(
    executable: str,
    args: Iterable[str],
    *,
    env: Optional[Mapping[str, str]] = None,
    dry: bool = False,
    module_name: Optional[str] = None,
    timeout: Optional[float] = None,
    minimal_env: bool = False,
) -> None:
    """
    Run an external tool securely.

    Security notes:
      - Never uses the shell (shell=False).
      - Validates and resolves the executable path (unless using -m).
      - Sanitizes the environment (removes LD_PRELOAD, PYTHONPATH, etc).
      - Optionally enforces a minimal env (PATH/HOME/TMP* and a few OS vars).
      - Prints a quoted command preview (safe for copy/paste).
      - Supports a timeout to avoid hanging indefinitely.

    If module_name is provided and equals `executable`, the tool is invoked as
    `sys.executable -m <module_name>` to guarantee the current interpreter.
    """
    # Validate args are strings/path-like and free of NULs
    safe_args: list[str] = []
    for a in args:
        s = os.fspath(a) if hasattr(a, "__fspath__") else str(a)
        if "\x00" in s:
            raise ValueError("Argument contains NUL byte")
        safe_args.append(s)

    # Build command
    if module_name:
        if module_name != executable:
            raise ValueError("module_name must match executable when provided")
        if not _MODNAME_RE.match(module_name):
            raise ValueError(f"Invalid module name: {module_name!r}")
        # Ensure the module exists (fail fast with a clear error)
        try:
            import importlib.util as _iu
            if _iu.find_spec(module_name) is None:
                # Might still be available at runtime via site-packages; don’t block if caller
                # explicitly wants to try. Comment the next line out if you prefer soft-check.
                pass
        except Exception:
            pass
        cmd = [sys.executable, "-m", module_name, *safe_args]
        display = cmd  # accurate preview
    else:
        # Resolve the executable to an absolute path via PATH
        exe_path = shutil.which(executable)
        if exe_path is None:
            raise FileNotFoundError(f"Executable not found on PATH: {executable!r}")
        cmd = [exe_path, *safe_args]
        display = cmd

    # Prepare environment
    if minimal_env:
        # Start fresh, keep only a small allowlist from the current env
        kept = {k: v for k, v in os.environ.items() if k in _MINIMAL_KEEP}
        if env:
            kept.update(env)
        env_final = _sanitize_env(kept)
    else:
        env_final = _sanitize_env(env)

    # Friendly preview (no execution)
    print("CMD:", _quote_cmd(display))
    if dry:
        return

    try:
        # Close file descriptors (POSIX default True), never use shell
        subprocess.run(
            cmd,
            env=env_final,
            check=True,
            timeout=timeout,
        )
    except FileNotFoundError as e:
        raise RuntimeError(f"Failed to execute: {display[0]!r}: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f"Command timed out after {timeout}s: {_quote_cmd(display)}") from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Command failed with exit code {e.returncode}: {_quote_cmd(display)}"
        ) from e

