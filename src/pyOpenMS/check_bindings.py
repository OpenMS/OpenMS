#!/usr/bin/env python3
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Verify that every expected pyOpenMS extension module imports.

Run before ``nanobind.stubgen``. Two reasons it exists:

1. ``pyopenms/__init__.py`` deliberately tolerates a partial import: it warns per
   module and only raises when *no* binding imports at all. stubgen would happily
   generate stubs from whatever did load, shipping a ``.pyi`` set that silently
   describes half an API. This script makes a missing module fail the build.

2. When a module fails to load on Windows, the loader says only
   ``ImportError: DLL load failed while importing _pyopenms_analysis: The
   specified module could not be found`` -- never which dependency is missing.
   That cost a lot of guessing in OpenMS/OpenMS#10141, where the answer was a
   ``libcurl.dll`` that no search directory pointed at. On failure this script
   prints the DLL search path actually in effect, what each directory holds, and
   (when ``dumpbin`` is available, as it is inside the MSVC environment the wheel
   is built in) the dependency names that resolve nowhere.

Exits non-zero when anything is missing, so it can be the first COMMAND of the
``pyopenms_stubs`` target.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

# The extensions built by src/pyOpenMS/CMakeLists.txt that are not domain
# modules. Keep in sync with the nanobind_add_module() calls there.
EXTRA_MODULES = ("_pyopenms", "_arrow_zerocopy")


def _search_dirs(package_dir: Path) -> list[Path]:
    """Directories the loader will actually look in, in order.

    Mirrors what ``pyopenms/__init__.py`` registers with ``os.add_dll_directory``
    plus the module's own directory, which Python always searches.
    """
    dirs = [package_dir]
    for entry in os.environ.get("PYOPENMS_DLL_PATH", "").split(os.pathsep):
        if entry:
            dirs.append(Path(entry))
    if sys.platform == "win32":
        # The application directory is part of the default search order, and it
        # is where python3.dll lives -- a load-time import of every stable-ABI
        # module. Without it the report would list python3.dll as unresolvable
        # fifteen times over and bury the real culprit.
        exe_dir = Path(sys.executable).resolve().parent
        dirs += [exe_dir, exe_dir / "DLLs"]
        system32 = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32"
        dirs.append(system32)
    return dirs


def _direct_dependencies(binary: Path) -> list[str] | None:
    """Names imported by ``binary``, via dumpbin. None when unavailable."""
    dumpbin = shutil.which("dumpbin")
    if dumpbin is None:
        return None
    try:
        out = subprocess.run(
            [dumpbin, "/dependents", str(binary)],
            capture_output=True,
            text=True,
            timeout=120,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    names = []
    # dumpbin prints one indented "name.dll" per line under a header; anything
    # else on the line disqualifies it.
    for line in out.splitlines():
        stripped = line.strip()
        if stripped.lower().endswith(".dll") and " " not in stripped:
            names.append(stripped)
    return names


def _report_windows_dependencies(package_dir: Path, dirs: list[Path]) -> None:
    """Name the unresolvable dependencies of the built binaries, if we can."""
    binaries = sorted(package_dir.glob("*.pyd")) + sorted(package_dir.glob("*.dll"))
    if not binaries:
        return

    resolvable = set()
    for directory in dirs:
        if directory.is_dir():
            for entry in directory.iterdir():
                if entry.suffix.lower() in (".dll", ".pyd"):
                    resolvable.add(entry.name.lower())

    checked: set[Path] = set()
    pending = list(binaries)
    unresolved: dict[str, set[str]] = {}
    dumpbin_missing = False

    while pending:
        binary = pending.pop()
        if binary in checked:
            continue
        checked.add(binary)
        names = _direct_dependencies(binary)
        if names is None:
            dumpbin_missing = True
            break
        for name in names:
            lowered = name.lower()
            # api-ms-* and ext-ms-* are API sets resolved by the loader itself and
            # never exist as files; they are not evidence of a missing library.
            if lowered.startswith(("api-ms-", "ext-ms-")):
                continue
            if lowered in resolvable:
                for directory in dirs:
                    candidate = directory / name
                    if candidate.is_file() and candidate not in checked:
                        pending.append(candidate)
                        break
            else:
                unresolved.setdefault(name, set()).add(binary.name)

    if dumpbin_missing:
        print(
            "  (dumpbin is not on PATH, so the dependency closure could not be "
            "walked; run this from a Visual Studio developer shell to get the "
            "missing DLL named)",
            file=sys.stderr,
        )
        return

    if unresolved:
        print("\nDependencies that resolve in NO search directory:", file=sys.stderr)
        for name in sorted(unresolved):
            needed_by = ", ".join(sorted(unresolved[name]))
            print(f"  {name}   (needed by {needed_by})", file=sys.stderr)
        print(
            "\nCMake can only list directories it was told about. Put the prefix "
            "holding the missing DLL on CMAKE_PREFIX_PATH (its bin/ and lib/ are "
            "then searched), or name the directory in PYOPENMS_DLL_PATH "
            "(os.pathsep-separated), or configure with -DPYOPENMS_GENERATE_STUBS=OFF.",
            file=sys.stderr,
        )
    else:
        print(
            "\nEvery direct dependency resolves in the search directories above; "
            "the failure is elsewhere (an ABI mismatch or a failing module init).",
            file=sys.stderr,
        )


def _report_environment(package_dir: Path, dirs: list[Path]) -> None:
    """Print what the loader had to work with, to stderr.

    Interpreter, search directories with the libraries each one holds, and on
    Windows the dependency names that resolve nowhere -- the information the
    loader's own error message leaves out.
    """
    print("\n--- pyOpenMS import diagnostics ---", file=sys.stderr)
    print(f"interpreter : {sys.executable}", file=sys.stderr)
    print(f"version     : {sys.version.splitlines()[0]}", file=sys.stderr)
    print(f"platform    : {sys.platform}", file=sys.stderr)
    print(f"package dir : {package_dir}", file=sys.stderr)

    raw = os.environ.get("PYOPENMS_DLL_PATH")
    # Printed from Python, not through the build tool's shell: the ';' separators
    # and any '(' in a path are literal here. Echoing this value from CMake is what
    # broke the stub target with "/bin/sh: Syntax error" once already.
    print(f"PYOPENMS_DLL_PATH: {raw if raw is not None else '(not set)'}", file=sys.stderr)

    print("\nDLL search directories (in order):", file=sys.stderr)
    for directory in dirs:
        if not directory.is_dir():
            print(f"  [missing] {directory}", file=sys.stderr)
            continue
        libs = sorted(
            entry.name
            for entry in directory.iterdir()
            if entry.suffix.lower() in (".dll", ".so", ".dylib", ".pyd")
        )
        shown = ", ".join(libs[:20]) + (", ..." if len(libs) > 20 else "")
        print(f"  [ok]      {directory}", file=sys.stderr)
        print(f"            {len(libs)} libraries: {shown or '(none)'}", file=sys.stderr)

    if sys.platform == "win32":
        _report_windows_dependencies(package_dir, dirs)


def main() -> int:
    """Import every expected module; return 0 if all load, 1 with diagnostics otherwise."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package-dir",
        required=True,
        type=Path,
        help="the built pyopenms package directory (containing the extensions)",
    )
    parser.add_argument(
        "--domains",
        required=True,
        nargs="+",
        help="domain names, i.e. the _pyopenms_<domain> modules that must import",
    )
    opt = parser.parse_args()

    package_dir = opt.package_dir.resolve()
    expected = [f"_pyopenms_{domain}" for domain in opt.domains] + list(EXTRA_MODULES)

    # Import the package the way stubgen will: by name, from its parent directory.
    sys.path.insert(0, str(package_dir.parent))

    failures: list[tuple[str, BaseException]] = []
    try:
        importlib.import_module(package_dir.name)
    except BaseException as exc:  # noqa: BLE001 - report anything, including SystemExit
        # __init__.py raises only when nothing imported; that is already fatal.
        print(
            f"error: `import {package_dir.name}` failed: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        traceback.print_exc()
        failures.append((package_dir.name, exc))

    if not failures:
        for name in expected:
            try:
                importlib.import_module(f"{package_dir.name}.{name}")
            except BaseException as exc:  # noqa: BLE001
                failures.append((name, exc))

    if not failures:
        print(f"pyOpenMS: all {len(expected)} extension modules import correctly")
        return 0

    print(
        f"\nerror: {len(failures)} of {len(expected)} expected pyOpenMS extension "
        f"modules failed to import; refusing to generate stubs from an "
        f"incomplete package.",
        file=sys.stderr,
    )
    for name, exc in failures:
        print(f"  {name}: {type(exc).__name__}: {exc}", file=sys.stderr)

    _report_environment(package_dir, _search_dirs(package_dir))
    return 1


if __name__ == "__main__":
    sys.exit(main())
