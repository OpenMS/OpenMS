"""Verify that every compiled pyOpenMS extension module imports.

``nanobind.stubgen`` generates stubs by *importing* the freshly built package.
``pyopenms/__init__.py`` deliberately swallows per-module import errors and only
raises when no module at all could be loaded, so stubgen can silently succeed on
a partial import and emit stubs that are missing whole domains.  This script is
run before stubgen and fails the build unless all expected modules are there.

On Windows it also explains *why* an import failed.  The loader only ever says

    ImportError: DLL load failed while importing _pyopenms_analysis:
                 The specified module could not be found.

which never names the DLL that is actually missing.  We therefore walk the
import table of the extension module and resolve each dependency against the
same directories the loader uses, and report the ones that cannot be found.

Usage:
    check_module_imports.py <package-dir> <module> [<module> ...]
"""

from __future__ import annotations

import importlib
import os
import struct
import sys
from importlib.machinery import EXTENSION_SUFFIXES

# DLLs that always come from the OS; never worth reporting as missing.
_SYSTEM_PREFIXES = ("api-ms-", "ext-ms-")
_SYSTEM_NAMES = frozenset(
    (
        "advapi32.dll", "bcrypt.dll", "cfgmgr32.dll", "combase.dll", "crypt32.dll",
        "gdi32.dll", "iphlpapi.dll", "kernel32.dll", "kernelbase.dll", "msvcrt.dll",
        "normaliz.dll", "ntdll.dll", "ole32.dll", "oleaut32.dll", "rpcrt4.dll",
        "secur32.dll", "shell32.dll", "shlwapi.dll", "user32.dll", "userenv.dll",
        "version.dll", "ws2_32.dll", "wldap32.dll",
    )
)


def _is_system_dll(name: str) -> bool:
    lowered = name.lower()
    return lowered in _SYSTEM_NAMES or lowered.startswith(_SYSTEM_PREFIXES)


# --------------------------------------------------------------------------
# Minimal PE import-table reader (stdlib only; pefile is not a build dependency)
# --------------------------------------------------------------------------
def _pe_imported_dlls(path: str) -> list[str]:
    """Return the DLL names in the import and delay-import tables of a PE file."""
    with open(path, "rb") as handle:
        data = handle.read()

    if data[:2] != b"MZ":
        return []
    pe_off = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe_off : pe_off + 4] != b"PE\0\0":
        return []

    coff = pe_off + 4
    n_sections, = struct.unpack_from("<H", data, coff + 2)
    opt_size, = struct.unpack_from("<H", data, coff + 16)
    opt = coff + 20
    magic, = struct.unpack_from("<H", data, opt)
    if magic == 0x20B:  # PE32+
        image_base, = struct.unpack_from("<Q", data, opt + 24)
        dir_off = opt + 112
    elif magic == 0x10B:  # PE32
        image_base, = struct.unpack_from("<I", data, opt + 28)
        dir_off = opt + 96
    else:
        return []

    sections = []
    sec_off = opt + opt_size
    for i in range(n_sections):
        base = sec_off + i * 40
        virt_size, virt_addr, raw_size, raw_ptr = struct.unpack_from("<IIII", data, base + 8)
        sections.append((virt_addr, max(virt_size, raw_size), raw_ptr))

    def to_offset(rva: int) -> int | None:
        for virt_addr, size, raw_ptr in sections:
            if virt_addr <= rva < virt_addr + size:
                return raw_ptr + (rva - virt_addr)
        return None

    def read_name(rva: int) -> str | None:
        # Delay-import descriptors predating v2 store virtual addresses, not RVAs.
        if rva > image_base:
            rva -= image_base
        off = to_offset(rva)
        if off is None or off >= len(data):
            return None
        end = data.find(b"\0", off)
        return data[off:end].decode("ascii", "replace") if end != -1 else None

    names: list[str] = []
    # Data directory 1 = import table (20-byte entries, name RVA at +12),
    # data directory 13 = delay-import table (32-byte entries, name RVA at +4).
    for index, entry_size, name_field in ((1, 20, 12), (13, 32, 4)):
        table_rva, table_size = struct.unpack_from("<II", data, dir_off + index * 8)
        if not table_rva or not table_size:
            continue
        cursor = to_offset(table_rva)
        if cursor is None:
            continue
        while cursor + entry_size <= len(data):
            entry = data[cursor : cursor + entry_size]
            if entry == b"\0" * entry_size:
                break
            name = read_name(struct.unpack_from("<I", entry, name_field)[0])
            if name:
                names.append(name)
            cursor += entry_size
    return names


def _search_directories(module_path: str) -> list[str]:
    """The directories the Windows loader consults for an extension's dependencies.

    CPython loads extension modules with LOAD_LIBRARY_SEARCH_DEFAULT_DIRS |
    LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR, i.e. the module's own directory, the
    directories registered with os.add_dll_directory() (which pyopenms fills
    from PYOPENMS_DLL_PATH), the application directory and the system
    directories -- but *not* PATH.
    """
    dirs = [os.path.dirname(os.path.abspath(module_path))]
    dirs += [d for d in os.environ.get("PYOPENMS_DLL_PATH", "").split(os.pathsep) if d]
    dirs.append(os.path.dirname(sys.executable))
    windir = os.environ.get("SystemRoot", r"C:\Windows")
    dirs += [os.path.join(windir, "System32"), os.path.join(windir, "SysWOW64"), windir]

    seen, unique = set(), []
    for directory in dirs:
        key = os.path.normcase(os.path.abspath(directory))
        if key not in seen and os.path.isdir(directory):
            seen.add(key)
            unique.append(directory)
    return unique


def _unresolved_dlls(module_path: str) -> tuple[list[str], list[tuple[str, str]]]:
    """Return the searched directories and the dependencies that are not in them."""
    directories = _search_directories(module_path)

    def locate(name: str) -> str | None:
        for directory in directories:
            candidate = os.path.join(directory, name)
            if os.path.isfile(candidate):
                return candidate
        return None

    # Breadth-first over the dependency graph: a module fails to load when *any*
    # DLL in its transitive closure is missing, not just a direct dependency.
    pending = [(os.path.basename(module_path), module_path)]
    visited = {os.path.normcase(os.path.basename(module_path))}
    unresolved: list[tuple[str, str]] = []
    while pending:
        owner, path = pending.pop(0)
        for name in _pe_imported_dlls(path):
            key = os.path.normcase(name)
            if key in visited or _is_system_dll(name):
                continue
            visited.add(key)
            found = locate(name)
            if found is None:
                unresolved.append((name, owner))
            else:
                pending.append((name, found))
    return directories, unresolved


def _module_path(package_dir: str, name: str) -> str | None:
    for suffix in EXTENSION_SUFFIXES:
        candidate = os.path.join(package_dir, name + suffix)
        if os.path.isfile(candidate):
            return candidate
    return None


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    package_dir = os.path.abspath(argv[1])
    expected = argv[2:]

    # Import the package the way stubgen will: by name, from its parent directory.
    sys.path.insert(0, os.path.dirname(package_dir))
    package_name = os.path.basename(package_dir)

    failure = None
    try:
        __import__(package_name)
    except Exception as error:  # noqa: BLE001 - reported below, not handled
        failure = f"{type(error).__name__}: {error}"

    # __init__.py only pulls in the _pyopenms_<domain> modules; the main module
    # and _arrow_zerocopy are imported lazily by callers, so ask for each one
    # explicitly rather than trusting what ended up in sys.modules.
    missing = {}
    for name in expected:
        full = f"{package_name}.{name}"
        if full in sys.modules:
            continue
        try:
            importlib.import_module(full)
        except Exception as error:  # noqa: BLE001 - reported below, not handled
            missing[name] = f"{type(error).__name__}: {error}"

    if not missing and failure is None:
        print(f"pyOpenMS: all {len(expected)} extension modules imported")
        return 0

    report = [
        "=" * 78,
        f"{package_name}: the freshly built package cannot be imported in full, so",
        "stub generation would produce incomplete or no type information.",
        "=" * 78,
    ]
    if failure is not None:
        report.append(f"  importing {package_name} raised {failure}")
    for name, reason in missing.items():
        report.append(f"  {name}: {reason}")
        path = _module_path(package_dir, name)
        if path is None:
            report.append(f"      no extension module file for {name} in {package_dir}")
        elif sys.platform == "win32":
            # "DLL load failed ... The specified module could not be found" never
            # names the DLL that is missing. Work it out from the import table.
            directories, unresolved = _unresolved_dlls(path)
            if unresolved:
                report.append("      unresolved dependencies:")
                report += [f"        {dll}  (needed by {owner})" for dll, owner in unresolved]
            else:
                report.append("      every dependency in the import table was located, "
                              "so this is not a missing DLL")
            report.append("      searched, in loader order:")
            report += [f"        {directory}" for directory in directories]
    report.append(f"  python: {sys.executable}")
    report.append(f"  PYOPENMS_DLL_PATH={os.environ.get('PYOPENMS_DLL_PATH', '<unset>')}")
    print("\n".join(report), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
