#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Copy the license files of the shared libraries a pyOpenMS wheel bundles into the
OpenMS installation the wheel is built from, before the wheel is built.

The wheel repair tools copy every shared library that libOpenMS needs and the target
system does not provide into the wheel: auditwheel on Linux, delocate on macOS. Their
licenses have to travel with them. This script walks the dependencies of the shared
libraries in <prefix>/lib (and lib64) as the loader resolves them, finds the package
that installed each one, and copies that package's license files to
<prefix>/share/OpenMS/LICENSES, which pyOpenMS puts into the wheel:

- Linux: the RPM package that owns the library (the manylinux images are RPM based),
  and its license files (rpm --licensefiles, else its COPYING/LICENSE/NOTICE files),
  to LICENSES/system/<package>/. A package without any, such as libgomp, takes those of
  an installed package built from the same source package (libgcc, from gcc).
- macOS: the Homebrew keg that holds the library (Cellar/<formula>/<version>), and its
  COPYING/LICENSE/NOTICE files, to LICENSES/homebrew/<formula>/.

--license-file <name>=<file> gives the license text for a package or formula that ships
none: Homebrew copies only files named exactly LICENSE, COPYING and the like into a keg,
so boost's LICENSE_1_0.txt does not reach it.

Libraries of the installation itself are skipped, and so are the C and C++ runtime
(glibc, libgcc_s, libstdc++; on macOS everything under /usr/lib and /System), which the
repair tools leave out as well. Libraries under a --skip-prefix are skipped, too: the
OpenMS build installs the licenses of the vcpkg ports itself
(cmake/third_party_licenses.cmake).

Exits with 1, naming the library, if a bundled library has no package or no license file.

Usage: collect_wheel_licenses.py --prefix <install prefix> [--skip-prefix <dir> ...]
                                 [--homebrew-formula <header-only formula> ...]
                                 [--license-file <package or formula>=<file> ...]
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

# The C and C++ runtime of the manylinux platform. auditwheel never bundles these.
LINUX_SYSTEM_LIBRARIES = re.compile(
    r"^(ld-linux.*|linux-vdso|linux-gate|libc|libm|libdl|librt|libpthread|libutil|libresolv"
    r"|libnsl|libgcc_s|libstdc\+\+|libmvec)\.so(\.|$)")
LICENSE_NAME = re.compile(r"^(licen[cs]e|copying|copyright|notice)([-_.0-9].*)?$", re.IGNORECASE)


def run(*command):
    return subprocess.run(command, capture_output=True, text=True, check=False)


def own_libraries(prefix):
    """The shared libraries the installation provides, i.e. the roots of the walk."""
    roots = []
    for libdir in ("lib", "lib64"):
        folder = os.path.join(prefix, libdir)
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            path = os.path.join(folder, name)
            if os.path.isfile(path) and not os.path.islink(path) and \
               (re.search(r"\.so(\.\d+)*$", name) or name.endswith(".dylib")):
                roots.append(path)
    return roots


def under(path, prefixes):
    real = os.path.realpath(path)
    return any(real == p or real.startswith(p + os.sep) for p in prefixes)


# ------------------------------------------------------------------------------ Linux
def linux_dependencies(roots):
    """{soname: resolved path} of everything the roots load, as ldd reports it."""
    found, missing = {}, []
    for root in roots:
        for line in run("ldd", root).stdout.splitlines():
            match = re.match(r"\s*(\S+)\s+=>\s+(\S+)", line)
            if not match:
                continue
            soname, path = match.groups()
            if path == "not":
                missing.append(soname)
            else:
                found.setdefault(soname, path)
    return found, missing


def rpm_own_license_files(package):
    files = [f for f in run("rpm", "-q", "--licensefiles", package).stdout.splitlines()
             if f.startswith("/") and os.path.isfile(f)]
    if not files:  # packages that predate %license list them as documentation
        files = [f for f in run("rpm", "-q", "--docfiles", package).stdout.splitlines()
                 if f.startswith("/") and os.path.isfile(f) and LICENSE_NAME.match(os.path.basename(f))]
    return files


def rpm_license_files(package):
    """(license files, the package they belong to) for an RPM package.

    Subpackages of one source package can leave the license files to one of them: gcc's
    libgomp has none, libgcc has gcc's. Such a package takes those of an installed
    package built from the same source package.
    """
    files = rpm_own_license_files(package)
    if files:
        return files, package
    source = run("rpm", "-q", "--queryformat", "%{SOURCERPM}", package).stdout.strip()
    if source:
        installed = run("rpm", "-qa", "--queryformat", "%{NAME}\t%{SOURCERPM}\n").stdout.splitlines()
        for sibling in sorted(line.split("\t")[0] for line in installed if line.endswith("\t" + source)):
            files = rpm_own_license_files(sibling) if sibling != package else []
            if files:
                return files, sibling
    return [], None


def collect_linux(roots, skip, dest, supplied):
    dependencies, missing = linux_dependencies(roots)
    problems = [f"{soname}: not found by the loader" for soname in missing]
    packages = {}
    for soname, path in sorted(dependencies.items()):
        if LINUX_SYSTEM_LIBRARIES.match(soname) or under(path, skip):
            continue
        owner = run("rpm", "-qf", "--queryformat", "%{NAME}\t%{VERSION}-%{RELEASE}\n",
                    os.path.realpath(path))
        if owner.returncode != 0 or "\t" not in owner.stdout:
            problems.append(f"{soname} ({path}): no RPM package owns it")
            continue
        name, version = owner.stdout.splitlines()[0].split("\t")
        packages.setdefault((name, version), []).append(soname)
    index = []
    for (name, version), sonames in sorted(packages.items()):
        files, owner = rpm_license_files(name)
        note = f" (license files of {owner}, built from the same source package)" if owner not in (name, None) else ""
        if not files and name in supplied:
            files, note = [supplied[name]], " (license text supplied by the build; the package has none)"
        if not files:
            problems.append(f"{', '.join(sonames)}: package {name} has no license file")
            continue
        copy_files(files, os.path.join(dest, "system", name))
        index.append(f"{', '.join(sorted(sonames))}\t{name} {version}{note}")
    write_index(os.path.join(dest, "system"), index,
                "Libraries of the build system that the wheel bundles, and the RPM package each comes from.")
    return problems


# ------------------------------------------------------------------------------ macOS
def macho_load_commands(path):
    """(dependencies, rpaths) of a Mach-O file, from otool -l."""
    output = run("otool", "-l", path).stdout.splitlines()
    deps, rpaths = [], []
    for i, line in enumerate(output):
        command = line.strip()
        if command in ("cmd LC_LOAD_DYLIB", "cmd LC_LOAD_WEAK_DYLIB", "cmd LC_REEXPORT_DYLIB",
                       "cmd LC_RPATH"):
            for follow in output[i + 1:i + 4]:
                match = re.match(r"\s*(name|path) (\S+) \(offset", follow)
                if match:
                    (rpaths if command == "cmd LC_RPATH" else deps).append(match.group(2))
                    break
    loader = os.path.dirname(path)
    rpaths = [r.replace("@loader_path", loader).replace("@executable_path", loader) for r in rpaths]
    return deps, rpaths


def resolve_macho(reference, loader_path, rpaths):
    if reference.startswith("@rpath/"):
        # The loader's own folder last: the libraries of the installation sit side by side.
        for rpath in rpaths + [os.path.dirname(loader_path)]:
            candidate = os.path.join(rpath, reference[len("@rpath/"):])
            if os.path.exists(candidate):
                return candidate
        return None
    if reference.startswith(("@loader_path/", "@executable_path/")):
        candidate = os.path.join(os.path.dirname(loader_path), reference.split("/", 1)[1])
        return candidate if os.path.exists(candidate) else None
    return reference if os.path.exists(reference) else None


def macos_dependencies(roots):
    """Resolved paths of the non-system dylibs the roots load, directly or indirectly.

    dyld resolves an @rpath reference against the run paths of the whole chain of images
    that loads the library: its own LC_RPATHs, then those of the images that led to it.
    A library reached again with run paths it has not been tried with is walked again.
    """
    load_commands = {}
    seen, stack = set(), [(root, ()) for root in roots]
    tried = {}                 # path -> the run paths it was walked with
    failed, resolved = set(), set()
    while stack:
        current, inherited = stack.pop()
        if current not in load_commands:
            load_commands[current] = macho_load_commands(current)
        deps, own = load_commands[current]
        rpaths = tuple(dict.fromkeys(list(own) + list(inherited)))
        if current in tried and set(rpaths) <= tried[current]:
            continue
        tried.setdefault(current, set()).update(rpaths)
        for reference in deps:
            if reference.startswith(("/usr/lib/", "/System/")):
                continue
            path = resolve_macho(reference, current, list(rpaths))
            if path is None:
                failed.add((reference, current))
                continue
            resolved.add((reference, current))
            real = os.path.realpath(path)
            seen.add(real)
            stack.append((real, rpaths))
    unresolved = [f"{reference} (needed by {loader})" for reference, loader in sorted(failed - resolved)]
    return sorted(seen), unresolved


def keg_of(path):
    match = re.match(r"^(.*/Cellar/([^/]+)/[^/]+)/", path)
    return (match.group(1), match.group(2)) if match else (None, None)


def keg_license_files(keg, formula):
    for folder in (keg, os.path.join(keg, "share", "doc", formula), os.path.join(keg, "share", "licenses", formula)):
        if os.path.isdir(folder):
            files = [os.path.join(folder, f) for f in sorted(os.listdir(folder))
                     if LICENSE_NAME.match(f) and os.path.isfile(os.path.join(folder, f))]
            if files:
                return files
    return []


def collect_macos(roots, skip, dest, compiled_in, supplied):
    dependencies, unresolved = macos_dependencies(roots)
    problems = [f"{reference}: cannot be resolved" for reference in unresolved]
    kegs = {}
    # Header-only formulae leave no dylib to follow; their code is compiled into libOpenMS.
    for formula in compiled_in:
        opt = run("brew", "--prefix", formula).stdout.strip()
        keg, name = keg_of(os.path.realpath(opt) + os.sep) if opt else (None, None)
        if keg is None:
            problems.append(f"{formula}: no Homebrew keg found")
        else:
            kegs.setdefault((keg, name), []).append("(compiled in)")
    for path in dependencies:
        if under(path, skip) or path in roots:
            continue
        keg, formula = keg_of(path)
        if keg is None:
            problems.append(f"{path}: not in a Homebrew keg")
            continue
        kegs.setdefault((keg, formula), []).append(os.path.basename(path))
    index = []
    for (keg, formula), libraries in sorted(kegs.items()):
        files, note = keg_license_files(keg, formula), ""
        if not files and formula in supplied:
            files, note = [supplied[formula]], " (license text supplied by the build; the keg has none)"
        if not files:
            problems.append(f"{', '.join(libraries)}: the keg {keg} has no license file")
            continue
        copy_files(files, os.path.join(dest, "homebrew", formula))
        index.append(f"{', '.join(sorted(libraries))}\t{formula} {os.path.basename(keg)}{note}")
    write_index(os.path.join(dest, "homebrew"), index,
                "Homebrew libraries that the wheel bundles, and the formula each comes from.")
    return problems


# ------------------------------------------------------------------------------ shared
def copy_files(files, folder):
    os.makedirs(folder, exist_ok=True)
    for source in files:
        target = os.path.join(folder, os.path.basename(source))
        if os.path.exists(target):  # two license files of one package with the same name
            target = os.path.join(folder, os.path.basename(os.path.dirname(source)) + "-" + os.path.basename(source))
        shutil.copyfile(source, target)


def write_index(folder, lines, header):
    if not lines:
        return
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "INDEX.txt"), "w", encoding="utf-8") as index:
        index.write(header + "\n\n" + "\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prefix", required=True, help="the OpenMS installation the wheel is built from")
    parser.add_argument("--skip-prefix", action="append", default=[],
                        help="libraries below this directory are skipped (their licenses come from elsewhere)")
    parser.add_argument("--homebrew-formula", action="append", default=[],
                        help="macOS: a header-only formula whose code libOpenMS compiles in, e.g. eigen")
    parser.add_argument("--license-file", action="append", default=[], metavar="NAME=FILE",
                        help="the license text for a package or formula that ships none")
    args = parser.parse_args()
    supplied = {}
    for entry in args.license_file:
        name, _, path = entry.partition("=")
        if not name or not os.path.isfile(path):
            parser.error(f"--license-file {entry}: expected <package or formula>=<existing file>")
        supplied[name] = path

    prefix = os.path.realpath(args.prefix)
    skip = [prefix] + [os.path.realpath(p) for p in args.skip_prefix]
    dest = os.path.join(prefix, "share", "OpenMS", "LICENSES")
    roots = own_libraries(prefix)
    if not roots:
        print(f"collect_wheel_licenses: no shared libraries in {prefix}/lib", file=sys.stderr)
        return 1
    if sys.platform == "darwin":
        problems = collect_macos(roots, skip, dest, args.homebrew_formula, supplied)
    elif sys.platform.startswith("linux"):
        problems = collect_linux(roots, skip, dest, supplied)
    else:
        print(f"collect_wheel_licenses: nothing to do on {sys.platform}")
        return 0
    for folder in ("system", "homebrew"):
        index = os.path.join(dest, folder, "INDEX.txt")
        if os.path.exists(index):
            print(open(index, encoding="utf-8").read())
    if problems:
        print("collect_wheel_licenses: no license found for:", file=sys.stderr)
        for problem in problems:
            print("  " + problem, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
