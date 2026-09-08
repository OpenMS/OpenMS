#!/usr/bin/env python3
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

"""Check the proposed core's textual include boundary without configuring OpenMS.

The core is derived, not hand-classified: starting from the seed headers, follow
header includes, each header's mirrored implementation file, and that file's
includes, stopping at the edges recorded in baseline.json. The result is
snapshotted in core_files.txt so that CI can name every crossing edge; regenerate
the snapshot with --update-core after intentionally changing the boundary.
"""

import argparse
import json
import posixpath
import re
import shlex
import sys
from pathlib import Path


LIBRARIES = ("src/openms", "src/openswathalgo")
HEADER_SUFFIXES = {".h", ".hpp", ".hxx", ".ipp", ".tpp"}
SOURCE_SUFFIXES = {".cpp", ".cc", ".cxx", ".c"}
DEFAULT_DIRECTORY = Path(__file__).resolve().parent / "core_dependencies"
LISTED_FILES = 10  # cap per-file diagnostics; the crossing edges are the actionable part
# Installed value-type headers the core grows from.
SEEDS = (
    "src/openms/include/OpenMS/KERNEL/Peak1D.h",
    "src/openms/include/OpenMS/KERNEL/ChromatogramPeak.h",
    "src/openms/include/OpenMS/KERNEL/MSSpectrum.h",
    "src/openms/include/OpenMS/KERNEL/MSChromatogram.h",
    "src/openms/include/OpenMS/KERNEL/MSExperiment.h",
)
# Headers produced by CMake; they need not exist in a source checkout.
GENERATED_HEADERS = {
    "OpenMS/config.h": "src/openms/configh.cmake",
    "OpenMS/OpenMSConfig.h": "cmake/add_library_macros.cmake",
    "OpenMS/openms_package_version.h": "src/openms/configh.cmake",
    "OpenMS/openms_data_path.h": "src/openms/configh.cmake",
    "OpenMS/build_config.h": "src/openms/configh.cmake",
    "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h": "cmake/add_library_macros.cmake",
}
# Preserve ordinary string/character literals, but mask comments and raw strings.
# Masking keeps newlines, so include locations remain useful in diagnostics.
TOKENS = re.compile(
    r'//[^\n]*|/\*.*?\*/|(?:u8|u|U|L)?R"([^ ()\\\t\r\n]{0,16})\(.*?\)\1"'
    r'|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
    re.DOTALL,
)


def mask_comments(text, raw_strings=False):
    def mask(match):
        token = match.group()
        if token.startswith(("//", "/*")) or (raw_strings and match.group(1) is not None):
            return re.sub(r"[^\n]", " ", token)
        return token

    return TOKENS.sub(mask, text)


def logical_lines(text):
    """Apply backslash-newline splicing before tokenization, retaining line numbers."""
    lines, numbers, pending = [], [], ""
    first = 1
    for number, line in enumerate(text.splitlines(keepends=True), 1):
        if not pending:
            first = number
        if line.endswith("\\\n"):
            pending += line[:-2]
        else:
            lines.append(pending + line)
            numbers.append(first)
            pending = ""
    if pending:
        lines.append(pending)
        numbers.append(first)
    return "".join(lines), numbers


def includes(text):
    spliced, numbers = logical_lines(text)
    clean = mask_comments(spliced, raw_strings=True)
    for number, line in zip(numbers, clean.splitlines()):
        directive = re.match(r"^\s*#\s*include\s*(.*)$", line)
        if not directive:
            continue
        value = directive[1].strip()
        literal = re.fullmatch(r'<([^>]+)>|"([^"]+)"', value)
        if literal:
            yield number, literal[1] or literal[2], bool(literal[2])
        else:
            yield number, value, None


def cmake_values(text, variable):
    """Union literal set/list(APPEND) values across all conditional branches.

    This deliberately supports the repository's install-list idiom, not CMake
    execution. Reject computed values/other mutations instead of guessing.
    """
    text = re.sub(r"#[^\n]*", "", text)
    values = set()
    for command, body in re.findall(r"\b(set|list|unset)\s*\(([^)]*)\)", text, re.I):
        tokens = shlex.split(body, posix=True)
        if not tokens:
            continue
        if command.lower() == "list":
            if len(tokens) < 2 or tokens[1] != variable:
                continue
            if tokens[0].upper() != "APPEND":
                raise ValueError(f"Unsupported CMake mutation of {variable}: {body.strip()}")
            items = tokens[2:]
        else:
            if tokens[0] != variable:
                continue
            if command.lower() != "set":
                raise ValueError(f"Unsupported CMake mutation of {variable}")
            items = tokens[1:]
        for item in items:
            for value in item.split(";"):
                if not re.fullmatch(r"[A-Za-z0-9_./+-]+", value):
                    raise ValueError(f"Nonliteral CMake value for {variable}: {value}")
                values.add(value)
    return values


def installed_headers(root):
    """Read the lists actually passed as HEADER_FILES for the two libraries."""
    installed = set()
    entry = root / "src/openms/includes.cmake"
    text = re.sub(r"#[^\n]*", "", entry.read_text())
    lists = re.findall(r"\binclude\s*\(\s*(include/OpenMS/[^\s)]+/sources\.cmake)\s*\)", text)
    if not lists:
        raise ValueError(f"No header install lists found in {entry}")
    for name in lists:
        path = root / "src/openms" / name
        content = path.read_text()
        directories = cmake_values(content, "directory")
        if len(directories) != 1:
            raise ValueError(f"Expected one literal header directory in {path}")
        directory = root / "src/openms" / directories.pop()
        for value in cmake_values(content, "sources_list_h"):
            installed.add((directory / value).relative_to(root).as_posix())

    path = root / "src/openswathalgo/source/OPENSWATHALGO/OpenSwathAlgoFiles.cmake"
    text = path.read_text()
    directories = cmake_values(text, "header_directory")
    if len(directories) != 1:
        raise ValueError("Expected one literal OpenSwathAlgo header_directory")
    directory = root / "src/openswathalgo" / directories.pop()
    for variable in ("header_algo_list", "header_dataaccess_list"):
        values = cmake_values(text, variable)
        if not values:
            raise ValueError(f"Empty OpenSwathAlgo install list: {variable}")
        installed.update((directory / value).relative_to(root).as_posix() for value in values)
    for library, variable in zip(LIBRARIES, ("OpenMS_sources_h", "OpenSwathAlgoHeaders")):
        text = (root / library / "CMakeLists.txt").read_text()
        if not re.search(r"HEADER_FILES\s+\$\{" + variable + r"\}", text):
            raise ValueError(f"Header install convention changed in {library}/CMakeLists.txt")
    return installed


def implementation_of(header):
    """The translation unit mirroring a header, e.g. include/OpenMS/X/A.h -> source/X/A.cpp."""
    match = re.fullmatch(r"(src/[^/]+)/include/OpenMS/(.+)\.[^.]+", header)
    return f"{match[1]}/source/{match[2]}.cpp" if match else None


def reachable(graph, seeds):
    seen, pending = set(), list(seeds)
    while pending:
        node = pending.pop()
        if node not in seen:
            seen.add(node)
            pending.extend(graph.get(node, set()) - seen)
    return sorted(seen)


def derive_core(graph, files, stop):
    """Closure over includes and mirrored implementations, not crossing baselined edges."""
    core, pending = set(), list(SEEDS)
    while pending:
        node = pending.pop()
        if node in core:
            continue
        core.add(node)
        implementation = implementation_of(node)
        if implementation in files:
            pending.append(implementation)
        pending.extend(target for target in graph.get(node, set()) if (node, target) not in stop)
    return sorted(core)


def cycles(graph):
    """Iterative Kosaraju SCCs; do not merge a header and its translation units."""
    seen, order, reverse = set(), [], {node: set() for node in graph}
    for source, targets in graph.items():
        for target in targets:
            reverse[target].add(source)
    for start in sorted(graph):
        stack = [(start, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                order.append(node)
            elif node not in seen:
                seen.add(node)
                stack.append((node, True))
                stack.extend((target, False) for target in sorted(graph[node], reverse=True))
    seen, result = set(), []
    for start in reversed(order):
        if start in seen:
            continue
        component, stack = [], [start]
        while stack:
            node = stack.pop()
            if node not in seen:
                seen.add(node)
                component.append(node)
                stack.extend(reverse[node] - seen)
        if len(component) > 1 or start in graph[start]:
            result.append(sorted(component))
    return sorted(result)


def edge_key(edge):
    return edge["rule"], edge["source"], edge["target"]


def baseline_edges(baseline):
    """Validate baseline.json and return its exception entries."""
    if baseline.get("schema_version") != 1:
        raise ValueError("Unsupported baseline schema_version")
    entries = baseline["edges"]
    keys = [edge_key(edge) for edge in entries]
    if len(keys) != len(set(keys)) or any(not edge.get("reason", "").strip() for edge in entries):
        raise ValueError("Baseline edges must be unique and each needs a reason")
    return entries


def check_baseline(violations, entries):
    keys = {edge_key(edge) for edge in entries}
    current = {edge_key(edge): edge for edge in violations}
    return (
        [current[key] for key in sorted(current.keys() - keys)],
        [edge for edge in entries if edge_key(edge) not in current],
    )


def read_core(path):
    """Read the snapshot: one file per line, blank lines and # comments ignored."""
    lines = (line.strip() for line in path.read_text(encoding="utf-8").splitlines())
    return {line for line in lines if line and not line.startswith("#")}


def write_core(path, core):
    path.write_text(
        "# Generated by tools/ci/check_core_dependencies.py --update-core; do not edit.\n"
        "# Core closure from the seed headers, stopping at baseline.json edges.\n"
        + "".join(f"{name}\n" for name in sorted(core)),
        encoding="utf-8",
    )


def scan(root, core, stop):
    """Inventory includes and check the boundary of `core`, a set of project files."""
    files = {}
    for library in LIBRARIES:
        for part in ("include", "source"):
            directory = root / library / part
            if not directory.is_dir():
                raise ValueError(f"Missing scan directory: {directory}")
            for path in sorted(directory.rglob("*")):
                if path.suffix in HEADER_SUFFIXES | SOURCE_SUFFIXES:
                    files[path.relative_to(root).as_posix()] = path
    installed = installed_headers(root)
    if missing := installed - files.keys():
        raise ValueError("Installed headers missing from checkout: " + ", ".join(sorted(missing)))
    if missing := set(SEEDS) - files.keys():
        raise ValueError("Seed headers missing from checkout: " + ", ".join(sorted(missing)))
    if missing := core - files.keys():
        raise ValueError("Core files missing from checkout: " + ", ".join(sorted(missing)))

    graph = {name: set() for name in files}
    edges, resources, violations = [], [], {}
    for source, path in sorted(files.items()):
        content = path.read_text(encoding="utf-8")
        for line, include, quoted in includes(content):
            normalized = posixpath.normpath(include.removeprefix("include/"))
            candidates = []
            if quoted:
                candidates.append(posixpath.normpath(posixpath.join(posixpath.dirname(source), include)))
            candidates.extend(f"{library}/include/{normalized}" for library in LIBRARIES)
            target = next((candidate for candidate in candidates if candidate in files), None)
            if quoted is None:
                target, status = include, "macro"
            elif target:
                status = "internal"
                graph[source].add(target)
            elif normalized in GENERATED_HEADERS:
                target, status = normalized, "generated"
            elif normalized.startswith("OpenMS/") or quoted:
                target, status = normalized, "unresolved"
            else:
                target, status = normalized, "external"
            edge = {"source": source, "target": target, "line": line, "include": include, "status": status}
            edges.append(edge)
            rules = []
            if source in core:
                if status == "internal" and target not in core:
                    rules.append("core-to-noncore")
                elif status in {"unresolved", "macro"}:
                    rules.append("core-unresolved")
            if source in installed and status == "internal" and target not in installed:
                rules.append("installed-to-private")
            for rule in rules:
                violation = {"rule": rule, **edge}
                violations.setdefault(edge_key(violation), violation)
        # Runtime data is a separate inventory, never an include-graph edge.
        clean = mask_comments(content, raw_strings=True)
        for match in re.finditer(r'\bFile\s*::\s*find\s*\(\s*"([^"\n]+)"', clean):
            resources.append({"source": source, "line": clean.count("\n", 0, match.start()) + 1,
                              "resource": match[1], "evidence": "literal File::find argument"})

    metadata = {
        name: {"kind": "header" if path.suffix in HEADER_SUFFIXES else "implementation",
               "installed": name in installed, "core": name in core}
        for name, path in sorted(files.items())
    }
    return {
        "schema_version": 1,
        "files": metadata,
        "header_includes": [edge for edge in edges if metadata[edge["source"]]["kind"] == "header"],
        "implementation_includes": [edge for edge in edges if metadata[edge["source"]]["kind"] == "implementation"],
        "seed_header_reachability": {seed: reachable(graph, [seed]) for seed in SEEDS},
        "core_files": derive_core(graph, files.keys(), stop),
        "cycles": cycles(graph),
        "runtime_resources": resources,
        "violations": [violations[key] for key in sorted(violations)],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--baseline", type=Path, default=DEFAULT_DIRECTORY / "baseline.json")
    parser.add_argument("--core", type=Path, default=DEFAULT_DIRECTORY / "core_files.txt",
                        help="Snapshot of the derived core files")
    parser.add_argument("--update-core", action="store_true",
                        help="Regenerate the snapshot from the seeds and baseline, then check it")
    parser.add_argument("--report", type=Path, help="Write the complete deterministic JSON inventory")
    args = parser.parse_args(argv)
    try:
        root = args.root.resolve()
        entries = baseline_edges(json.loads(args.baseline.read_text()))
        stop = {(edge["source"], edge["target"]) for edge in entries}
        snapshot = read_core(args.core) if args.core.is_file() else set()
        if args.update_core:
            snapshot = {name for name in snapshot if (root / name).is_file()}
            report = scan(root, snapshot, stop)
            derived = set(report["core_files"])
            for name in sorted(derived - snapshot):
                print(f"core snapshot: + {name}")
            for name in sorted(snapshot - derived):
                print(f"core snapshot: - {name}")
            write_core(args.core, derived)
            snapshot = derived
        elif not args.core.is_file():
            raise ValueError(f"Missing core snapshot {args.core}; run with --update-core")
        report = scan(root, snapshot, stop)
        new, stale = check_baseline(report["violations"], entries)
        derived = set(report["core_files"])
        report["new_violations"], report["stale_baseline"] = new, stale
        report["new_core_files"] = sorted(derived - snapshot)
        report["stale_core_files"] = sorted(snapshot - derived)
        if args.report:
            args.report.parent.mkdir(parents=True, exist_ok=True)
            args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        edges = report["header_includes"] + report["implementation_includes"]
        counts = {status: sum(edge["status"] == status for edge in edges)
                  for status in ("internal", "generated", "unresolved", "macro", "external")}
        print(f"Scanned {len(report['files'])} files; {len(report['header_includes'])} header includes; "
              f"{len(report['implementation_includes'])} implementation includes.")
        print(f"Include status: {counts}; {len(report['cycles'])} file-level cycles; "
              f"{len(report['runtime_resources'])} literal runtime resource references.")
        print(f"Core: {len(snapshot)} files in snapshot; {len(report['new_core_files'])} new; "
              f"{len(report['stale_core_files'])} stale.")
        print(f"Boundary: {len(report['violations'])} existing violations; {len(new)} new; "
              f"{len(stale)} stale baseline entries.")
        for edge in new:
            print(f"NEW {edge['rule']}: {edge['source']}:{edge['line']} -> {edge['target']}", file=sys.stderr)
        for edge in stale:
            print(f"STALE {edge['rule']}: {edge['source']} -> {edge['target']} "
                  "(remove baseline entry)", file=sys.stderr)
        for label, names, hint in (
            ("NEW", report["new_core_files"], "baseline the crossing edge, or accept with --update-core"),
            ("STALE", report["stale_core_files"], "run --update-core"),
        ):
            for name in names[:LISTED_FILES]:
                print(f"{label} core file: {name} ({hint})", file=sys.stderr)
            if len(names) > LISTED_FILES:
                print(f"... and {len(names) - LISTED_FILES} more {label.lower()} core files", file=sys.stderr)
        failed = new or stale or report["new_core_files"] or report["stale_core_files"]
        return 1 if failed else 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Dependency check error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
