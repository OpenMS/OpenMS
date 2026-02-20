#!/usr/bin/env python3
"""
Mechanically transform .cpp files to use C++23 `import std;`.

Replaces unconditional standard library #include directives, converts
macro-providing headers to their C-style equivalents, and inserts
`import std;` after the last remaining #include.

Usage:
    python3 tools/import_std_transform.py [--dry-run] [--verbose] [path...]

    --dry-run   Report changes without modifying files
    --verbose   Print details per file
    path...     Files or directories to process (default: src/)
"""

import argparse
import os
import re
from typing import NamedTuple

# ---------------------------------------------------------------------------
# Standard library header whitelist
# ---------------------------------------------------------------------------

CPP_HEADERS = frozenset("""
    algorithm any array atomic barrier bit bitset
    charconv chrono codecvt compare complex concepts condition_variable coroutine
    deque
    exception execution expected
    filesystem flat_map flat_set format forward_list fstream functional future
    generator
    initializer_list iomanip ios iosfwd iostream istream iterator
    latch limits list locale
    map mdspan memory memory_resource mutex
    new numbers numeric
    optional ostream
    print
    queue
    random ranges ratio regex
    scoped_allocator semaphore set shared_mutex source_location span spanstream
    sstream stack stacktrace stdexcept stdfloat stop_token streambuf string
    string_view strstream syncstream system_error
    thread tuple type_traits typeindex typeinfo
    unordered_map unordered_set utility
    valarray variant vector version
""".split())

C_COMPAT_HEADERS = frozenset("""
    cassert cctype cerrno cfenv cfloat cinttypes ciso646 climits clocale cmath
    csetjmp csignal cstdarg cstddef cstdint cstdio cstdlib cstring ctime
    cuchar cwchar cwctype
""".split())

ALL_STD_HEADERS = CPP_HEADERS | C_COMPAT_HEADERS

# Headers whose macros are NOT fully exported by `import std;`.
# These must be converted to their C-style <xxx.h> equivalents.
MACRO_HEADER_MAP = {
    "cassert":  "assert.h",
    "cerrno":   "errno.h",
    "cfenv":    "fenv.h",
    "cfloat":   "float.h",
    "climits":  "limits.h",
    "clocale":  "locale.h",
    "csetjmp":  "setjmp.h",
    "csignal":  "signal.h",
    "cstdarg":  "stdarg.h",
    "cstdint":  "stdint.h",
    "cstdio":   "stdio.h",
}

# Headers that inject functions into the global namespace which `import std;`
# does NOT do. These must be preserved when the file uses unqualified calls.
# `<cmath>` is the main case: it injects sqrt, log, pow, etc. into the
# global namespace, while `import std;` only provides std::sqrt, etc.
GLOBAL_NS_HEADERS = frozenset(["cmath", "cstdlib"])

# Unqualified C math function names that require <cmath> at global scope.
UNQUALIFIED_MATH_FUNCS = frozenset([
    "sqrt", "log", "log2", "log10", "fabs", "pow", "exp", "ceil", "floor",
    "round", "abs", "sin", "cos", "tan", "atan", "atan2", "asin", "acos",
    "lgamma", "cbrt", "hypot", "erf", "erfc", "tgamma", "remainder",
    "copysign", "fdim", "fmax", "fmin", "fmod", "trunc", "nearbyint",
    "nextafter", "isnan", "isinf", "isfinite", "fpclassify", "signbit",
])

RE_UNQUALIFIED_MATH = re.compile(
    r'\b(' + '|'.join(UNQUALIFIED_MATH_FUNCS) + r')\s*\('
)

# ---------------------------------------------------------------------------
# Directories to exclude
# ---------------------------------------------------------------------------

EXCLUDE_DIRS = {"src/openms/extern", "src/pyOpenMS", ".venv", "__pycache__"}

# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# Matches: #include <header_name>  (angle-bracket includes only)
RE_INCLUDE = re.compile(r'^\s*#\s*include\s*<([^>]+)>\s*(?://.*)?$')

# Matches any #include directive (angle-bracket or quoted)
RE_ANY_INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]')

# Preprocessor directives that affect conditional nesting
RE_PP_IF     = re.compile(r'^\s*#\s*(?:if|ifdef|ifndef)\b')
RE_PP_ENDIF  = re.compile(r'^\s*#\s*endif\b')
RE_PP_ELIF   = re.compile(r'^\s*#\s*(?:elif|else)\b')

# Copyright header comment block detection (// style)
RE_COMMENT_LINE = re.compile(r'^\s*//')

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

class IncludeInfo(NamedTuple):
    """Information about an #include directive."""
    line_idx: int        # 0-based line index
    header: str          # header name (e.g. "vector", "cassert")
    is_std: bool         # whether it's a standard library header
    is_conditional: bool # whether it's inside a conditional block
    is_macro_header: bool  # whether it needs C-style conversion


class FileResult(NamedTuple):
    """Summary of changes for a single file."""
    path: str
    removed: list       # headers removed
    converted: list     # headers converted (macro headers)
    preserved: list     # headers preserved (conditional)
    modified: bool


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def file_has_unqualified_math(lines: list[str]) -> bool:
    """
    Check if a file uses unqualified C math function calls (e.g. sqrt()
    instead of std::sqrt()). Such files need to keep #include <cmath>.
    """
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('//') or stripped.startswith('#') or stripped.startswith('/*'):
            continue
        for m in RE_UNQUALIFIED_MATH.finditer(line):
            start = m.start()
            prefix = line[:start]
            rprefix = prefix.rstrip()
            # Skip if preceded by qualifier (std::, boost::, etc.)
            if rprefix.endswith('::') or rprefix.endswith('.') or rprefix.endswith('->'):
                continue
            # Skip if inside string literal (rough check)
            quote_count = prefix.count('"') - prefix.count('\\"')
            if quote_count % 2 == 1:
                continue
            return True
    return False


def is_in_conditional_block(lines: list[str], target_idx: int) -> bool:
    """
    Walk backward from target_idx to determine if the line is inside
    an #if/#ifdef/#ifndef block.

    Track nesting depth going backward:
    - #endif increments depth (we're entering a completed block)
    - #if/#ifdef/#ifndef at depth 0 means we're inside a conditional
    - #elif/#else at depth 0 means we're in a conditional branch
    """
    depth = 0
    for i in range(target_idx - 1, -1, -1):
        line = lines[i]
        if RE_PP_ENDIF.match(line):
            depth += 1
        elif RE_PP_IF.match(line):
            if depth == 0:
                return True
            depth -= 1
        elif RE_PP_ELIF.match(line):
            if depth == 0:
                return True
    return False


def find_copyright_end(lines: list[str]) -> int:
    """
    Find the line index after the copyright header block.
    The copyright block is a contiguous run of // comment lines
    at the start of the file (possibly with blank lines interspersed
    within the block).

    Returns the index of the first non-comment, non-blank line,
    or 0 if no copyright block found.
    """
    in_comment_block = False
    last_comment_line = -1

    for i, line in enumerate(lines):
        stripped = line.strip()
        if RE_COMMENT_LINE.match(line):
            in_comment_block = True
            last_comment_line = i
        elif stripped == '' and in_comment_block:
            # Blank line within/after comment block — continue
            continue
        else:
            break

    if last_comment_line >= 0:
        return last_comment_line + 1
    return 0


def transform_file(filepath: str, dry_run: bool = False) -> FileResult | None:
    """
    Transform a single .cpp file. Returns FileResult or None if skipped.
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    # Find all #include directives and classify them
    includes: list[IncludeInfo] = []
    has_std_include = False

    for idx, line in enumerate(lines):
        m = RE_INCLUDE.match(line)
        if not m:
            continue
        header = m.group(1)
        is_std = header in ALL_STD_HEADERS
        if is_std:
            has_std_include = True
            conditional = is_in_conditional_block(lines, idx)
            is_macro = header in MACRO_HEADER_MAP
            includes.append(IncludeInfo(
                line_idx=idx,
                header=header,
                is_std=True,
                is_conditional=conditional,
                is_macro_header=is_macro,
            ))
        else:
            includes.append(IncludeInfo(
                line_idx=idx,
                header=header,
                is_std=False,
                is_conditional=False,
                is_macro_header=False,
            ))

    if not has_std_include:
        return None  # No standard library includes — skip

    # Check if the file uses unqualified C math functions.
    # If so, <cmath> (and <cstdlib>) must be preserved to keep those
    # functions available in the global namespace.
    needs_global_ns = file_has_unqualified_math(lines)

    # Classify std includes into actions
    to_remove = []
    to_convert = []
    to_preserve = []

    for inc in includes:
        if not inc.is_std:
            continue
        if inc.is_conditional:
            to_preserve.append(inc)
        elif inc.is_macro_header:
            to_convert.append(inc)
        elif needs_global_ns and inc.header in GLOBAL_NS_HEADERS:
            to_preserve.append(inc)
        else:
            to_remove.append(inc)

    if not to_remove and not to_convert:
        # All std includes are conditional — nothing to do
        return FileResult(
            path=filepath,
            removed=[],
            converted=[],
            preserved=[inc.header for inc in to_preserve],
            modified=False,
        )

    # Build the result for reporting
    result = FileResult(
        path=filepath,
        removed=[inc.header for inc in to_remove],
        converted=[inc.header for inc in to_convert],
        preserved=[inc.header for inc in to_preserve],
        modified=True,
    )

    if dry_run:
        return result

    # --- Apply transformations ---

    # Collect line indices to remove (unconditional std includes that
    # are NOT macro headers)
    lines_to_remove = set()
    for inc in to_remove:
        lines_to_remove.add(inc.line_idx)

    # Convert macro headers in-place
    for inc in to_convert:
        c_header = MACRO_HEADER_MAP[inc.header]
        # Replace the line, preserving any leading whitespace
        old_line = lines[inc.line_idx]
        indent = old_line[:len(old_line) - len(old_line.lstrip())]
        lines[inc.line_idx] = f'{indent}#include <{c_header}>\n'

    # Remove lines marked for removal
    new_lines = []
    for idx, line in enumerate(lines):
        if idx in lines_to_remove:
            continue
        new_lines.append(line)

    # Find insertion point for `import std;`
    # Place it after the VERY LAST unconditional #include directive of any
    # kind (angle-bracket or quoted). This ensures `import std;` never
    # appears before a textual #include that might transitively pull in
    # standard library headers — all three compilers (GCC, Clang, MSVC)
    # only support #include-before-import, not the reverse.
    last_include_idx = -1
    for idx, line in enumerate(new_lines):
        if RE_ANY_INCLUDE.match(line) and not is_in_conditional_block(new_lines, idx):
            last_include_idx = idx

    if last_include_idx >= 0:
        insert_idx = last_include_idx + 1
        # Check if next line is blank; if not, add one
        import_block = []
        # Add blank line before import std; if the previous line isn't blank
        if insert_idx > 0 and new_lines[insert_idx - 1].strip() != '':
            import_block.append('\n')
        import_block.append('import std;\n')
        # Add blank line after if the next line isn't blank and exists
        if insert_idx < len(new_lines) and new_lines[insert_idx].strip() != '':
            import_block.append('\n')

        for i, bl in enumerate(import_block):
            new_lines.insert(insert_idx + i, bl)
    else:
        # No remaining #include directives at all — place after copyright header
        copyright_end = find_copyright_end(new_lines)
        insert_idx = copyright_end
        import_block = []
        # Ensure blank line before
        if insert_idx > 0 and new_lines[insert_idx - 1].strip() != '':
            import_block.append('\n')
        import_block.append('import std;\n')
        # Ensure blank line after
        if insert_idx < len(new_lines) and new_lines[insert_idx].strip() != '':
            import_block.append('\n')

        for i, bl in enumerate(import_block):
            new_lines.insert(insert_idx + i, bl)

    # Clean up: remove excessive blank lines (more than 2 consecutive)
    cleaned = []
    blank_count = 0
    for line in new_lines:
        if line.strip() == '':
            blank_count += 1
            if blank_count <= 2:
                cleaned.append(line)
        else:
            blank_count = 0
            cleaned.append(line)

    # Write the file
    with open(filepath, 'w', encoding='utf-8', newline='') as f:
        f.writelines(cleaned)

    return result


# ---------------------------------------------------------------------------
# File collection
# ---------------------------------------------------------------------------

def should_exclude(filepath: str, base_dir: str) -> bool:
    """Check if a file path should be excluded."""
    rel = os.path.relpath(filepath, base_dir)
    for excl in EXCLUDE_DIRS:
        if rel.startswith(excl + os.sep) or rel.startswith(excl + '/') or rel == excl:
            return True
    # Also check normalized with forward slashes
    rel_fwd = rel.replace(os.sep, '/')
    for excl in EXCLUDE_DIRS:
        if rel_fwd.startswith(excl + '/') or rel_fwd == excl:
            return True
    return False


def collect_cpp_files(paths: list[str], base_dir: str) -> list[str]:
    """Collect all .cpp files from the given paths, respecting exclusions."""
    cpp_files = []
    for p in paths:
        p = os.path.abspath(p)
        if os.path.isfile(p):
            if p.endswith('.cpp') and not should_exclude(p, base_dir):
                cpp_files.append(p)
        elif os.path.isdir(p):
            for root, dirs, files in os.walk(p):
                # Prune excluded directories
                dirs[:] = [
                    d for d in dirs
                    if d not in {'__pycache__', '.venv'}
                    and not should_exclude(os.path.join(root, d), base_dir)
                ]
                for fname in sorted(files):
                    if fname.endswith('.cpp'):
                        fpath = os.path.join(root, fname)
                        if not should_exclude(fpath, base_dir):
                            cpp_files.append(fpath)
    return sorted(set(cpp_files))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Transform .cpp files to use C++23 import std;'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Report changes without modifying files',
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print details per file',
    )
    parser.add_argument(
        'paths',
        nargs='*',
        help='Files or directories to process (default: src/)',
    )
    args = parser.parse_args()

    # Determine project root (script lives in tools/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)

    # Default paths
    if not args.paths:
        search_paths = [os.path.join(project_root, 'src')]
    else:
        search_paths = args.paths

    cpp_files = collect_cpp_files(search_paths, project_root)

    mode_label = "[DRY RUN] " if args.dry_run else ""
    print(f"{mode_label}Processing {len(cpp_files)} .cpp files...\n")

    # Counters
    files_processed = 0
    files_modified = 0
    files_skipped = 0
    total_removed = 0
    total_converted = 0
    total_preserved = 0

    for fpath in cpp_files:
        files_processed += 1
        rel_path = os.path.relpath(fpath, project_root)

        result = transform_file(fpath, dry_run=args.dry_run)

        if result is None:
            files_skipped += 1
            if args.verbose:
                print(f"  {rel_path}: (skipped, no std includes)")
            continue

        if not result.modified:
            files_skipped += 1
            if args.verbose:
                print(f"  {rel_path}: (skipped, all conditional)")
            continue

        files_modified += 1
        total_removed += len(result.removed)
        total_converted += len(result.converted)
        total_preserved += len(result.preserved)

        if args.verbose:
            print(f"  {rel_path}:")
            if result.removed:
                print(f"    Remove: {', '.join(result.removed)}")
            if result.converted:
                for h in result.converted:
                    print(f"    Convert: <{h}> -> <{MACRO_HEADER_MAP[h]}>")
            if result.preserved:
                print(f"    Preserve (conditional): {', '.join(result.preserved)}")
            print()

    # Summary
    print("=" * 60)
    summary_label = f"Summary ({mode_label.strip()}):" if args.dry_run else "Summary:"
    print(summary_label)
    print(f"  Files processed:    {files_processed}")
    print(f"  Files modified:     {files_modified}")
    print(f"  Files skipped:      {files_skipped}")
    print(f"  Includes removed:   {total_removed}")
    print(f"  Includes converted: {total_converted} (macro headers -> C-style)")
    print(f"  Includes preserved: {total_preserved} (inside #ifdef blocks)")


if __name__ == '__main__':
    main()
