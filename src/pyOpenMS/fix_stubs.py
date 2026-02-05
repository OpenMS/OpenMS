#!/usr/bin/env python3
"""Post-process nanobind-generated .pyi stubs for pyOpenMS.

Fixes issues that nanobind's stubgen cannot handle:
1. Python keywords used as parameter names (from C++ parameter names)
2. Malformed NDArray type annotations
"""

import re
import sys
from pathlib import Path

# Python keywords that appear as C++ parameter names
PYTHON_KEYWORDS = {
    "is", "from", "in", "and", "or", "not", "class", "def", "return",
    "import", "with", "as", "for", "while", "if", "else", "elif",
    "try", "except", "finally", "raise", "pass", "break", "continue",
    "lambda", "yield", "global", "nonlocal", "del", "assert",
    "True", "False", "None",
}


def fix_keyword_params(line: str) -> str:
    """Rename Python keyword parameter names by appending underscore."""
    # Match 'keyword:' or 'keyword =' in function signatures
    for kw in PYTHON_KEYWORDS:
        # Match keyword as a parameter name (preceded by comma/paren and space, followed by colon or =)
        line = re.sub(
            rf'(\b){kw}(\s*[:=])',
            rf'\g<1>{kw}_\2',
            line,
        )
    return line


def fix_ndarray_annotations(line: str) -> str:
    """Fix malformed NDArray annotations from nanobind stubgen.

    nanobind emits: Annotated[NDArray, dict(numpy.float64[Any, 2)]]
    Should be:      numpy.ndarray
    """
    line = re.sub(
        r'Annotated\[NDArray,\s*dict\(numpy\.\w+\[[^\]]*\)\]\]',
        'numpy.ndarray',
        line,
    )
    return line


def fix_singleton_types(line: str) -> str:
    """Fix CallableSingleton type annotations from singleton addon.

    nanobind emits: ElementDB: pyopenms.addons.singletons.make_singleton_callable.<locals>.CallableSingleton = ...
    Should be a variable with the actual class type.
    """
    line = re.sub(
        r': pyopenms\.addons\.singletons\.make_singleton_callable\.<locals>\.CallableSingleton = \.\.\.',
        ': ... # singleton (use getInstance())',
        line,
    )
    return line


def fix_stub_file(path: Path) -> bool:
    """Fix a single .pyi file. Returns True if modified."""
    content = path.read_text()
    lines = content.splitlines(keepends=True)
    modified = False

    new_lines = []
    for line in lines:
        new_line = line
        new_line = fix_keyword_params(new_line)
        new_line = fix_ndarray_annotations(new_line)
        new_line = fix_singleton_types(new_line)
        if new_line != line:
            modified = True
        new_lines.append(new_line)

    if modified:
        path.write_text("".join(new_lines))
    return modified


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <pyopenms_dir>", file=sys.stderr)
        sys.exit(1)

    pyopenms_dir = Path(sys.argv[1])
    if not pyopenms_dir.is_dir():
        print(f"Error: {pyopenms_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    fixed = 0
    for pyi_file in sorted(pyopenms_dir.glob("**/*.pyi")):
        if fix_stub_file(pyi_file):
            fixed += 1
            print(f"  Fixed: {pyi_file.name}")

    print(f"Fixed {fixed} stub files")


if __name__ == "__main__":
    main()
