#!/usr/bin/env python3
"""Post-process nanobind-generated .pyi stubs for pyOpenMS.

Fixes issues that nanobind's stubgen cannot handle:
1. Python keywords used as parameter names (from C++ parameter names)
2. Malformed NDArray type annotations
3. nanobind emits ``None_`` instead of ``None`` for return types
4. RST ``.. code-block::`` directives need blank line + indented body
5. C++ type annotations leaked as string literals
6. Duplicate identical overloads from const/non-const C++ methods
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
    """Rename Python keyword parameter names by appending underscore.

    Uses negative lookbehind ``(?<!:)`` to avoid mangling RST directives
    like ``:return:`` or ``:param from:``.
    """
    for kw in PYTHON_KEYWORDS:
        line = re.sub(
            rf'(?<!:)(\b){kw}(\s*[:=])',
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


def fix_none_type(line: str) -> str:
    r"""Replace nanobind's ``None_`` with ``None`` in type annotations.

    nanobind stubgen emits ``None_`` (with trailing underscore) instead of
    the standard Python ``None`` type.  This is invalid in .pyi stubs since
    ``None_`` is never imported or defined.
    """
    return re.sub(r'\bNone_\b', 'None', line)


def fix_cpp_type_annotations(line: str) -> str:
    """Replace C++ type string annotations with ``Any``.

    nanobind stubgen emits C++ types as string-literal annotations when it
    cannot resolve a type to a Python equivalent, e.g.::

        def foo(self, arg: "OpenMS::SomeType") -> None: ...

    These are replaced with ``Any`` since the C++ types are not valid Python.
    """
    return re.sub(
        r'"(?:OpenMS::|std::|__gnu_cxx::)[^"]*"',
        'Any',
        line,
    )


def fix_capitalized_builtins(line: str) -> str:
    """Convert capitalized typing forms to PEP 585 lowercase builtins.

    nanobind's built-in type casters emit ``Set[T]``, ``List[T]``, etc.
    which require ``from typing import Set, List, ...``.  The lowercase
    ``set[T]``, ``list[T]`` forms work natively in Python 3.9+.
    """
    line = re.sub(r'\bList\[', 'list[', line)
    line = re.sub(r'\bDict\[', 'dict[', line)
    line = re.sub(r'\bSet\[', 'set[', line)
    line = re.sub(r'\bTuple\[', 'tuple[', line)
    line = re.sub(r'\bFrozenSet\[', 'frozenset[', line)
    return line


def fix_code_blocks(content: str) -> str:
    """Fix RST ``.. code-block::`` directives in docstrings.

    nanobind copies docstrings verbatim from C++ bindings where the code
    examples lack the blank line and extra indentation required by RST::

        Before (broken):
            .. code-block:: python
            exp = MSExperiment()

        After (valid RST):
            .. code-block:: python

                exp = MSExperiment()
    """
    lines = content.split('\n')
    result: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r'^(\s*)\.\. code-block::\s*\w+\s*$', line)
        if m:
            directive_indent = m.group(1)
            code_indent = directive_indent + '    '
            result.append(line)

            # Check if the next line is already blank (already formatted)
            next_is_blank = (i + 1 < len(lines) and lines[i + 1].strip() == '')
            if not next_is_blank:
                result.append('')  # insert required blank line after directive
            i += 1

            # Indent subsequent lines until end of docstring
            while i < len(lines):
                code_line = lines[i]
                # Stop at docstring closing quotes
                if code_line.strip() in ('"""', "'''"):
                    result.append(code_line)
                    i += 1
                    break
                # Blank lines pass through unchanged
                if code_line.strip() == '':
                    result.append('')
                    i += 1
                    continue
                # Only add indentation if content is not already indented
                # beyond the directive level
                stripped = code_line.lstrip()
                current_indent = len(code_line) - len(stripped)
                directive_indent_len = len(directive_indent)
                if current_indent <= directive_indent_len:
                    result.append(code_indent + stripped)
                else:
                    result.append(code_line)
                i += 1
        else:
            result.append(line)
            i += 1
    return '\n'.join(result)


def fix_duplicate_overloads(content: str) -> str:
    """Remove consecutive duplicate ``@overload`` definitions.

    nanobind may emit identical overloads for const/non-const C++ methods.
    This pass removes exact consecutive duplicates.
    """
    lines = content.split('\n')
    result: list[str] = []
    i = 0
    last_overload_block: str | None = None

    while i < len(lines):
        stripped = lines[i].strip()

        if stripped == '@overload':
            # Collect this overload block (until blank line or next @overload)
            block = [lines[i]]
            i += 1

            while i < len(lines):
                s = lines[i].strip()
                if s == '':
                    block.append(lines[i])
                    i += 1
                    break
                if s == '@overload':
                    break  # next overload starts, don't consume
                block.append(lines[i])
                i += 1

            block_text = '\n'.join(block)

            if block_text == last_overload_block:
                continue  # skip duplicate

            last_overload_block = block_text
            result.extend(block)
        else:
            if stripped:  # non-blank, non-overload line resets tracking
                last_overload_block = None
            result.append(lines[i])
            i += 1

    return '\n'.join(result)


def ensure_any_import(content: str) -> str:
    """Add ``Any`` to typing imports if ``Any`` is used but not imported."""
    if 'Any' not in content:
        return content

    # Already imported
    if re.search(r'from typing import.*\bAny\b', content):
        return content

    # Add to existing typing import
    new_content = re.sub(
        r'(from typing import )',
        r'\1Any, ',
        content,
        count=1,
    )
    if new_content != content:
        return new_content

    # No typing import exists — add one after the module docstring
    new_content = re.sub(
        r'("""[^"]*""")\n',
        r'\1\nfrom typing import Any\n',
        content,
        count=1,
    )
    return new_content


def fix_stub_file(path: Path) -> bool:
    """Fix a single .pyi file. Returns True if modified."""
    content = path.read_text()

    # Multi-line fixes (operate on full content)
    new_content = fix_code_blocks(content)

    # Line-by-line fixes (run before dedup since type replacements
    # like C++ types -> Any can create new duplicates)
    lines = new_content.splitlines(keepends=True)
    new_lines = []
    for line in lines:
        new_line = line
        new_line = fix_keyword_params(new_line)
        new_line = fix_ndarray_annotations(new_line)
        new_line = fix_singleton_types(new_line)
        new_line = fix_none_type(new_line)
        new_line = fix_cpp_type_annotations(new_line)
        new_line = fix_capitalized_builtins(new_line)
        new_lines.append(new_line)

    new_content = "".join(new_lines)

    # Dedup after line-level fixes (C++ type -> Any may create duplicates)
    new_content = fix_duplicate_overloads(new_content)

    # Ensure imports for added types
    new_content = ensure_any_import(new_content)

    modified = new_content != content
    if modified:
        path.write_text(new_content)
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
