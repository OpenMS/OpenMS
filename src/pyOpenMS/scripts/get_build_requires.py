#!/usr/bin/env python
"""Output build dependencies from pyproject.toml for CI workflows.

Uses a simple regex parser - no external dependencies (tomli/tomllib) needed.
This works on any Python 3.9+ without requiring conditional package installs.
"""

import re
import sys
from pathlib import Path


def parse_build_requires(content: str) -> list:
    """Extract [build-system].requires from TOML content using regex.

    This is intentionally simple - it only handles the specific TOML structure
    we use in pyproject.toml, not arbitrary TOML files.
    """
    # Find [build-system] section and its requires array
    # Match: [build-system] ... requires = [ ... ]
    match = re.search(
        r'\[build-system\].*?requires\s*=\s*\[(.*?)\]',
        content,
        re.DOTALL
    )
    if not match:
        return []

    # Extract all double-quoted strings from the array
    array_content = match.group(1)
    return re.findall(r'"([^"]*)"', array_content)


def main():
    # Find pyproject.toml relative to this script
    script_dir = Path(__file__).parent
    pyproject_path = script_dir.parent / "pyproject.toml"

    # Allow override via command line argument
    if len(sys.argv) > 1:
        pyproject_path = Path(sys.argv[1])

    if not pyproject_path.exists():
        print(f"Error: {pyproject_path} not found", file=sys.stderr)
        sys.exit(1)

    content = pyproject_path.read_text()
    requires = parse_build_requires(content)

    # Output space-separated unquoted strings for uv consumption
    print(" ".join(requires))


if __name__ == "__main__":
    main()
