#!/usr/bin/env python
"""Output build dependencies from pyproject.toml for CI workflows."""

import sys
from pathlib import Path

try:
    import tomllib
except ImportError:
    import tomli as tomllib


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

    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)

    requires = data.get("build-system", {}).get("requires", [])
    # Output space-separated quoted strings for shell consumption
    print(" ".join(repr(dep) for dep in requires))


if __name__ == "__main__":
    main()
