"""
Check if all Python modules required for pyOpenMS build are installed.

Reads build dependencies from pyproject.toml [build-system].requires section
using a simple regex parser (no external TOML library needed).
"""

import argparse
import re
import sys
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet


def parse_build_requires(content: str) -> list:
    """Extract [build-system].requires from TOML content using regex.

    This is intentionally simple - it only handles the specific TOML structure
    used in pyproject.toml, not arbitrary TOML files.
    """
    match = re.search(
        r'\[build-system\].*?requires\s*=\s*\[(.*?)\]',
        content,
        re.DOTALL
    )
    if not match:
        return []
    return re.findall(r'"([^"]*)"', match.group(1))


def check_dependencies(pyproject_path: str) -> int:
    """
    Checks if the Python build dependencies are fulfilled.

    Args:
        pyproject_path: Path to pyproject.toml file

    Returns:
        0 if all dependencies satisfied, 1 otherwise
    """
    pyproject_file = Path(pyproject_path)

    if not pyproject_file.exists():
        print(f"Error: {pyproject_path} not found")
        return 1

    content = pyproject_file.read_text()
    build_deps = parse_build_requires(content)

    if not build_deps:
        print("Warning: No build dependencies found in [build-system].requires")
        return 0

    print(f"Found {len(build_deps)} build requirements")

    for dep_str in build_deps:
        try:
            parsed_req = Requirement(dep_str)
            installed_ver = version(parsed_req.name)
            specifier_set = SpecifierSet(str(parsed_req.specifier))

            if not parsed_req.specifier or installed_ver in specifier_set:
                print(f" + {parsed_req.name}=={installed_ver} is installed (required: {parsed_req})")
            else:
                print(f" - {parsed_req.name}=={installed_ver} is installed but does not match {parsed_req}")
                return 1
        except PackageNotFoundError:
            print(f" - {parsed_req.name} is not installed (required: {parsed_req})")
            return 1
        except Exception as e:
            print(f"Error checking {dep_str}: {e}")
            return 1

    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Check if all modules required for pyOpenMS build are installed."
    )
    parser.add_argument(
        'pyproject_path',
        nargs='?',
        default='pyproject.toml',
        help="Path to pyproject.toml (default: pyproject.toml)"
    )

    args = parser.parse_args()
    ret = check_dependencies(args.pyproject_path)
    sys.exit(ret)


if __name__ == "__main__":
    main()
