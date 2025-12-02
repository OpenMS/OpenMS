"""
Check if all Python modules required for pyOpenMS build are installed.

Reads build dependencies from pyproject.toml [build-system].requires section.
"""

import argparse
import sys
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

try:
    import tomllib
except ImportError:
    import tomli as tomllib

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet


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

    with open(pyproject_file, "rb") as f:
        pyproject = tomllib.load(f)

    # Get build dependencies from [build-system].requires (PEP 517 standard)
    build_system = pyproject.get("build-system", {})
    build_deps = build_system.get("requires", [])

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
