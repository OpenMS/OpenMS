import pip
from pip._internal.req.req_file import parse_requirements
from importlib.metadata import version, PackageNotFoundError
from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version
import argparse
import sys

def check_dependencies(requirement_file_name):
    """
    Checks to see if the python dependencies are fulfilled.
    If check passes return 0. Otherwise print error and return 1.
    """
    requirements = list(parse_requirements(requirement_file_name, session=False))
    print("Found {} requirements".format(len(requirements)))

    for req in requirements:
        requirement_str = str(req.requirement)

        try:
            parsed_req = Requirement(requirement_str)
            installed_version = version(parsed_req.name)
            specifier_set = SpecifierSet(str(parsed_req.specifier))

            if installed_version in specifier_set:
                print(f" + {parsed_req.name}=={installed_version} is installed (required: {parsed_req})")
            else:
                print(f" - {parsed_req.name}=={installed_version} is installed but does not match the requirement {parsed_req}.")
                return 1
        except PackageNotFoundError:
            print(f"{parsed_req.name} is not installed.")
            return 1
        except Exception as e:
            print(f"An error occurred: {e}")
            return 1

    return 0


def main():
    parser = argparse.ArgumentParser(description="Check if all modules (with correct version) required for pyOpenMS are installed.")
    parser.add_argument('filename', help="The path where to find `requirements.txt`")
    
    args = parser.parse_args()
    
    if not args.filename:
        print("Error: Filename providing a `requirements.txt` is required.")
        sys.exit(1)
    
    ret = check_dependencies(args.filename)
    sys.exit(ret)

if __name__ == "__main__":
    main()
