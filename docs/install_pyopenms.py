#!/usr/bin/env python
"""
Script to install pyopenms with the version from ReadTheDocs environment variable.
This is used by ReadTheDocs build process to ensure the correct version is installed.
"""
import subprocess
import sys
import os
import re

def get_version_from_rtd():
    """
    Extract version from READTHEDOCS_VERSION_NAME environment variable.
    
    The variable will be:
    - "latest" for latest builds
    - "release/X.Y.Z" for release branches
    - Tag names for tag builds
    
    Returns the version to install (e.g., "3.5.0" or "latest")
    """
    rtd_version = os.environ.get('READTHEDOCS_VERSION_NAME', '')
    
    if not rtd_version:
        print("Warning: READTHEDOCS_VERSION_NAME not set, defaulting to 'latest'")
        return 'latest'
    
    print(f"READTHEDOCS_VERSION_NAME: {rtd_version}")
    
    # If it's "latest", install latest version
    if rtd_version == "latest":
        return 'latest'
    
    # If it's a release branch like "release/3.5.0" or "Release/3.5.0", extract the version
    if rtd_version.lower().startswith('release/'):
        version = rtd_version.split('/', 1)[1]
        print(f"Extracted version from release branch: {version}")
        return version
    
    # If it looks like a version tag (e.g., "v3.5.0", "3.5.0"), use it directly
    # Strip leading 'v' if present
    version = rtd_version[1:] if rtd_version.startswith('v') else rtd_version
    
    # Validate that it looks like a version number (basic check)
    # Accept versions like "3.5.0", "3.5.0rc1", "3.5.0.dev1", etc.
    if not re.match(r'^\d+\.\d+', version):
        print(f"Warning: '{version}' does not look like a valid version, using 'latest'")
        return 'latest'
    
    print(f"Using version from tag: {version}")
    return version

def install_pyopenms(version):
    """Install pyopenms with specific version or latest"""
    
    # Build pip install command
    cmd = [
        sys.executable, "-m", "pip", "install",
        "--extra-index-url", "https://pypi.cs.uni-tuebingen.de/simple/"
    ]
    
    if version == 'latest':
        print("Installing latest version of pyopenms (including pre-releases)")
        cmd.append("--pre")  # Allow pre-release versions only for latest
        cmd.append("pyopenms")
    else:
        print(f"Installing pyopenms=={version}")
        cmd.append(f"pyopenms=={version}")
    
    subprocess.check_call(cmd)
    
    print(f"Successfully installed pyopenms")

if __name__ == "__main__":
    try:
        version = get_version_from_rtd()
        install_pyopenms(version)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
