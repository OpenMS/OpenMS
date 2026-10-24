#!/usr/bin/env python3
"""Reject a reintroduced Xerces dependency in product sources and linked OpenMS."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--build-dir', type=Path)
    args = parser.parse_args()
    source = args.source_dir
    forbidden = re.compile(r'#\s*include\s*[<"]xercesc/|\bxercesc\s*::|\bXERCES_CPP_NAMESPACE')
    failures = []
    manifest = json.loads((source / 'vcpkg.json').read_text())
    names = [item if isinstance(item, str) else item['name'] for item in manifest['dependencies']]
    if 'xerces-c' in names:
        failures.append('vcpkg.json: xerces-c dependency')
    cmake_dependency = re.compile(r'\bXercesC\s*::|\bfind_(?:package|dependency)\s*\(\s*"?XercesC\b', re.IGNORECASE)
    package_dependency = re.compile(r'\b(?:lib)?xerces-c(?:[0-9][A-Za-z0-9_.-]*)?(?:[-_.][A-Za-z0-9]+)*\b', re.IGNORECASE)
    excluded = {'.git', 'CMakeFiles', '_deps', 'node_modules', '.venv'}
    excluded_root = {'build', 'vcpkg', 'vcpkg_installed'}
    excluded_paths = {('src', 'openms', 'extern'), ('src', 'openms', 'thirdparty'),
                      ('tools', 'xml-parser-benchmark')}
    for directory, directories, filenames in os.walk(source):
        parent = Path(directory)
        # Prune before descent: CI can contain very large dependency/build trees.
        directories[:] = [name for name in directories if name not in excluded
                          and not (parent == source and name in excluded_root)
                          and (parent / name).relative_to(source).parts not in excluded_paths
                          and not name.startswith('cmake-build-')
                          and not (parent / name / 'CMakeCache.txt').is_file()]
        for name in filenames:
            path = parent / name
            relative = path.relative_to(source)
            is_product_source = relative.parts[0] == 'src' and path.suffix in ('.h', '.hpp', '.cpp') and 'data' not in relative.parts
            is_cmake = name == 'CMakeLists.txt' or path.suffix == '.cmake'
            is_package_recipe = (name == 'Dockerfile' or name.startswith('Dockerfile.')
                                 or (relative.parts[0] == '.github' and path.suffix in ('.yml', '.yaml'))
                                 or (relative.parts[:2] == ('tools', 'ci') and path.suffix == '.sh')
                                 or relative.as_posix() == 'src/pyOpenMS/pyproject.toml')
            if not (is_product_source or is_cmake or is_package_recipe):
                continue
            content = path.read_text(encoding='utf-8', errors='replace')
            if is_product_source and forbidden.search(content):
                failures.append(str(relative))
            # Match dependency commands/packages, not removal notes or disable flags.
            content = re.sub(r'#[^\n]*', '', content)
            if ((is_cmake and cmake_dependency.search(content))
                    or (is_package_recipe and package_dependency.search(content))):
                failures.append(str(relative))
    if args.build_dir:
        libraries = set()
        for pattern in ('**/libOpenMS*.dylib', '**/libOpenMS.so*', '**/OpenMS.dll'):
            libraries.update(path.resolve() for path in args.build_dir.glob(pattern) if path.is_file())
        if not libraries:
            failures.append('No built OpenMS library found')
        for library in sorted(libraries):
            if sys.platform == 'darwin':
                command = ['otool', '-L', str(library)]
            elif sys.platform == 'win32':
                command = ['dumpbin', '/DEPENDENTS', str(library)]
            else:
                command = ['ldd', str(library)]
            try:
                links = subprocess.check_output(command, text=True, stderr=subprocess.STDOUT)
            except (OSError, subprocess.CalledProcessError, UnicodeError) as error:
                failures.append(f'{library}: could not inspect dependencies with {command[0]}: {error}')
                continue
            if re.search(r'\b(?:lib)?xerces[-_A-Za-z0-9.]*\.(?:dylib|so(?:\.\d+)*|dll)\b', links, re.IGNORECASE):
                failures.append(str(library) + ': linked Xerces dependency')
    if failures:
        parser.exit(1, 'Xerces removal check failed:\n' + '\n'.join(failures) + '\n')
    print('No Xerces dependency in product sources, manifest or checked OpenMS libraries.')


if __name__ == '__main__':
    main()
