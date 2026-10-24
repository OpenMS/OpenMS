"""Exercise the CI removal guard against source and binary regressions."""
import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('check_no_xerces', Path(__file__).resolve().parents[1] / 'check_no_xerces.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)


class RemovalGuardTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix='XERCES Replacement ')
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        for name in ['src/openms/CMakeLists.txt', 'cmake/cmake_findExternalLibs.cmake']:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('')
        (self.root / 'vcpkg.json').write_text('{"dependencies":["libxml2"]}')

    def run_guard(self, *args):
        with patch.object(sys, 'argv', ['guard', '--source-dir', str(self.root), *args]), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                check.main()
            except SystemExit as error:
                return error.code
        return 0

    def test_product_source_and_manifest_regressions(self):
        self.assertEqual(self.run_guard(), 0)
        header = self.root / 'src/consumer.h'
        header.write_text('#include <xercesc/dom/DOM.hpp>\n')
        self.assertEqual(self.run_guard(), 1)
        header.unlink()
        (self.root / 'vcpkg.json').write_text('{"dependencies":[{"name":"xerces-c"}]}')
        self.assertEqual(self.run_guard(), 1)

    def test_nested_cmake_dependencies(self):
        for name, content in [
            ('src/topp/CMakeLists.txt', 'find_package ( \n XercesC REQUIRED)'),
            ('cmake/nested/dependencies.cmake', 'target_link_libraries(tool XercesC::XercesC)'),
            ('CMakeLists.txt', 'FIND_PACKAGE(\nXercesC)'),
            ('cmake/dependencies.cmake', 'find_dependency ( "XercesC" REQUIRED)'),
        ]:
            with self.subTest(name=name):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
                self.assertEqual(self.run_guard(), 1)
                path.unlink()

    def test_removal_notes_and_vendored_cmake_are_allowed(self):
        (self.root / 'CMakeLists.txt').write_text(
            '# Removed find_package(XercesC)\nset(CMAKE_DISABLE_FIND_PACKAGE_XercesC ON)\n')
        for name in ['src/openms/extern/CMakeLists.txt',
                     'vcpkg/ports/xerces-c/portfile.cmake',
                     'tools/xml-parser-benchmark/CMakeLists.txt']:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('find_package(XercesC REQUIRED)')
        self.assertEqual(self.run_guard(), 0)

    def test_package_reintroduction_in_build_recipes(self):
        for name, content in [
            ('dockerfiles/Dockerfile', 'RUN apt-get install libxerces-c-dev'),
            ('.github/workflows/build.yml', 'run: brew install xerces-c'),
            ('.github/actions/build/action.yml', 'run: brew install xerces-c'),
            ('tools/ci/deps.sh', 'apt-get install libxerces-c3.2'),
            ('src/pyOpenMS/pyproject.toml', 'before-all = "brew install xerces-c"'),
        ]:
            with self.subTest(name=name):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
                self.assertEqual(self.run_guard(), 1)
                path.write_text('# Removed xerces-c and libxerces-c-dev\n'
                                'cmake -DCMAKE_DISABLE_FIND_PACKAGE_XercesC=ON\n')
                self.assertEqual(self.run_guard(), 0)
                path.unlink()

    def test_nested_source_directories_are_not_build_trees(self):
        for name in ['src/foo/build/CMakeLists.txt', 'src/foo/vcpkg/CMakeLists.txt']:
            with self.subTest(name=name):
                path = self.root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('find_package(XercesC REQUIRED)')
                self.assertEqual(self.run_guard(), 1)
                path.unlink()

    def test_generated_and_vendor_trees_are_pruned(self):
        for name in ['build', 'custom-build', 'vcpkg', 'src/openms/extern']:
            path = self.root / name
            path.mkdir(parents=True, exist_ok=True)
            (path / 'CMakeCache.txt').touch()
            (path / 'CMakeLists.txt').write_text('find_package(XercesC REQUIRED)')
        visited = []
        walk = check.os.walk

        def record_walk(*args, **kwargs):
            for item in walk(*args, **kwargs):
                visited.append(Path(item[0]).relative_to(self.root).as_posix())
                yield item

        with patch.object(check.os, 'walk', side_effect=record_walk):
            self.assertEqual(self.run_guard(), 0)
        for excluded in ['build', 'custom-build', 'vcpkg', 'src/openms/extern']:
            self.assertNotIn(excluded, visited)

    def test_non_utf8_source_has_no_decode_traceback(self):
        path = self.root / 'src/consumer.cpp'
        path.write_bytes(b'// Latin-1 caf\xe9\n')
        self.assertEqual(self.run_guard(), 0)
        path.write_bytes(b'// Latin-1 caf\xe9\n#include <xercesc/dom/DOM.hpp>\n')
        self.assertEqual(self.run_guard(), 1)

    def test_dependency_names_not_workspace_directory(self):
        build = self.root / 'build'
        build.mkdir()
        library = build / 'libOpenMS.dylib'
        library.touch()
        clean = str(library) + ':\n /opt/lib/libxml2.16.dylib (compatibility version 16.0)\n'
        with patch.object(check.subprocess, 'check_output', return_value=clean):
            self.assertEqual(self.run_guard('--build-dir', str(build)), 0)
        for dependency in ['libxerces-c-3.3.so', 'libxerces-c.dylib', 'xerces-c_3_3.dll']:
            with self.subTest(dependency=dependency), patch.object(check.subprocess, 'check_output', return_value=clean + '/lib/' + dependency):
                self.assertEqual(self.run_guard('--build-dir', str(build)), 1)

    def test_missing_build_output_fails(self):
        self.assertEqual(self.run_guard('--build-dir', str(self.root / 'missing')), 1)

    def test_checkout_parent_does_not_exclude_product_sources(self):
        original = self.root
        for parent in ['data', 'extern', 'thirdparty']:
            with self.subTest(parent=parent):
                self.root = original / parent / 'OpenMS'
                (self.root / 'src/openms').mkdir(parents=True)
                (self.root / 'cmake').mkdir()
                (self.root / 'src/openms/CMakeLists.txt').write_text('')
                (self.root / 'cmake/cmake_findExternalLibs.cmake').write_text('')
                (self.root / 'vcpkg.json').write_text('{"dependencies":["libxml2"]}')
                (self.root / 'src/consumer.h').write_text('#include <xercesc/dom/DOM.hpp>\n')
                self.assertEqual(self.run_guard(), 1)
        self.root = original

    def test_dependency_inspector_failure_is_reported(self):
        build = self.root / 'build'
        build.mkdir()
        (build / 'libOpenMS.dylib').touch()
        for error in [FileNotFoundError('missing dependency inspector'),
                      check.subprocess.CalledProcessError(1, ['inspector']),
                      UnicodeDecodeError('utf-8', b'\xff', 0, 1, 'invalid inspector output')]:
            with self.subTest(error=error), patch.object(check.subprocess, 'check_output', side_effect=error):
                self.assertEqual(self.run_guard('--build-dir', str(build)), 1)


if __name__ == '__main__':
    unittest.main()
