"""Dependency-free regression checks: python -m unittest discover -s tests -p test_library_packaging.py."""
import ast
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch


SOURCE = Path(__file__).resolve().parents[1]


class LibraryPackagingTest(unittest.TestCase):
    def test_cleanup_retains_current_copies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            current = ["libOpenMS.so.3.7.0", "libOpenSwathAlgo.so.3.7.0"]
            stale = ["libOpenMS.so", "libOpenSwathAlgo.so", "libOpenMS.so.3.6.0",
                     "libOpenSwathAlgo.so.3.6.0", "libOpenMS.dylib", "libOpenSwathAlgo.dylib"]
            for name in current + stale:
                (root / name).write_text(name)
            timestamps = [(root / name).stat().st_mtime_ns for name in current]
            for _ in range(2):
                subprocess.run(["cmake", f"-DDIR={root}", f"-DCURRENT_OPENMS={current[0]}",
                                f"-DCURRENT_SWATH={current[1]}", "-P",
                                str(SOURCE / "pyopenms_clean_old_libs.cmake")], check=True)
                self.assertEqual(sorted(p.name for p in root.iterdir()), sorted(current))
                self.assertEqual(timestamps, [(root / name).stat().st_mtime_ns for name in current])

    def test_preload_uses_exact_pair(self):
        tree = ast.parse((SOURCE / "pyopenms/__init__.py").read_text())
        block = next(node for node in tree.body if isinstance(node, ast.If)
                     and ast.unparse(node.test) == "sys.platform.startswith('linux')")
        code = compile(ast.Module(body=[block], type_ignores=[]), "preload", "exec")
        names = types.ModuleType("pyopenms._library_names")
        names.OPENMS_LIBRARY = "libOpenMS.so.3.7.0"
        names.OPENSWATH_LIBRARY = "libOpenSwathAlgo.so.3.7.0"
        with tempfile.TemporaryDirectory(prefix="openms[paths]") as directory:
            root = Path(directory)
            for name in [names.OPENMS_LIBRARY, names.OPENSWATH_LIBRARY,
                         "libOpenMS.so.3.6.0", "libOpenSwathAlgo.so"]:
                (root / name).touch()
            scope = dict(os=os, sys=types.SimpleNamespace(platform="linux"), here=directory,
                         __package__="pyopenms")
            with patch.dict(sys.modules, {"pyopenms._library_names": names}), \
                    patch("ctypes.cdll.LoadLibrary") as load:
                exec(code, scope.copy())
                self.assertEqual([call.args[0] for call in load.call_args_list],
                                 [str(root / names.OPENSWATH_LIBRARY), str(root / names.OPENMS_LIBRARY)])
                (root / names.OPENSWATH_LIBRARY).unlink()
                load.reset_mock()
                with self.assertRaisesRegex(ImportError, "Missing matching"):
                    exec(code, scope.copy())
                load.assert_not_called()
                (root / names.OPENMS_LIBRARY).unlink()
                exec(code, scope.copy())
                load.assert_not_called()
