"""Regression tests for generated stub syntax and memory-query fallbacks.

$Maintainer: Timo Sachsenberg $
"""

import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

import pytest


# The stub post-processor is build tooling, not part of the package: it sits
# next to tests/ in the source tree, the in-tree build copies it next to its
# copy of tests/, and the wheel test job checks it out alongside tests/.
FIX_STUBS = Path(__file__).resolve().parents[1] / "fix_stubs.py"


def load_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def package_dir():
    """Directory of the pyopenms package under test, located without importing it.

    The wheel test job runs this suite from a sparse checkout that holds only
    tests/ next to an installed wheel, so package files are resolved through the
    import system -- wherever the wheel or the in-tree build put them -- rather
    than through a source tree that may not be checked out.
    """
    spec = importlib.util.find_spec("pyopenms")
    assert spec is not None and spec.submodule_search_locations, "pyopenms is not importable"
    return Path(next(iter(spec.submodule_search_locations)))


@pytest.mark.parametrize("platform", ["win32", "unsupported"])
def test_free_mem_fallback_is_a_named_function(platform):
    # Load only the helper, without importing native bindings. Force the
    # Windows fallback even on a machine that happens to have pywin32 installed.
    helper = package_dir() / "_sysinfo.py"
    with patch.object(sys, "platform", platform), patch.dict(sys.modules, {"win32api": None}):
        sysinfo = load_source("_sysinfo_test", helper)
    assert sysinfo.free_mem() == 0
    assert sysinfo.free_mem.__name__ == "free_mem"
    assert sysinfo.free_mem.__annotations__["return"] is int


@pytest.mark.parametrize("content", [
    "free_mem = <lambda>\n",
    "from pyopenms._sysinfo import <lambda> as free_mem\n",
])
def test_postprocessing_rejects_invalid_stub_syntax(tmp_path, content):
    fixes = load_source("_fix_stubs_test", FIX_STUBS)
    stub = tmp_path / "broken.pyi"
    stub.write_text(content, encoding="utf-8")
    with pytest.raises(SyntaxError) as error:
        fixes.fix_stub_file(stub)
    assert error.value.filename == str(stub)


def test_postprocessing_validates_after_repairs(tmp_path):
    fixes = load_source("_fix_stubs_test", FIX_STUBS)
    stub = tmp_path / "repaired.pyi"
    stub.write_text('"""測定"""\ndef example(from: int) -> None_: ...\n', encoding="utf-8")
    assert fixes.fix_stub_file(stub)
    tree = ast.parse(stub.read_text(encoding="utf-8"))
    function = tree.body[1]
    assert function.args.args[0].arg == "from_"
    assert isinstance(function.returns, ast.Constant) and function.returns.value is None


def test_installed_stubs_are_valid_python():
    import pyopenms

    package = Path(pyopenms.__file__).resolve().parent
    if not (package / "py.typed").is_file():
        pytest.skip("This build has stub generation disabled")
    stubs = sorted(package.rglob("*.pyi"))
    assert stubs, "py.typed must not advertise a package without generated stubs"
    for stub in stubs:
        ast.parse(stub.read_text(encoding="utf-8"), filename=str(stub))


def exported_names(tree):
    """Names a stub module exports to a type checker."""
    names = set()
    for node in tree.body:
        if isinstance(node, (ast.ClassDef, ast.FunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(target.id for target in node.targets if isinstance(target, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.ImportFrom):
            # A stub re-exports an imported name only in the form "import X as X".
            names.update(alias.name for alias in node.names if alias.asname == alias.name)
    return names


def fresh_namespace():
    """{public name: module name, or None} of pyopenms right after "import pyopenms".

    Checked in a new interpreter, as stubgen imports the package: within the test
    session other tests import submodules such as pyopenms.plotting, which the import
    system then binds as attributes of the package.
    """
    code = (
        "import inspect, json, pyopenms\n"
        "names = {name: getattr(pyopenms, name).__name__ if inspect.ismodule(getattr(pyopenms, name)) else None\n"
        "         for name in dir(pyopenms) if not name.startswith('_')}\n"
        "print('NAMESPACE ' + json.dumps(names))\n")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, f"import pyopenms failed:\n{result.stderr}\n{result.stdout}"
    line = [line for line in result.stdout.splitlines() if line.startswith("NAMESPACE ")][-1]
    return json.loads(line[len("NAMESPACE "):])


def test_installed_stubs_export_the_package_namespace():
    package = package_dir()
    if not (package / "py.typed").is_file():
        pytest.skip("This build has stub generation disabled")
    init = ast.parse((package / "__init__.pyi").read_text(encoding="utf-8"))

    # Every public name is declared and exported. An alias such as PeakMap was
    # imported under another name, which is no export, and a stray import such as
    # ctypes was not declared at all.
    public = set(fresh_namespace())
    assert sorted(public - exported_names(init)) == []

    # Every name imported from a module of the package exists there. stubgen used
    # to import a nested enum such as FileTypes.FileType as a top-level name.
    unresolved = []
    for node in init.body:
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("pyopenms."):
            source = package.joinpath(*node.module.split(".")[1:])
            stub = source / "__init__.pyi" if source.is_dir() else source.with_suffix(".pyi")
            names = exported_names(ast.parse(stub.read_text(encoding="utf-8")))
            unresolved += [f"{node.module}.{alias.name}" for alias in node.names
                           if alias.name not in names]
    assert unresolved == []


def test_package_namespace_holds_no_stray_imports():
    namespace = fresh_namespace()
    stray = [name for name, module in namespace.items()
             if module is not None and not module.startswith("pyopenms.")]
    # "from __future__ import annotations" binds a name as well.
    if "annotations" in namespace:
        stray.append("annotations")
    assert stray == []
