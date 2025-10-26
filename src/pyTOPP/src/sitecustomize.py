# Auto-imported by Python's 'site' module if this directory is on sys.path.
# Extends/creates the 'pyopenms' namespace to include the vendored 'pytopp'.

import importlib
import importlib.machinery
import os
import sys
import types

def _ensure_pyopenms_namespace():
    site = os.path.dirname(os.path.abspath(__file__))
    pyo_dir = os.path.join(site, "pyopenms")

    try:
        m = importlib.import_module("pyopenms")
    except Exception:
        m = None

    if m is not None:
        try:
            path_list = list(m.__path__)  # may raise if not a pkg
        except Exception:
            path_list = []
        if pyo_dir not in path_list:
            try:
                m.__path__.append(pyo_dir)
            except Exception:
                m.__path__ = path_list + [pyo_dir]
        return

    m = types.ModuleType("pyopenms")
    m.__package__ = "pyopenms"
    m.__path__ = [pyo_dir]
    m.__spec__ = importlib.machinery.ModuleSpec("pyopenms", loader=None, is_package=True)
    sys.modules["pyopenms"] = m

_ensure_pyopenms_namespace()
