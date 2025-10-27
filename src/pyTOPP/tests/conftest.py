import os
import sys
from pathlib import Path
import importlib
import importlib.machinery
import types
import pytest

# ── Put source tree on sys.path so we import the in-tree pyTOPP during tests
#    and make sure the 'pyopenms' namespace points at it.
SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Ensure 'pyopenms' is importable as a namespace that includes SRC_DIR/pyopenms
try:
    m = importlib.import_module("pyopenms")
except Exception:
    m = types.ModuleType("pyopenms")
    m.__package__ = "pyopenms"
    m.__path__ = [str(SRC_DIR / "pyopenms")]
    m.__spec__ = importlib.machinery.ModuleSpec("pyopenms", loader=None, is_package=True)
    sys.modules["pyopenms"] = m
else:
    # Extend namespace with our in-tree pyopenms if not already present
    try:
        pth = list(m.__path__)
    except Exception:
        pth = []
    want = str(SRC_DIR / "pyopenms")
    if want not in pth:
        try:
            m.__path__.append(want)
        except Exception:
            m.__path__ = pth + [want]


# ── Tiny test doubles to stand in for CTDopts Parameter/Model where needed
class _P:
    """Minimal parameter stub used by tests."""
    def __init__(self, name, typ, *, required=False, description="", choices=None, restrictions=None, file_formats=None, is_list=False, tags=None, default=None):
        self.name = name
        self.type = typ
        self.required = required
        self.description = description
        self.choices = choices
        self.restrictions = restrictions
        self.file_formats = file_formats
        self.is_list = is_list
        self.tags = tags or []
        self.default = default

class _Model:
    """Minimal CTD model stub with the APIs used by ctdsupport helpers."""
    def __init__(self, params, defaults=None, *, name="Tool", desc="", docurl=""):
        self._params = list(params)
        self._defaults = self._as_nested(defaults or {})
        self.name = name
        self.opt_attribs = {"description": desc, "docurl": docurl}

    def list_parameters(self):
        return list(self._params)

    def get_defaults(self):
        return self._defaults

    @staticmethod
    def _as_nested(flat):
        """Turn {'a:b:c': 1, 'x': 2} into {'a': {'b': {'c': 1}}, 'x': 2}."""
        nested = {}
        for k, v in (flat or {}).items():
            cur = nested
            if isinstance(k, str) and ":" in k:
                parts = k.split(":")
                for p in parts[:-1]:
                    cur = cur.setdefault(p, {})
                cur[parts[-1]] = v
            else:
                cur[k] = v
        return nested


@pytest.fixture
def param_stub():
    return _P

@pytest.fixture
def model_stub():
    return _Model
