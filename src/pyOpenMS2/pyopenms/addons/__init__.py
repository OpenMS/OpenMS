"""
pyOpenMS2 Addon System

This module provides additional Python methods that are injected into
nanobind-wrapped C++ classes at import time.

The addon system allows extending C++ classes with pure Python methods
without modifying the C++ bindings.
"""

from __future__ import annotations

import importlib
import pkgutil
from pathlib import Path
from typing import Any, Callable, Dict, Type

# Registry of addon functions: {class_name: {method_name: function}}
_addon_registry: Dict[str, Dict[str, Callable]] = {}


def addon(class_name: str, method_name: str | None = None):
    """
    Decorator to register an addon method for a class.

    Parameters
    ----------
    class_name : str
        Name of the class to add the method to.
    method_name : str, optional
        Name of the method. If None, uses the function name.

    Examples
    --------
    @addon("MSSpectrum")
    def __repr__(self):
        return f"MSSpectrum(n_peaks={len(self)}, rt={self.getRT():.2f})"

    @addon("MSSpectrum", "custom_method")
    def my_custom_method(self, arg):
        return self.getRT() + arg
    """

    def decorator(func: Callable) -> Callable:
        name = method_name if method_name is not None else func.__name__
        if class_name not in _addon_registry:
            _addon_registry[class_name] = {}
        _addon_registry[class_name][name] = func
        return func

    return decorator


def apply_addons(namespace: Dict[str, Any]) -> None:
    """
    Apply all registered addons to classes in the given namespace.

    Parameters
    ----------
    namespace : dict
        The namespace containing the classes (typically globals() from __init__.py).
    """
    # First, import all addon modules to populate the registry
    _import_addon_modules()

    # Then apply addons to classes
    for class_name, methods in _addon_registry.items():
        if class_name in namespace:
            cls = namespace[class_name]
            for method_name, method in methods.items():
                setattr(cls, method_name, method)


def _import_addon_modules() -> None:
    """Import all addon modules in this package."""
    package_path = Path(__file__).parent

    for module_info in pkgutil.iter_modules([str(package_path)]):
        if module_info.name.startswith("_"):
            continue
        try:
            importlib.import_module(f".{module_info.name}", __package__)
        except ImportError as e:
            import warnings

            warnings.warn(f"Failed to import addon module {module_info.name}: {e}")


# Import standard addon modules
# These provide common functionality like __repr__, __iter__, etc.
# for core classes.
