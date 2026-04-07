"""
Tests that copy.copy(), copy.deepcopy(), and ClassName(obj) work for all
pyopenms classes that support copying.

Catches regressions where nanobind dispatches copy-construction to a string
constructor instead of the copy constructor due to registration order.
"""

import copy
import enum
import inspect

import pytest

import pyopenms


def _get_all_classes():
    """Return all public non-enum classes exported by pyopenms."""
    return sorted(
        name for name in dir(pyopenms)
        if not name.startswith("_")
        and inspect.isclass(getattr(pyopenms, name, None))
        and not issubclass(getattr(pyopenms, name), enum.Enum)
    )


ALL_CLASSES = _get_all_classes()


def _try_default_construct(cls):
    """Try to default-construct a class instance. Returns instance or None."""
    try:
        return cls()
    except (TypeError, RuntimeError):
        return None


def _has_own_eq(cls):
    """Check if cls defines __eq__ (not inherited from object)."""
    return any(
        "__eq__" in base.__dict__
        for base in cls.__mro__
        if base is not object
    )


@pytest.mark.parametrize("class_name", ALL_CLASSES)
def test_copy_copy(class_name):
    cls = getattr(pyopenms, class_name)
    if not hasattr(cls, "__copy__"):
        pytest.skip(f"{class_name} does not define __copy__")
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.copy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", ALL_CLASSES)
def test_deepcopy(class_name):
    cls = getattr(pyopenms, class_name)
    if not hasattr(cls, "__deepcopy__"):
        pytest.skip(f"{class_name} does not define __deepcopy__")
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.deepcopy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", ALL_CLASSES)
def test_copy_constructor(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    try:
        obj_copy = cls(obj)
    except TypeError:
        pytest.skip(f"{class_name} does not support copy construction")

    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj
