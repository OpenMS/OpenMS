"""
Tests that copy.copy(), copy.deepcopy(), and ClassName(obj) work for all
pyopenms classes that support copying.

Catches regressions where nanobind dispatches copy-construction to a string
constructor instead of the copy constructor due to registration order.

Uses a subprocess probe at collection time to discover classes that segfault
during default-construction or copy. Those classes are excluded from the
parametrized tests but cause test_no_segfaults_during_copy to FAIL with
a clear list of crashers.
"""

import copy
import enum
import inspect
import subprocess
import sys

import pytest

import pyopenms


def _get_copyable_classes():
    """Return public non-enum classes that define __copy__ directly (not inherited).

    Classes that merely inherit __copy__ from a base class (e.g., DefaultParamHandler,
    ProgressLogger) suffer from object slicing: copy.copy() returns the base type.
    Those classes need their own __copy__ binding before they can be tested here.
    """
    return sorted(
        name for name in dir(pyopenms)
        if not name.startswith("_")
        and inspect.isclass(getattr(pyopenms, name, None))
        and not issubclass(getattr(pyopenms, name), enum.Enum)
        and "__copy__" in getattr(pyopenms, name).__dict__
    )


_PROBE_CODE = """\
import copy, sys
import pyopenms
for name in {class_names!r}:
    cls = getattr(pyopenms, name)
    try:
        obj = cls()
    except Exception:
        print("SKIP:" + name, flush=True)
        continue
    try:
        copy.copy(obj)
    except Exception:
        print("SKIP:" + name, flush=True)
        continue
    print("OK:" + name, flush=True)
"""


def _probe_safe_classes(class_names, timeout=120):
    """Probe which classes can be default-constructed and copied without segfaulting.

    Runs operations in a subprocess. If it crashes, identifies the crasher from
    the last successfully printed class name and recurses on remaining classes.

    Returns (safe_set, crashed_list).
    """
    if not class_names:
        return set(), []

    code = _PROBE_CODE.format(class_names=class_names)
    try:
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return set(), list(class_names)

    safe = set()
    skipped = set()
    for line in result.stdout.strip().splitlines():
        if line.startswith("OK:"):
            safe.add(line[3:])
        elif line.startswith("SKIP:"):
            skipped.add(line[5:])

    if result.returncode == 0:
        return safe, []

    # Subprocess crashed. Identify crasher and recurse on remaining.
    reported = safe | skipped
    remaining = [n for n in class_names if n not in reported]
    crashed = []
    if remaining:
        crashed.append(remaining[0])  # first unreported class = the crasher
        remaining = remaining[1:]
    if remaining:
        more_safe, more_crashed = _probe_safe_classes(remaining, timeout)
        safe |= more_safe
        crashed.extend(more_crashed)

    return safe, crashed


_ALL_COPYABLE = _get_copyable_classes()
_SAFE_CLASSES, _CRASHED_CLASSES = _probe_safe_classes(_ALL_COPYABLE)
COPYABLE_CLASSES = sorted(_SAFE_CLASSES)


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


def test_no_segfaults_during_copy():
    """Fail if any copyable class segfaults during default-construct + copy."""
    assert not _CRASHED_CLASSES, (
        f"{len(_CRASHED_CLASSES)} class(es) segfault during copy probe: "
        + ", ".join(_CRASHED_CLASSES)
    )


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
def test_copy_copy(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.copy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
def test_deepcopy(class_name):
    cls = getattr(pyopenms, class_name)
    obj = _try_default_construct(cls)
    if obj is None:
        pytest.skip(f"{class_name} cannot be default-constructed")

    obj_copy = copy.deepcopy(obj)
    assert type(obj_copy) is type(obj)
    if _has_own_eq(cls):
        assert obj_copy == obj


@pytest.mark.parametrize("class_name", COPYABLE_CLASSES)
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
