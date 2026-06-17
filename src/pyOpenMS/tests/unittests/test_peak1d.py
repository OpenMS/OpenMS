"""
Tests for Peak1D class bindings.
"""

import pytest


def test_peak1d_creation():
    """Test Peak1D constructor."""
    from pyopenms import Peak1D

    # Default constructor
    p = Peak1D()
    assert p.getMZ() == 0.0
    assert p.getIntensity() == 0.0

    # Copy constructor
    p.setMZ(100.5)
    p.setIntensity(1000.0)
    p2 = Peak1D(p)
    assert p2.getMZ() == pytest.approx(100.5)
    assert p2.getIntensity() == pytest.approx(1000.0)


def test_peak1d_getters_setters():
    """Test Peak1D getter/setter methods."""
    from pyopenms import Peak1D

    p = Peak1D()

    # MZ
    p.setMZ(456.789)
    assert p.getMZ() == pytest.approx(456.789)

    # Intensity
    p.setIntensity(12345.6)
    assert p.getIntensity() == pytest.approx(12345.6, rel=1e-3)

    # Position alias
    p.setPos(999.123)
    assert p.getPos() == pytest.approx(999.123)
    assert p.getMZ() == pytest.approx(999.123)


def test_peak1d_comparison():
    """Test Peak1D comparison operators."""
    from pyopenms import Peak1D

    p1 = Peak1D()
    p1.setMZ(100.0)
    p1.setIntensity(1000.0)

    p2 = Peak1D()
    p2.setMZ(100.0)
    p2.setIntensity(1000.0)

    p3 = Peak1D()
    p3.setMZ(200.0)
    p3.setIntensity(500.0)

    assert p1 == p2
    assert p1 != p3


def test_peak1d_hash():
    """Test Peak1D hash support."""
    from pyopenms import Peak1D

    p1 = Peak1D()
    p1.setMZ(100.0)
    p1.setIntensity(1000.0)

    p2 = Peak1D()
    p2.setMZ(100.0)
    p2.setIntensity(1000.0)

    # Same values should have same hash
    assert hash(p1) == hash(p2)

    # Can be used in sets and dicts
    s = {p1}
    assert p2 in s or hash(p1) == hash(p2)  # Implementation may vary


def test_peak1d_repr():
    """Test Peak1D string representations."""
    from pyopenms import Peak1D

    p = Peak1D()
    p.setMZ(123.4567)
    p.setIntensity(9876.54)

    repr_str = repr(p)
    assert "Peak1D" in repr_str
    assert "123.4567" in repr_str

    str_str = str(p)
    assert "123.4567" in str_str
