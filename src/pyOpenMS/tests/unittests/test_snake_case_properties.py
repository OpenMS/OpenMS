"""
Tests for the Pythonic snake_case properties on core classes (issue #9760).

These properties are additive: the existing getX()/setX() methods are unchanged.
Only clean, no-arg scalar get/set pairs are exposed as properties.
"""

import pytest


def test_peak1d_scalar_properties_roundtrip():
    """Peak1D.mz / .intensity read and write, and stay in sync with the getters."""
    from pyopenms import Peak1D

    pk = Peak1D()
    pk.mz = 250.5
    pk.intensity = 900.0

    # property read
    assert pk.mz == pytest.approx(250.5)
    assert pk.intensity == pytest.approx(900.0)
    # property write is reflected in the legacy getters (additive, same storage)
    assert pk.getMZ() == pytest.approx(250.5)
    assert pk.getIntensity() == pytest.approx(900.0)
    # legacy setter is reflected in the property
    pk.setMZ(111.0)
    assert pk.mz == pytest.approx(111.0)


def test_peak1d_pos_alias_is_not_a_property():
    """'pos' is an alias of mz; it must NOT be exposed as a property (avoid alias pairs)."""
    from pyopenms import Peak1D

    assert not hasattr(Peak1D, "pos")
    # but the legacy method alias still works
    pk = Peak1D()
    pk.mz = 42.0
    assert pk.getPos() == pytest.approx(42.0)


def test_msspectrum_scalar_properties_roundtrip():
    """MSSpectrum scalar properties read/write and stay in sync with getters."""
    from pyopenms import MSSpectrum

    s = MSSpectrum()
    s.rt = 205.2
    s.ms_level = 2
    s.drift_time = 25.5
    s.name = "scan=1"
    s.comment = "hello"
    s.native_id = "controllerType=0 scan=1"

    assert s.rt == pytest.approx(205.2)
    assert s.ms_level == 2
    assert s.drift_time == pytest.approx(25.5)
    assert s.name == "scan=1"
    assert s.comment == "hello"
    assert s.native_id == "controllerType=0 scan=1"

    # agreement with the legacy getters
    assert s.getRT() == pytest.approx(205.2)
    assert s.getMSLevel() == 2
    assert s.getName() == "scan=1"
    assert s.getNativeID() == "controllerType=0 scan=1"

    # legacy setter reflected in the property
    s.setRT(9.0)
    assert s.rt == pytest.approx(9.0)


def test_msspectrum_drift_time_unit_enum_property():
    """drift_time_unit is an enum-valued property backed by get/setDriftTimeUnit."""
    from pyopenms import MSSpectrum, DriftTimeUnit

    s = MSSpectrum()
    s.drift_time_unit = DriftTimeUnit.MILLISECOND
    assert s.drift_time_unit == DriftTimeUnit.MILLISECOND
    assert s.getDriftTimeUnit() == DriftTimeUnit.MILLISECOND


def test_getitem_returns_reference_inplace_mutation():
    """spectrum[i] returns a reference: spectrum[i].mz = x mutates the real peak (#9760)."""
    from pyopenms import MSSpectrum

    s = MSSpectrum()
    s.set_peaks(([100.0, 200.0], [10.0, 20.0]))

    # mutate through the indexed peak
    s[0].mz = 999.0
    s[1].intensity = 55.0

    # the change is visible on a fresh index read ...
    assert s[0].mz == pytest.approx(999.0)
    assert s[1].intensity == pytest.approx(55.0)
    # ... and in the underlying peak data
    mz, intensity = s.get_peaks()
    assert mz[0] == pytest.approx(999.0)
    assert intensity[1] == pytest.approx(55.0)


def test_getitem_reference_invalidation_is_documented_behavior():
    """
    A reference obtained from spectrum[i] dangles if the spectrum reallocates
    (push_back / resize / set_peaks), mirroring C++ iterator invalidation.

    We do NOT dereference a stale reference here (that would be undefined behavior).
    Instead we assert the safe pattern: re-index after mutation always sees the
    current data.
    """
    from pyopenms import MSSpectrum, Peak1D

    s = MSSpectrum()
    s.set_peaks(([100.0], [10.0]))

    # Reallocating the peak container; do not reuse any previously held reference.
    s.push_back(Peak1D())
    s[1].mz = 300.0

    # Safe: fresh indexing reflects the new state.
    assert s.size() == 2
    assert s[1].mz == pytest.approx(300.0)


def test_getters_setters_still_present_additive():
    """The legacy methods must remain (properties are purely additive)."""
    from pyopenms import Peak1D, MSSpectrum

    for name in ("getMZ", "setMZ", "getIntensity", "setIntensity", "getPos", "setPos"):
        assert hasattr(Peak1D, name)
    for name in ("getRT", "setRT", "getMSLevel", "setMSLevel", "getName", "setName",
                 "getNativeID", "setNativeID", "getComment", "setComment"):
        assert hasattr(MSSpectrum, name)
