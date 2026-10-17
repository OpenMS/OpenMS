"""Tests for the ion mobility *unit* accessors on MSSpectrum.

OpenMS stores ion mobility two different ways, and a spectrum may carry either -- or, on some
vendor readers, both:

* **scalar** -- the whole spectrum was acquired at a single drift time, in ``drift_time_`` /
  ``drift_time_unit_``, reached through ``getDriftTime()`` / ``getDriftTimeUnit()`` and the
  ``drift_time`` / ``drift_time_unit`` properties;
* **per-peak** -- one mobility per peak, in a ``FloatDataArray`` whose *name* is a PSI-MS CV term,
  reached through ``containsIMData()`` / ``getIMData()`` / ``get_drift_time_array()``.

``get_drift_time_unit()`` guarded on the per-peak array but returned the scalar, so it mixed the
two.  It answered ``DriftTimeUnit.NONE`` for a frame whose array was perfectly well labelled, and
``None`` -- "no ion mobility here" -- for a spectrum whose scalar unit was set, which is the very
field it returned.  It is replaced by ``get_drift_time_array_unit()``, which is array-derived
throughout; the old spelling remains for one release as a deprecated alias with the corrected
behaviour.

``test_unit_is_stable_across_an_mzml_roundtrip`` is the one that pins down why this mattered in
practice rather than only on paper.  ``BrukerTimsFile::frameToSpectrum`` sets the scalar *unit* to
VSSC and attaches a VSSC array but never sets a scalar drift *time*, and mzML only writes the
scalar unit when the time is set.  The old accessor therefore reported VSSC for a TIMS frame in
memory and ``NONE`` for the same frame after a save/load, with the array byte-identical
throughout.
"""

import warnings

import numpy as np
import pytest

import pyopenms


# CV term name -> the unit it encodes
_VSSC_ARRAY = "raw inverse reduced ion mobility array"
_MS_ARRAY = "mean ion mobility array"


def _peaks(spec, n=2):
    spec.set_peaks(np.arange(1.0, n + 1.0), (np.arange(1.0, n + 1.0) * 10.0).astype(np.float32))
    return spec


def _with_array(name, n=2):
    """A frame: per-peak ion mobility array, scalars left unset."""
    spec = _peaks(pyopenms.MSSpectrum(), n)
    fda = pyopenms.FloatDataArray()
    fda.setName(name)
    fda.set_data((np.arange(1.0, n + 1.0) / 10.0).astype(np.float32))
    spec.setFloatDataArrays([fda])
    return spec


def _with_scalar(unit):
    """An IM_SPECTRUM: one drift time for the whole spectrum, no array."""
    spec = _peaks(pyopenms.MSSpectrum())
    spec.setDriftTime(0.5)
    spec.setDriftTimeUnit(unit)
    return spec


def _with_both(array_name, scalar_unit):
    spec = _with_array(array_name)
    spec.setDriftTimeUnit(scalar_unit)
    return spec


# The five states the accessor has to distinguish, and what it must report for each.
_STATES = [
    pytest.param(lambda: _with_array(_VSSC_ARRAY), pyopenms.DriftTimeUnit.VSSC, id="array_only"),
    pytest.param(lambda: _with_scalar(pyopenms.DriftTimeUnit.MILLISECOND), None, id="scalar_only"),
    pytest.param(lambda: _with_both(_VSSC_ARRAY, pyopenms.DriftTimeUnit.VSSC),
                 pyopenms.DriftTimeUnit.VSSC, id="both_agreeing"),
    pytest.param(lambda: _with_both(_VSSC_ARRAY, pyopenms.DriftTimeUnit.MILLISECOND),
                 pyopenms.DriftTimeUnit.VSSC, id="both_disagreeing_array_wins"),
    pytest.param(lambda: _peaks(pyopenms.MSSpectrum()), None, id="neither"),
]


@pytest.mark.parametrize("build,expected", _STATES)
def test_array_unit_is_array_derived_in_every_state(build, expected):
    """None means "no ion mobility array"; it never means "the scalar was unset"."""
    assert build().get_drift_time_array_unit() == expected


@pytest.mark.parametrize("build,expected", _STATES)
def test_the_scalar_family_is_unaffected(build, expected):
    """The fix must not disturb the accessors that legitimately report the scalar."""
    spec = build()
    # getDriftTimeUnit() and the property are the scalar, whatever the array says
    assert spec.getDriftTimeUnit() == spec.drift_time_unit


def test_array_unit_and_scalar_unit_can_disagree_and_both_are_reported_honestly():
    spec = _with_both(_VSSC_ARRAY, pyopenms.DriftTimeUnit.MILLISECOND)
    assert spec.get_drift_time_array_unit() == pyopenms.DriftTimeUnit.VSSC
    assert spec.getDriftTimeUnit() == pyopenms.DriftTimeUnit.MILLISECOND
    assert spec.getDriftTimeUnitAsString() == "ms"


@pytest.mark.parametrize("name,unit", [(_VSSC_ARRAY, pyopenms.DriftTimeUnit.VSSC),
                                       (_MS_ARRAY, pyopenms.DriftTimeUnit.MILLISECOND),
                                       ("Ion Mobility", pyopenms.DriftTimeUnit.MILLISECOND)])
def test_getIMData_returns_the_unit_as_an_enum(name, unit):
    """It used to downcast to a bare int, which left the array's unit the only one not
    reachable as an enum -- and is why callers reached for the scalar accessor instead."""
    spec = _with_array(name)
    index, found = spec.getIMData()

    assert isinstance(found, pyopenms.DriftTimeUnit)
    assert found == unit
    assert spec.getFloatDataArrays()[index].getName() == name


def test_getIMData_unit_still_compares_and_converts_like_an_int():
    """DriftTimeUnit is bound with nb::is_arithmetic(), so the enum change is not a break
    for equality, int(), isinstance or dict lookup -- only for str()/type() identity."""
    _, unit = _with_array(_VSSC_ARRAY).getIMData()
    assert unit == int(pyopenms.DriftTimeUnit.VSSC)
    assert int(unit) == 2
    assert isinstance(unit, int)
    assert {2: "ok"}[unit] == "ok"


def test_deprecated_alias_delegates_and_warns():
    spec = _with_array(_VSSC_ARRAY)
    with pytest.deprecated_call():
        value = spec.get_drift_time_unit()
    assert value == spec.get_drift_time_array_unit() == pyopenms.DriftTimeUnit.VSSC


def test_deprecated_alias_reports_the_corrected_value():
    """The alias is not merely renamed: the old implementation returned NONE here."""
    spec = _with_array(_VSSC_ARRAY)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert spec.get_drift_time_unit() == pyopenms.DriftTimeUnit.VSSC


@pytest.mark.parametrize("name,expected", [(_VSSC_ARRAY, "1/K0"), (_MS_ARRAY, "ms")])
def test_dataframe_unit_column_comes_from_the_array(name, expected):
    """to_df() took its values from the array but its unit string from the scalar, so a frame
    -- which leaves the scalar unset -- reported '<NONE>' for data that had a unit."""
    spec = _with_array(name)
    df = spec.to_df(columns=["mz", "intensity", "ion_mobility", "ion_mobility_unit"])

    assert list(df["ion_mobility_unit"]) == [expected] * len(df)
    assert np.allclose(df["ion_mobility"], [0.1, 0.2])


def test_dataframe_unit_column_follows_the_array_when_the_two_disagree():
    spec = _with_both(_VSSC_ARRAY, pyopenms.DriftTimeUnit.MILLISECOND)
    df = spec.to_df(columns=["mz", "ion_mobility_unit"])
    assert list(df["ion_mobility_unit"]) == ["1/K0"] * len(df)


def test_unit_is_stable_across_an_mzml_roundtrip(tmp_path):
    """The TIMS shape: scalar unit set, scalar time unset, per-peak array present.

    mzML only writes the scalar unit when the scalar time is set, so the scalar is dropped by
    the round-trip while the array survives. An array-derived accessor is therefore the only
    one that answers the same before and after.
    """
    spec = _with_both(_VSSC_ARRAY, pyopenms.DriftTimeUnit.VSSC)
    spec.setRT(1.0)
    spec.setMSLevel(1)
    assert spec.getDriftTime() == -1.0  # scalar time deliberately unset, as BrukerTimsFile leaves it

    exp = pyopenms.MSExperiment()
    exp.addSpectrum(spec)
    path = str(tmp_path / "tims.mzML")
    pyopenms.MzMLFile().store(path, exp)
    loaded = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(path, loaded)
    back = loaded[0]

    assert np.allclose(back.get_drift_time_array(), spec.get_drift_time_array())
    assert back.get_drift_time_array_unit() == pyopenms.DriftTimeUnit.VSSC
    assert back.get_drift_time_array_unit() == spec.get_drift_time_array_unit()
