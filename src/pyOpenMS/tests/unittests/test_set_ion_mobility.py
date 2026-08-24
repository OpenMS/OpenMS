"""Tests for writing the per-peak ion mobility dimension from Python (issue #9792).

Before this, pyOpenMS could *read* the ion mobility array five ways -- ``containsIMData()``,
``getIMData()``, ``get_drift_time_array()``, ``drift_time_array_view()``, ``rasterizeIMFrame()``
-- and could not write it at all.  The only route was the generic data-array API with the name
spelled by hand::

    fda = FloatDataArray()
    fda.set_data(im)
    fda.setName("Ion Mobility")     # this string is the entire contract
    spec.setFloatDataArrays([fda])

which is easy to get wrong in a way that silently changes the *unit*: a name starting with
"Ion Mobility" and carrying no CV accession is read back as milliseconds, so labelling TIMS 1/K0
data that way misreports it.  ``set_drift_time_array(values, unit)`` and
``set_peaks(..., ion_mobility=...)`` now name the array through ``IMDataConverter::setIMUnit()``,
i.e. from the PSI-MS controlled vocabulary, so what is written is what ``getIMData()`` reads back
and what mzML round-trips (``test_mzml_roundtrip_preserves_ion_mobility``).

``DriftTimeUnit.MILLISECOND`` and ``DriftTimeUnit.VSSC`` are the only units the PSI-MS CV has an
ion-mobility-array term for, so the core rejects the rest; that is inherited here rather than
worked around, see ``test_units_without_a_cv_term_are_rejected``.
"""

import numpy as np
import pytest

import pyopenms


# unit -> the CV term name IMDataConverter writes for it
_SUPPORTED_UNITS = [
    pytest.param(pyopenms.DriftTimeUnit.MILLISECOND, "mean ion mobility array", id="MILLISECOND"),
    pytest.param(pyopenms.DriftTimeUnit.VSSC, "raw inverse reduced ion mobility array", id="VSSC"),
]

_CALL_FORMS = [
    pytest.param(lambda s, p, i, **kw: s.set_peaks(p, i, **kw), id="two_positional"),
    pytest.param(lambda s, p, i, **kw: s.set_peaks((p, i), **kw), id="tuple"),
    pytest.param(lambda s, p, i, **kw: s.set_peaks([p, i], **kw), id="list"),
]


def _peaks(n):
    return np.arange(1.0, n + 1.0), (np.arange(1.0, n + 1.0) * 10.0).astype(np.float32)


def _im(n):
    return (np.arange(1.0, n + 1.0) / 10.0).astype(np.float32)


def _spectrum(n, im_unit=None):
    s = pyopenms.MSSpectrum()
    s.set_peaks(*_peaks(n))
    if im_unit is not None:
        s.set_drift_time_array(_im(n), im_unit)
    return s


@pytest.mark.parametrize("unit,cv_name", _SUPPORTED_UNITS)
@pytest.mark.parametrize("call", _CALL_FORMS)
def test_set_peaks_installs_a_discoverable_ion_mobility_array(unit, cv_name, call):
    """The array written by set_peaks() is the one containsIMData()/getIMData() find."""
    spec = pyopenms.MSSpectrum()
    pos, inten = _peaks(3)
    call(spec, pos, inten, ion_mobility=_im(3), ion_mobility_unit=unit)

    assert spec.containsIMData()
    index, found_unit = spec.getIMData()
    assert found_unit == int(unit)
    assert spec.getFloatDataArrays()[index].getName() == cv_name
    assert np.allclose(spec.get_drift_time_array(), _im(3))


@pytest.mark.parametrize("unit,cv_name", _SUPPORTED_UNITS)
def test_set_drift_time_array_standalone(unit, cv_name):
    spec = _spectrum(4)
    spec.set_drift_time_array(_im(4), unit)

    assert spec.containsIMData()
    assert spec.getFloatDataArrays()[spec.getIMData()[0]].getName() == cv_name
    assert np.allclose(spec.get_drift_time_array(), _im(4))


@pytest.mark.parametrize("unit,_cv", _SUPPORTED_UNITS)
def test_setting_twice_replaces_rather_than_duplicates(unit, _cv):
    """A second write must not leave two ion mobility arrays, which getIMData() cannot disambiguate."""
    spec = _spectrum(4, unit)
    spec.set_drift_time_array(_im(4) * 2.0, unit)

    assert len(spec.getFloatDataArrays()) == 1
    assert np.allclose(spec.get_drift_time_array(), _im(4) * 2.0)


@pytest.mark.parametrize("unit,_cv", _SUPPORTED_UNITS)
def test_unit_is_inherited_from_an_existing_array(unit, _cv):
    """Editing the mobilities of a spectrum that already has them needs no unit argument."""
    spec = _spectrum(4, unit)
    spec.set_drift_time_array(_im(4) * 3.0)

    assert spec.getIMData()[1] == int(unit)
    assert np.allclose(spec.get_drift_time_array(), _im(4) * 3.0)


def test_unit_is_required_when_there_is_nothing_to_inherit_from():
    spec = _spectrum(3)
    with pytest.raises(ValueError, match="ion_mobility_unit is required"):
        spec.set_drift_time_array(_im(3))
    assert not spec.containsIMData()


@pytest.mark.parametrize("unit", [pyopenms.DriftTimeUnit.CCS,
                                 pyopenms.DriftTimeUnit.FAIMS_COMPENSATION_VOLTAGE])
def test_units_without_a_cv_term_are_rejected(unit):
    """The PSI-MS CV has no ion-mobility-array term for these, so the core refuses to name them."""
    spec = _spectrum(3)
    with pytest.raises(RuntimeError):
        spec.set_drift_time_array(_im(3), unit)
    assert not spec.containsIMData()
    assert len(spec.getFloatDataArrays()) == 0


@pytest.mark.parametrize("n_im", [2, 7])
def test_wrong_length_is_rejected_and_changes_nothing(n_im):
    spec = _spectrum(3)
    with pytest.raises(RuntimeError, match="ion mobility array has"):
        spec.set_drift_time_array(_im(n_im), pyopenms.DriftTimeUnit.VSSC)

    assert spec.size() == 3
    assert len(spec.getFloatDataArrays()) == 0
    assert np.allclose(spec.get_peaks()[0], _peaks(3)[0])


@pytest.mark.parametrize("call", _CALL_FORMS)
def test_set_peaks_rejects_a_mismatched_ion_mobility_array_atomically(call):
    """The ion mobility array is built before anything is modified, so this is a no-op."""
    spec = _spectrum(3)
    with pytest.raises(RuntimeError, match="ion mobility array has"):
        call(spec, *_peaks(5), ion_mobility=_im(4), ion_mobility_unit=pyopenms.DriftTimeUnit.VSSC)

    assert spec.size() == 3
    assert np.allclose(spec.get_peaks()[0], _peaks(3)[0])
    assert len(spec.getFloatDataArrays()) == 0


@pytest.mark.parametrize("call", _CALL_FORMS)
def test_replacing_the_ion_mobility_array_exempts_it_from_the_metadata_guard(call):
    """Changing the peak count *and* the mobilities in one call is the point of the keyword.

    The old array is mis-sized against the new count, but it is being overwritten, so refusing
    would leave no way to do this in one step. The unit is inherited from the array being replaced.
    """
    spec = _spectrum(100, pyopenms.DriftTimeUnit.VSSC)
    call(spec, *_peaks(5), ion_mobility=_im(5))

    assert spec.size() == 5
    assert spec.getIMData()[1] == int(pyopenms.DriftTimeUnit.VSSC)
    assert np.allclose(spec.get_drift_time_array(), _im(5))
    assert len(spec.getFloatDataArrays()) == 1


@pytest.mark.parametrize("call", _CALL_FORMS)
def test_not_replacing_it_still_trips_the_metadata_guard(call):
    """The exemption is scoped to the array actually being replaced -- R-6.2 is unaffected."""
    spec = _spectrum(100, pyopenms.DriftTimeUnit.VSSC)
    with pytest.raises(RuntimeError, match="does not match the new spectrum size"):
        call(spec, *_peaks(5))
    assert spec.size() == 100


def test_other_arrays_are_still_governed_by_the_policy_when_ion_mobility_is_replaced():
    """Exempting the ion mobility array must not exempt everything else."""
    spec = _spectrum(100, pyopenms.DriftTimeUnit.VSSC)
    labels = pyopenms.StringDataArray()
    labels.set_data(["y%d" % i for i in range(100)])
    labels.setName("Peak annotation")
    spec.setStringDataArrays([labels])

    with pytest.raises(RuntimeError, match="StringDataArray"):
        spec.set_peaks(*_peaks(5), ion_mobility=_im(5))

    spec.set_peaks(*_peaks(5), ion_mobility=_im(5), metadata="clear")
    assert spec.size() == 5
    assert len(spec.getStringDataArrays()) == 0
    assert np.allclose(spec.get_drift_time_array(), _im(5))


def test_clear_policy_does_not_drop_the_array_installed_in_the_same_call():
    """"clear" runs before the new array is installed, so the two must not fight."""
    spec = _spectrum(100, pyopenms.DriftTimeUnit.MILLISECOND)
    spec.set_peaks(*_peaks(6), ion_mobility=_im(6),
                   ion_mobility_unit=pyopenms.DriftTimeUnit.MILLISECOND, metadata="clear")

    assert spec.containsIMData()
    assert len(spec.getFloatDataArrays()) == 1
    assert np.allclose(spec.get_drift_time_array(), _im(6))


@pytest.mark.parametrize("unit,_cv", _SUPPORTED_UNITS)
def test_mzml_roundtrip_preserves_ion_mobility(unit, _cv, tmp_path):
    """The whole point of naming from the CV: the unit has to survive a file round-trip."""
    spec = pyopenms.MSSpectrum()
    spec.setRT(1.0)
    spec.setMSLevel(1)
    spec.set_peaks(*_peaks(5), ion_mobility=_im(5), ion_mobility_unit=unit)
    exp = pyopenms.MSExperiment()
    exp.addSpectrum(spec)

    path = str(tmp_path / "im.mzML")
    pyopenms.MzMLFile().store(path, exp)
    loaded = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(path, loaded)

    back = loaded[0]
    assert back.containsIMData()
    assert back.getIMData()[1] == int(unit)
    assert np.allclose(back.get_drift_time_array(), _im(5))


def test_written_array_is_readable_through_the_zero_copy_view():
    spec = _spectrum(5, pyopenms.DriftTimeUnit.VSSC)
    view = spec.drift_time_array_view()
    assert np.allclose(view, _im(5))
    view[0] = 99.0
    assert spec.get_drift_time_array()[0] == pytest.approx(99.0)


def test_sortByPosition_carries_the_written_array_along():
    """An array written this way is a normal aligned float array, so sorting permutes it."""
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(np.array([300.0, 100.0, 200.0]), np.array([3.0, 1.0, 2.0], dtype=np.float32),
                   ion_mobility=np.array([0.3, 0.1, 0.2], dtype=np.float32),
                   ion_mobility_unit=pyopenms.DriftTimeUnit.VSSC)
    spec.sortByPosition()

    assert np.allclose(spec.get_peaks()[0], [100.0, 200.0, 300.0])
    assert np.allclose(spec.get_drift_time_array(), [0.1, 0.2, 0.3])


def test_replacing_the_array_does_not_free_a_live_view():
    """Replacing an equal-length IM array must not deallocate the buffer under a live view.

    ``FloatDataArray`` derives from ``std::vector<float>``, so move-assigning a replacement over
    it freed the destination's storage -- while ``drift_time_array_view()`` had already handed
    numpy a pointer into exactly that storage, keeping only the spectrum alive and not the buffer.
    Reading the view afterwards was a use-after-free with the ownership graph fully intact.
    """
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(np.array([100.0, 200.0, 300.0]), np.array([1.0, 2.0, 3.0], dtype=np.float32))
    spec.set_drift_time_array(np.array([1.0, 2.0, 3.0], dtype=np.float32),
                              pyopenms.DriftTimeUnit.MILLISECOND)

    view = spec.drift_time_array_view()
    spec.set_drift_time_array(np.array([4.0, 5.0, 6.0], dtype=np.float32),
                              pyopenms.DriftTimeUnit.MILLISECOND)

    # the view stays valid and observes the new values, as it would for any in-place write
    assert np.allclose(view, [4.0, 5.0, 6.0])
    assert np.allclose(spec.get_drift_time_array(), [4.0, 5.0, 6.0])


def test_replacing_the_array_still_replaces_the_unit():
    """The in-place path must carry the metadata across, not just the values."""
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(np.array([100.0, 200.0]), np.array([1.0, 2.0], dtype=np.float32))
    spec.set_drift_time_array(np.array([1.0, 2.0], dtype=np.float32),
                              pyopenms.DriftTimeUnit.MILLISECOND)
    assert spec.get_drift_time_array_unit() == pyopenms.DriftTimeUnit.MILLISECOND

    spec.set_drift_time_array(np.array([3.0, 4.0], dtype=np.float32),
                              pyopenms.DriftTimeUnit.VSSC)
    assert spec.get_drift_time_array_unit() == pyopenms.DriftTimeUnit.VSSC
    assert np.allclose(spec.get_drift_time_array(), [3.0, 4.0])


def test_replacing_with_a_different_length_still_works():
    """A length change must replace the storage; only the equal-length path can reuse it."""
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(np.array([100.0, 200.0, 300.0]), np.array([1.0, 2.0, 3.0], dtype=np.float32))
    spec.set_drift_time_array(np.array([1.0, 2.0, 3.0], dtype=np.float32),
                              pyopenms.DriftTimeUnit.MILLISECOND)

    spec.set_peaks(np.array([100.0, 200.0]), np.array([1.0, 2.0], dtype=np.float32),
                   ion_mobility=np.array([7.0, 8.0], dtype=np.float32),
                   ion_mobility_unit=pyopenms.DriftTimeUnit.MILLISECOND)

    assert np.allclose(spec.get_drift_time_array(), [7.0, 8.0])
