# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
"""The explicit `_view` family: opt-in aliasing access (see OWNERSHIP.md).

Every plain getter returns a copy; the `_view`-suffixed accessors are the one
deliberate way to alias an object's storage, and the name itself carries the
warning. Naming scheme under test:

* singular + ``_view``            -> one live element   (``get_spectrum_view(i)``)
* plural + ``_view``              -> list of live elements (``get_spectra_view()``)
* ``iter_<singular>_views()``     -> iterator of live elements

Contract: edits land immediately; each view keeps its parent alive; a view is
invalidated by resizing or reordering the underlying list (that hazard is
undefined behaviour and deliberately not exercised here).
"""

import gc

import pytest

import pyopenms


def _experiment(rts=(1.0, 2.0)):
    exp = pyopenms.MSExperiment()
    for rt in rts:
        spec = pyopenms.MSSpectrum()
        spec.setRT(rt)
        exp.addSpectrum(spec)
    return exp


def test_spectrum_view_edit_lands_without_write_back():
    exp = _experiment()
    exp.get_spectrum_view(0).setRT(50.0)
    assert exp[0].getRT() == 50.0
    assert exp[1].getRT() == 2.0


def test_spectrum_view_identity_is_stable():
    """Asking twice for the same element yields the same wrapper."""
    exp = _experiment()
    view = exp.get_spectrum_view(0)
    assert exp.get_spectrum_view(0) is view


def test_spectrum_view_keeps_parent_alive():
    exp = _experiment(rts=(7.0,))
    view = exp.get_spectrum_view(0)
    del exp
    gc.collect()
    assert view.getRT() == 7.0


def test_spectrum_view_bounds_checked():
    with pytest.raises(IndexError):
        _experiment().get_spectrum_view(99)


def test_spectra_view_and_iterator_edit_all_elements():
    exp = _experiment()
    for view in exp.get_spectra_view():
        view.setRT(view.getRT() + 1.0)
    assert [exp[i].getRT() for i in range(2)] == [2.0, 3.0]

    for view in exp.iter_spectrum_views():
        view.setRT(0.0)
    assert [exp[i].getRT() for i in range(2)] == [0.0, 0.0]


def test_chromatogram_views():
    exp = _experiment()
    chrom = pyopenms.MSChromatogram()
    chrom.setNativeID("c1")
    exp.addChromatogram(chrom)

    exp.get_chromatogram_view(0).setNativeID("edited")
    assert exp.getChromatogram(0).getNativeID() == "edited"

    for view in exp.get_chromatograms_view():
        view.setNativeID("again")
    for view in exp.iter_chromatogram_views():
        assert view.getNativeID() == "again"
    with pytest.raises(IndexError):
        exp.get_chromatogram_view(99)


def test_feature_views():
    fm = pyopenms.FeatureMap()
    feature = pyopenms.Feature()
    feature.setRT(5.0)
    fm.push_back(feature)

    fm.get_feature_view(0).setRT(9.0)
    assert fm[0].getRT() == 9.0
    fm.get_features_view()[0].setIntensity(3.0)
    assert fm[0].getIntensity() == 3.0
    for view in fm.iter_feature_views():
        view.setRT(1.0)
    assert fm[0].getRT() == 1.0
    with pytest.raises(IndexError):
        fm.get_feature_view(99)


def test_views_do_not_change_the_copy_getters():
    """The `_view` family is additive: everything unsuffixed still copies."""
    exp = _experiment()
    held = exp[0]
    held.setRT(999.0)
    assert exp[0].getRT() == 1.0
    held = exp.getSpectrum(0)
    held.setRT(999.0)
    assert exp.getSpectrum(0).getRT() == 1.0
