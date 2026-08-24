# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
"""The explicit `_view` family: opt-in aliasing access (see OWNERSHIP.md).

Every plain getter returns a copy; the `_view`-suffixed accessors are the one
deliberate way to alias an object's storage, and the name itself carries the
warning. Naming scheme under test:

* singular + ``_view``            -> one live element   (``spectrum_view(i)``)
* plural + ``_view``              -> list of live elements (``spectrum_views()``)
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
    exp.spectrum_view(0).setRT(50.0)
    assert exp[0].getRT() == 50.0
    assert exp[1].getRT() == 2.0


def test_spectrum_view_identity_is_stable():
    """Asking twice for the same element yields the same wrapper."""
    exp = _experiment()
    view = exp.spectrum_view(0)
    assert exp.spectrum_view(0) is view


def test_spectrum_view_keeps_parent_alive():
    exp = _experiment(rts=(7.0,))
    view = exp.spectrum_view(0)
    del exp
    gc.collect()
    assert view.getRT() == 7.0


def test_spectrum_view_bounds_checked():
    with pytest.raises(IndexError):
        _experiment().spectrum_view(99)


def test_spectra_view_and_iterator_edit_all_elements():
    exp = _experiment()
    for view in exp.spectrum_views():
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

    exp.chromatogram_view(0).setNativeID("edited")
    assert exp.getChromatogram(0).getNativeID() == "edited"

    for view in exp.chromatogram_views():
        view.setNativeID("again")
    for view in exp.iter_chromatogram_views():
        assert view.getNativeID() == "again"
    with pytest.raises(IndexError):
        exp.chromatogram_view(99)


def test_feature_views():
    fm = pyopenms.FeatureMap()
    feature = pyopenms.Feature()
    feature.setRT(5.0)
    fm.push_back(feature)

    fm.feature_view(0).setRT(9.0)
    assert fm[0].getRT() == 9.0
    fm.feature_views()[0].setIntensity(3.0)
    assert fm[0].getIntensity() == 3.0
    for view in fm.iter_feature_views():
        view.setRT(1.0)
    assert fm[0].getRT() == 1.0
    with pytest.raises(IndexError):
        fm.feature_view(99)


def test_views_do_not_change_the_copy_getters():
    """The `_view` family is additive: everything unsuffixed still copies."""
    exp = _experiment()
    held = exp[0]
    held.setRT(999.0)
    assert exp[0].getRT() == 1.0
    held = exp.getSpectrum(0)
    held.setRT(999.0)
    assert exp.getSpectrum(0).getRT() == 1.0


def test_consensus_feature_views():
    cm = pyopenms.ConsensusMap()
    cf = pyopenms.ConsensusFeature()
    cf.setRT(1.0)
    cm.push_back(cf)

    cm.consensus_feature_view(0).setRT(9.0)
    assert cm[0].getRT() == 9.0
    for view in cm.iter_consensus_feature_views():
        view.setIntensity(5.0)
    assert cm.consensus_feature_views()[0].getIntensity() == 5.0
    with pytest.raises(IndexError):
        cm.consensus_feature_view(99)


def test_peptide_identification_views():
    peps = pyopenms.PeptideIdentificationList()
    pid = pyopenms.PeptideIdentification()
    pid.setRT(1.0)
    peps.push_back(pid)

    peps.peptide_identification_view(0).setRT(9.0)
    assert peps[0].getRT() == 9.0
    assert len(peps.peptide_identification_views()) == 1
    for view in peps.iter_peptide_identification_views():
        view.setScoreType("s")
    assert peps[0].getScoreType() == "s"
    with pytest.raises(IndexError):
        peps.peptide_identification_view(99)


def test_data_array_views_zero_copy_chain():
    """spec.float_data_array_view(i).data_view() writes spectrum storage."""
    import numpy as np

    spec = pyopenms.MSSpectrum()
    fda = pyopenms.FloatDataArray()
    fda.set_data(np.array([1.0, 2.0], dtype=np.float32))
    fda.setName("im")
    spec.setFloatDataArrays([fda])

    view = spec.float_data_array_view(0)
    view.setName("edited")
    assert spec.getFloatDataArrays()[0].getName() == "edited"

    arr = view.data_view()
    arr[0] = 42.0
    assert spec.getFloatDataArrays()[0].get_data()[0] == 42.0

    assert len(spec.float_data_array_views()) == 1
    assert spec._integer_data_array_count() == 0
    with pytest.raises(IndexError):
        spec.float_data_array_view(99)


def test_chromatogram_data_array_views():
    import numpy as np

    chrom = pyopenms.MSChromatogram()
    fda = pyopenms.FloatDataArray()
    fda.set_data(np.array([1.0], dtype=np.float32))
    chrom.setFloatDataArrays([fda])

    chrom.float_data_array_view(0).setName("edited")
    assert chrom.getFloatDataArrays()[0].getName() == "edited"
    for view in chrom.iter_float_data_array_views():
        assert view.getName() == "edited"


def test_mrm_transition_group_views():
    group = pyopenms.MRMTransitionGroupCP()
    feature = pyopenms.MRMFeature()
    feature.setMetaValue("k", 1)
    group.addFeature(feature)
    chrom = pyopenms.MSChromatogram()
    chrom.setNativeID("c1")
    group.addChromatogram(chrom, "c1")

    group.feature_view(0).setMetaValue("k", 2)
    assert group.getFeatures()[0].getMetaValue("k") == 2
    group.chromatogram_view(0).setNativeID("edited")
    assert group.getChromatogram("c1").getNativeID() == "edited"
    assert len(group.feature_views()) == 1
    assert len(group.chromatogram_views()) == 1
    with pytest.raises(IndexError):
        group.feature_view(99)
