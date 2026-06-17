"""
Regression tests for mutable reference returns and reference_internal rv_policy.

These tests verify that methods returning T& actually return mutable references
(not copies), and that __getitem__ on container classes returns references tied
to the container's lifetime.
"""

import pytest


def test_msexperiment_getitem_returns_mutable_reference():
    """__getitem__ on MSExperiment must return a reference, not a copy.

    Regression: non-vector container __getitem__ was missing
    rv_policy::reference_internal, causing operator[] to return a copy.
    Modifying the copy had no effect on the container.
    """
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    # Modify via __getitem__ — must mutate the actual spectrum in exp
    exp[0].setRT(99.0)
    assert exp[0].getRT() == pytest.approx(99.0), (
        "__getitem__ returned a copy instead of a mutable reference"
    )


def test_msexperiment_getspectra_returns_mutable_reference():
    """getSpectra() must return a mutable reference to the internal vector."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    spectra = exp.getSpectra()
    spectra[0].setRT(42.0)
    assert exp.getSpectra()[0].getRT() == pytest.approx(42.0), (
        "getSpectra() returned a copy instead of a mutable reference"
    )


def test_peptide_identification_gethits_mutable():
    """getHits() must return a mutable reference to the internal hits vector."""
    from pyopenms import PeptideIdentification, PeptideHit

    pid = PeptideIdentification()
    hit = PeptideHit()
    hit.setScore(1.0)
    pid.insertHit(hit)

    hits = pid.getHits()
    hits[0].setScore(99.0)
    assert pid.getHits()[0].getScore() == pytest.approx(99.0), (
        "getHits() returned a copy instead of a mutable reference"
    )


def test_protein_identification_gethits_mutable():
    """getHits() must return a mutable reference to the internal hits vector."""
    from pyopenms import ProteinIdentification, ProteinHit

    prot = ProteinIdentification()
    hit = ProteinHit()
    hit.setScore(1.0)
    prot.insertHit(hit)

    hits = prot.getHits()
    hits[0].setScore(99.0)
    assert prot.getHits()[0].getScore() == pytest.approx(99.0), (
        "getHits() returned a copy instead of a mutable reference"
    )


def test_featuremap_getitem_returns_mutable_reference():
    """__getitem__ on FeatureMap must return a reference, not a copy."""
    from pyopenms import FeatureMap, Feature

    fm = FeatureMap()
    f = Feature()
    f.setRT(1.0)
    fm.push_back(f)

    fm[0].setRT(77.0)
    assert fm[0].getRT() == pytest.approx(77.0), (
        "FeatureMap __getitem__ returned a copy instead of a mutable reference"
    )


def test_consensusmap_getitem_returns_mutable_reference():
    """__getitem__ on ConsensusMap must return a reference, not a copy."""
    from pyopenms import ConsensusMap, ConsensusFeature

    cm = ConsensusMap()
    cf = ConsensusFeature()
    cf.setRT(1.0)
    cm.push_back(cf)

    cm[0].setRT(55.0)
    assert cm[0].getRT() == pytest.approx(55.0), (
        "ConsensusMap __getitem__ returned a copy instead of a mutable reference"
    )
