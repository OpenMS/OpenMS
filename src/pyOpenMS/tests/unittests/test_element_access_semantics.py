"""
Regression tests for element-access ownership semantics.

Policy (see src/pyOpenMS/OWNERSHIP_PROPOSAL.md and issue #9792): indexing and
iterating a container returns an **owned value**, never a live alias into the
container's storage. Writes go back through ``__setitem__``.

This restores the semantics of the Cython bindings up to release 3.5, where the
wrapper representation forced a copy on every wrapped-type return. The nanobind
port made these paths aliases, which is what produced the use-after-free and
wrong-element-write failures catalogued in #9792.

The final section pins the surfaces that intentionally still alias, so the
boundary of the policy is explicit rather than accidental.
"""

import pytest


# --------------------------------------------------------------------------
# Element access returns owned values
# --------------------------------------------------------------------------

def test_msexperiment_getitem_returns_copy():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    exp[0].setRT(99.0)  # edits a detached copy
    assert exp[0].getRT() == pytest.approx(1.0), (
        "__getitem__ aliased the container instead of returning a copy"
    )


def test_msexperiment_iteration_returns_copies():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    for rt in (1.0, 2.0):
        spec = MSSpectrum()
        spec.setRT(rt)
        exp.addSpectrum(spec)

    for spec in exp:
        spec.setRT(99.0)

    assert [exp[i].getRT() for i in range(2)] == pytest.approx([1.0, 2.0]), (
        "iteration aliased the container instead of yielding copies"
    )


def test_msexperiment_getspectra_returns_copies():
    """getSpectra() returned a copy in the Cython bindings (declared by value)."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    spectra = exp.getSpectra()
    spectra[0].setRT(42.0)
    assert exp[0].getRT() == pytest.approx(1.0), (
        "getSpectra() handed out aliases instead of copies"
    )


def test_msspectrum_peak_access_is_consistent():
    """Indexing and iteration must agree: both copy."""
    from pyopenms import MSSpectrum

    spec = MSSpectrum()
    spec.set_peaks(([100.0, 200.0], [10.0, 20.0]))

    spec[0].setMZ(999.0)
    assert spec[0].getMZ() == pytest.approx(100.0), "spec[i] aliased"

    for peak in spec:
        peak.setMZ(999.0)
    assert spec[0].getMZ() == pytest.approx(100.0), "iteration aliased"


def test_featuremap_getitem_returns_copy():
    from pyopenms import FeatureMap, Feature

    fm = FeatureMap()
    f = Feature()
    f.setRT(1.0)
    fm.push_back(f)

    fm[0].setRT(77.0)
    assert fm[0].getRT() == pytest.approx(1.0), "FeatureMap __getitem__ aliased"


def test_consensusmap_getitem_returns_copy():
    from pyopenms import ConsensusMap, ConsensusFeature

    cm = ConsensusMap()
    cf = ConsensusFeature()
    cf.setRT(1.0)
    cm.push_back(cf)

    cm[0].setRT(55.0)
    assert cm[0].getRT() == pytest.approx(1.0), "ConsensusMap __getitem__ aliased"


def test_mschromatogram_getitem_returns_copy():
    from pyopenms import MSChromatogram, ChromatogramPeak

    chrom = MSChromatogram()
    p = ChromatogramPeak()
    p.setRT(1.0)
    p.setIntensity(10.0)
    chrom.push_back(p)

    chrom[0].setIntensity(999.0)
    assert chrom[0].getIntensity() == pytest.approx(10.0), "MSChromatogram __getitem__ aliased"


def test_mobilogram_getitem_returns_copy():
    from pyopenms import Mobilogram, MobilityPeak1D

    mob = Mobilogram()
    p = MobilityPeak1D()
    p.setMobility(1.0)
    p.setIntensity(10.0)
    mob.push_back(p)

    mob[0].setIntensity(999.0)
    assert mob[0].getIntensity() == pytest.approx(10.0), "Mobilogram __getitem__ aliased"


# --------------------------------------------------------------------------
# Writes go back through __setitem__
# --------------------------------------------------------------------------

def test_write_back_via_setitem():
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    edited = exp[0]
    edited.setRT(55.0)
    exp[0] = edited
    assert exp[0].getRT() == pytest.approx(55.0)


def test_write_back_via_setspectra():
    """The copy-edit-write-back idiom used by test_AcquisitionInfo."""
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(1.0)
    exp.addSpectrum(spec)

    spectra = exp.getSpectra()
    spectra[0].setRT(33.0)
    exp.setSpectra(spectra)
    assert exp[0].getRT() == pytest.approx(33.0)


# --------------------------------------------------------------------------
# The failures this policy exists to prevent (#9792)
# --------------------------------------------------------------------------

def test_held_spectrum_survives_container_growth():
    """#9792 L2: growing the experiment reallocates its storage.

    Under aliasing, a previously obtained spectrum pointed into the freed block
    and silently reported zero peaks. A copy is unaffected.
    """
    from pyopenms import MSExperiment, MSSpectrum

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(42.0)
    spec.set_peaks(([100.0, 200.0, 300.0], [1.0, 2.0, 3.0]))
    exp.addSpectrum(spec)

    held = exp[0]
    for i in range(2000):
        s = MSSpectrum()
        s.setRT(float(i))
        exp.addSpectrum(s)

    assert held.getRT() == pytest.approx(42.0)
    assert held.size() == 3, "held spectrum lost its peaks (use-after-free)"


def test_held_feature_unaffected_by_sort():
    """#9792 L4: sorting reorders slots, so an alias silently re-labels itself."""
    from pyopenms import FeatureMap, Feature

    fm = FeatureMap()
    for rt, mz in ((30.0, 300.0), (10.0, 100.0), (20.0, 200.0)):
        f = Feature()
        f.setRT(rt)
        f.setMZ(mz)
        fm.push_back(f)

    held = fm[0]              # the RT=30 feature
    assert held.getRT() == pytest.approx(30.0)
    fm.sortByRT()
    assert held.getRT() == pytest.approx(30.0), (
        "held feature followed the slot instead of staying itself"
    )

    held.setMZ(999.0)         # must not land on whatever now occupies slot 0
    assert [fm[i].getMZ() for i in range(3)] == pytest.approx([100.0, 200.0, 300.0])


# --------------------------------------------------------------------------
# Vector-returning getters also hand back owned values
# --------------------------------------------------------------------------

def test_peptide_identification_gethits_returns_copies():
    from pyopenms import PeptideIdentification, PeptideHit

    pid = PeptideIdentification()
    hit = PeptideHit()
    hit.setScore(1.0)
    pid.insertHit(hit)

    hits = pid.getHits()
    hits[0].setScore(99.0)
    assert pid.getHits()[0].getScore() == pytest.approx(1.0), "getHits() aliased"

    pid.setHits(hits)  # write-back is the supported route
    assert pid.getHits()[0].getScore() == pytest.approx(99.0)


def test_protein_identification_gethits_returns_copies():
    from pyopenms import ProteinIdentification, ProteinHit

    prot = ProteinIdentification()
    hit = ProteinHit()
    hit.setScore(1.0)
    prot.insertHit(hit)

    hits = prot.getHits()
    hits[0].setScore(99.0)
    assert prot.getHits()[0].getScore() == pytest.approx(1.0), "getHits() aliased"

    prot.setHits(hits)
    assert prot.getHits()[0].getScore() == pytest.approx(99.0)


def test_float_data_arrays_return_copies():
    from pyopenms import MSSpectrum, FloatDataArray

    spec = MSSpectrum()
    spec.set_peaks(([100.0, 200.0], [10.0, 20.0]))
    fda = FloatDataArray()
    fda.setName("sn")
    fda.push_back(1.0)
    fda.push_back(2.0)
    spec.setFloatDataArrays([fda])

    arrays = spec.getFloatDataArrays()
    arrays[0].setName("changed")
    assert spec.getFloatDataArrays()[0].getName() == "sn", "getFloatDataArrays() aliased"

    spec.setFloatDataArrays(arrays)
    assert spec.getFloatDataArrays()[0].getName() == "changed"


def test_precursors_return_copies():
    from pyopenms import MSSpectrum, Precursor

    spec = MSSpectrum()
    p = Precursor()
    p.setMZ(500.0)
    spec.setPrecursors([p])

    precs = spec.getPrecursors()
    precs[0].setMZ(999.0)
    assert spec.getPrecursors()[0].getMZ() == pytest.approx(500.0), "getPrecursors() aliased"

    spec.setPrecursors(precs)
    assert spec.getPrecursors()[0].getMZ() == pytest.approx(999.0)


def test_feature_subordinates_return_copies():
    from pyopenms import Feature

    f = Feature()
    sub = Feature()
    sub.setRT(1.0)
    f.setSubordinates([sub])

    subs = f.getSubordinates()
    subs[0].setRT(42.0)
    assert f.getSubordinates()[0].getRT() == pytest.approx(1.0), "getSubordinates() aliased"

    f.setSubordinates(subs)
    assert f.getSubordinates()[0].getRT() == pytest.approx(42.0)
