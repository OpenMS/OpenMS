"""
Regression tests for element-access ownership semantics.

Policy (see src/pyOpenMS/OWNERSHIP.md and issue #9792): indexing and
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


# --------------------------------------------------------------------------
# Containers found later: same policy, one with a sharper motivation
# --------------------------------------------------------------------------

def test_aasequence_residue_access_does_not_corrupt_residuedb():
    """AASequence stores const Residue* into the process-lifetime ResidueDB.

    A reference here aliased a database entry that nanobind exposes as mutable,
    so editing it changed the Residue for every sequence in the process --
    including ones constructed later. Copies remove the route entirely.
    """
    from pyopenms import AASequence

    seq = AASequence.fromString("PEPTIDE")
    res = seq[0]
    original = res.getName()
    res.setName("CORRUPTED")

    assert seq[0].getName() == original, "AASequence[i] aliased the ResidueDB entry"
    # a sequence built afterwards must be unaffected too -- this is the part that
    # proves the shared database itself was not edited
    assert AASequence.fromString("PEPTIDE")[0].getName() == original


def test_aasequence_iteration_returns_copies():
    from pyopenms import AASequence

    seq = AASequence.fromString("PEPTIDE")
    original = [r.getName() for r in seq]
    for r in seq:
        r.setName("CORRUPTED")
    assert [r.getName() for r in seq] == original


def test_isotopedistribution_getitem_returns_copy():
    from pyopenms import IsotopeDistribution

    dist = IsotopeDistribution()
    dist.insert(100.0, 1.0)
    dist.insert(101.0, 0.5)

    # a default-constructed distribution already holds one element, so the first
    # inserted peak is at index 1
    original = dist[1].getMZ()
    peak = dist[1]
    peak.setMZ(999.0)
    assert dist[1].getMZ() == pytest.approx(original), "IsotopeDistribution[i] aliased"


def test_acquisitioninfo_getitem_returns_copy():
    from pyopenms import AcquisitionInfo, Acquisition

    info = AcquisitionInfo()
    acq = Acquisition()
    acq.setMetaValue("key", 1)
    info.push_back(acq)

    got = info[0]
    got.setMetaValue("key", 42)
    assert info[0].getMetaValue("key") == 1, "AcquisitionInfo[i] aliased"

    info[0] = got  # write-back is the supported route
    assert info[0].getMetaValue("key") == 42


def test_peptide_identification_list_getitem_returns_copy():
    from pyopenms import PeptideIdentificationList, PeptideIdentification

    lst = PeptideIdentificationList()
    pid = PeptideIdentification()
    pid.setRT(1.0)
    lst.push_back(pid)

    got = lst[0]
    got.setRT(99.0)
    assert lst[0].getRT() == pytest.approx(1.0), "PeptideIdentificationList[i] aliased"

    lst[0] = got
    assert lst[0].getRT() == pytest.approx(99.0)


# ---------------------------------------------------------------------------
# Every getter that hands out a copy needs a route back in. These pin the
# write-back paths that were missing when the getters were converted.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cls_name", ["DTA2DFile", "ConsensusXMLFile", "MzMLFile",
                                      "MzXMLFile", "MzDataFile", "FileHandler"])
def test_file_options_round_trip(cls_name):
    import pyopenms

    handler = getattr(pyopenms, cls_name)()

    opts = handler.getOptions()
    opts.setMSLevels([2])
    handler.setOptions(opts)
    assert list(handler.getOptions().getMSLevels()) == [2]

    # the getter itself is a copy -- editing it must not reach the file object
    stray = handler.getOptions()
    stray.setMSLevels([3])
    assert list(handler.getOptions().getMSLevels()) == [2], \
        f"{cls_name}.getOptions() aliased"


@pytest.mark.parametrize("kind,ref", [("Protein", "P1"), ("Peptide", "PEP1"),
                                      ("Compound", "C1")])
def test_targeted_experiment_by_ref_returns_copy(kind, ref):
    import pyopenms

    exp = pyopenms.TargetedExperiment()
    entry = getattr(pyopenms, kind)()
    entry.id = ref
    getattr(exp, f"set{kind}s")([entry])

    got = getattr(exp, f"get{kind}ByRef")(ref)
    got.id = "MUTATED"
    assert getattr(exp, f"get{kind}ByRef")(ref).id == ref, \
        f"get{kind}ByRef aliased the stored entry"

    # write-back through the list setter is the supported route
    entries = getattr(exp, f"get{kind}s")()
    entries[0].id = "CHANGED"
    getattr(exp, f"set{kind}s")(entries)
    assert getattr(exp, f"get{kind}ByRef")("CHANGED").id == "CHANGED"


def test_mass_traces_supports_write_back():
    from pyopenms import MassTraces

    traces = MassTraces()
    assert len(traces) == 0

    # MassTrace (the FeatureFinderAlgorithmPicked helper struct) is not wrapped,
    # so no element can be handed in or out today. What this pins is that the
    # container protocol stays symmetric: a copy-returning __getitem__ is always
    # paired with a __setitem__, so a caller who can get an element can put it back.
    assert hasattr(MassTraces, "__getitem__")
    assert hasattr(MassTraces, "__setitem__")
    assert "trace" in MassTraces.__setitem__.__doc__
    with pytest.raises(IndexError):
        traces[0]


def test_data_processing_entries_are_shared_pointers():
    """The one documented aliasing carve-out (see OWNERSHIP.md).

    getDataProcessing() copies the list, but its entries are shared_ptrs, so an
    in-place edit of an entry is visible through the owner. Pinned here so the
    documented behaviour and the actual behaviour cannot drift apart.
    """
    from pyopenms import MSSpectrum, DataProcessing

    spec = MSSpectrum()
    dp = DataProcessing()
    dp.setMetaValue("note", "before")
    spec.setDataProcessing([dp])

    entries = spec.getDataProcessing()
    entries.append(DataProcessing())
    assert len(spec.getDataProcessing()) == 1, "the returned list is not a copy"

    entries = spec.getDataProcessing()
    entries[0].setMetaValue("note", "after")
    assert spec.getDataProcessing()[0].getMetaValue("note") == "after", \
        "shared_ptr entries no longer alias -- OWNERSHIP.md needs updating"


# ---------------------------------------------------------------------------
# Getters bound as &Class::method (rather than a lambda) were a separate shape
# that the first conversion sweeps did not reach. These pin the ones that used
# to hand out live aliases.
# ---------------------------------------------------------------------------

def test_mrm_transition_group_write_back():
    """getFeaturesMuteable() & friends copy; the matching setters write back."""
    from pyopenms import (MRMTransitionGroupCP, MRMFeature, MSChromatogram,
                          ReactionMonitoringTransition)

    group = MRMTransitionGroupCP()

    feature = MRMFeature()
    feature.setMetaValue("score", 0.0)
    group.addFeature(feature)
    features = group.getFeaturesMuteable()
    features[0].setMetaValue("score", 1.0)
    assert group.getFeaturesMuteable()[0].getMetaValue("score") == 0.0, \
        "getFeaturesMuteable() aliased the stored features"
    group.setFeatures(features)
    assert group.getFeaturesMuteable()[0].getMetaValue("score") == 1.0

    chrom = MSChromatogram()
    chrom.setNativeID("c1")
    group.addChromatogram(chrom, "c1")
    chroms = group.getChromatograms()
    chroms[0].setNativeID("edited")
    assert group.getChromatograms()[0].getNativeID() == "c1"
    group.setChromatograms(chroms)
    assert group.getChromatograms()[0].getNativeID() == "edited"

    transition = ReactionMonitoringTransition()
    transition.setName("t1")
    group.addTransition(transition, "t1")
    transitions = group.getTransitionsMuteable()
    transitions[0].setName("edited")
    assert group.getTransitionsMuteable()[0].getName() == "t1"
    group.setTransitions(transitions)
    assert group.getTransitionsMuteable()[0].getName() == "edited"


def test_mrm_transition_group_write_back_keeps_lookup_keys():
    """The keyed collections carry a key -> index map the setters must not invalidate."""
    from pyopenms import MRMTransitionGroupCP, MSChromatogram

    group = MRMTransitionGroupCP()
    chrom = MSChromatogram()
    chrom.setNativeID("native_id")
    group.addChromatogram(chrom, "custom_key")   # key deliberately differs from the native ID

    edited = group.getChromatograms()
    edited[0].setMetaValue("k", 1)
    group.setChromatograms(edited)

    assert group.hasChromatogram("custom_key"), "the caller's key was dropped"
    assert not group.hasChromatogram("native_id"), "the collection was silently re-keyed"
    assert group.getChromatogram("custom_key").getMetaValue("k") == 1

    # A size change would leave the key map pointing at the wrong (or no) element.
    for wrong_size in ([], list(edited) + list(edited)):
        with pytest.raises(ValueError):
            group.setChromatograms(wrong_size)
    assert len(group.getChromatograms()) == 1 and group.hasChromatogram("custom_key")


def test_terminal_modifications_do_not_alias_modifications_db():
    from pyopenms import AASequence

    plain = AASequence.fromString("PEPTIDE")
    assert plain.getNTerminalModification() is None
    assert plain.getCTerminalModification() is None

    seq = AASequence.fromString(".(Acetyl)PEPTIDE")
    mod = seq.getNTerminalModification()
    assert mod is not None
    original = mod.getDiffMonoMass()
    mod.setDiffMonoMass(999.0)
    assert AASequence.fromString(".(Acetyl)PEPTIDE").getNTerminalModification().getDiffMonoMass() \
        == pytest.approx(original), "editing the returned modification reached ModificationsDB"


def test_residue_modification_is_a_copy():
    from pyopenms import AASequence

    assert AASequence.fromString("PEPTIDE")[0].getModification() is None
    residue = AASequence.fromString("PEM(Oxidation)TIDE")[2]
    mod = residue.getModification()
    assert mod is not None

    original = mod.getDiffMonoMass()
    mod.setDiffMonoMass(999.0)
    assert AASequence.fromString("PEM(Oxidation)TIDE")[2].getModification().getDiffMonoMass() \
        == pytest.approx(original), "editing the returned modification reached ModificationsDB"


def test_nasequence_getsequence_does_not_alias_ribonucleotide_db():
    from pyopenms import NASequence

    seq = NASequence.fromString("AAUC")
    ribos = seq.getSequence()
    assert len(ribos) == 4
    original = NASequence.fromString("AAUC").getSequence()[0].getName()
    ribos[0].setName("CORRUPTED")
    assert NASequence.fromString("AAUC").getSequence()[0].getName() == original, \
        "editing a returned ribonucleotide reached RibonucleotideDB"


def test_param_entry_and_node_finders_return_copies():
    from pyopenms import Param, ParamNode

    param = Param()
    param.setValue("a:b", 1)
    entry = param.getEntry("a:b")
    entry.description = "changed"
    assert param.getEntry("a:b").description != "changed", "Param.getEntry() aliased"

    # both finders returned raw pointers before; a miss is None, not a dangling alias
    node = ParamNode()
    assert node.findEntryRecursive("missing") is None
    assert node.findParentOf("missing") is None


def test_imaging_geometry_getters_return_copies():
    from pyopenms import MSImagingGeometry, MSImagingRegion

    geom = MSImagingGeometry()
    geom.addRegion(MSImagingRegion.rectangle(1, "r", 0, 0, 1, 1))
    assert geom.getRegions()[0] is not geom.getRegions()[0]
    assert geom.getRegion(1) is not geom.getRegion(1)


# ---------------------------------------------------------------------------
# NASequence element access: the stored elements are pointers into the
# process-wide RibonucleotideDB, so every accessor has to copy, and the
# setters resolve copies back to database entries by code.
# ---------------------------------------------------------------------------

def test_nasequence_element_access_does_not_alias_ribonucleotide_db():
    from pyopenms import NASequence

    original = NASequence.fromString("AAUC")[0].getName()

    seq = NASequence.fromString("AAUC")
    seq[0].setName("CORRUPTED")                       # __getitem__
    seq.get(0).setName("CORRUPTED")                   # get()
    for ribo in seq:                                  # __iter__
        ribo.setName("CORRUPTED")
    assert NASequence.fromString("AAUC")[0].getName() == original, \
        "editing an accessed ribonucleotide reached RibonucleotideDB"
    assert seq[0].getName() == original


def test_nasequence_write_back_round_trip():
    """setSequence(getSequence()) must work: setters resolve copies by code."""
    from pyopenms import NASequence, Ribonucleotide

    seq = NASequence.fromString("AAUC")
    seq.setSequence(seq.getSequence())
    assert seq.toString() == "AAUC"

    seq.set(0, seq.getSequence()[2])                  # put a copied 'U' at position 0
    assert seq.toString() == "UAUC"

    # a ribonucleotide whose code is not in the database is still rejected
    with pytest.raises(ValueError):
        seq.set(0, Ribonucleotide())


def test_nasequence_terminal_mods_are_copies():
    from pyopenms import NASequence

    assert NASequence.fromString("AAUC").getFivePrimeMod() is None
    assert NASequence.fromString("AAUC").getThreePrimeMod() is None


# ---------------------------------------------------------------------------
# MSExperiment.getSpectrum(i)/getChromatogram(i) wrap an unchecked C++
# spectra_[i]; the copy conversion made an out-of-bounds read unconditional,
# so the binding has to bounds-check like __getitem__ does.
# ---------------------------------------------------------------------------

def test_msexperiment_named_getters_bounds_checked():
    from pyopenms import MSExperiment

    exp = MSExperiment()
    with pytest.raises(IndexError):
        exp.getSpectrum(0)
    with pytest.raises(IndexError):
        exp.getChromatogram(0)


# ---------------------------------------------------------------------------
# MSImagingExperiment: the per-pixel getter copies, so a per-pixel setter is
# the write-back route.
# ---------------------------------------------------------------------------

def test_imaging_experiment_per_pixel_write_back():
    from pyopenms import (MSExperiment, MSSpectrum, MSImagingExperiment,
                          MSImagingGeometry)

    exp = MSExperiment()
    for _ in range(3):
        spec = MSSpectrum()
        spec.setRT(1.0)
        exp.addSpectrum(spec)
    geom = MSImagingGeometry()
    geom.setDimensions(2, 2)
    geom.addPixel(0, 0, 0)
    geom.addPixel(1, 0, 1)
    geom.addPixel(0, 1, 2)
    mie = MSImagingExperiment(exp)
    mie.setGeometry(geom)

    got = mie.getSpectrum(0, 1)
    got.setRT(99.0)
    assert mie.getSpectrum(0, 1).getRT() == 1.0, "per-pixel getter aliases"
    mie.setSpectrum(0, 1, got)
    assert mie.getSpectrum(0, 1).getRT() == 99.0
    assert mie.getSpectrum(0, 0).getRT() == 1.0, "wrong pixel touched"
    with pytest.raises(Exception):
        mie.setSpectrum(1, 1, got)                    # pixel not in the geometry


# ---------------------------------------------------------------------------
# SimpleOpenMSSpectraFactory: probes the cached_data marker C++-side (no
# per-spectrum copies) and both branches resolve their class names. The
# factory was unreachable before -- both branches raised NameError -- so this
# is the first test to actually call it.
# ---------------------------------------------------------------------------

def test_spectrum_access_factory():
    from pyopenms import (SimpleOpenMSSpectraFactory, MSExperiment, MSSpectrum,
                          SpectrumAccessOpenMS, DataProcessing)

    exp = MSExperiment()
    exp.addSpectrum(MSSpectrum())
    assert exp._contains_cached_data_marker() is False
    access = SimpleOpenMSSpectraFactory.getSpectrumAccessOpenMSPtr(exp)
    assert isinstance(access, SpectrumAccessOpenMS)

    marked = MSSpectrum()
    processing = DataProcessing()
    processing.setMetaValue("cached_data", 1)
    marked.setDataProcessing([processing])
    exp.addSpectrum(marked)
    assert exp._contains_cached_data_marker() is True


def test_spectrum_access_factory_cached_branch(tmp_path):
    """The factory returns the disk-backed accessor for cached experiments.

    Full round trip: cache an experiment with CachedmzML.store (which stamps
    the 'cached_data' marker into the metadata file), reload the metadata,
    and check the factory hands back a working SpectrumAccessOpenMSCached
    that decodes the original peaks from the .cached side-car.
    """
    from pyopenms import (SimpleOpenMSSpectraFactory, MSExperiment, MSSpectrum,
                          SpectrumAccessOpenMSCached, CachedmzML, MzMLFile)

    exp = MSExperiment()
    spec = MSSpectrum()
    spec.setRT(41.5)
    spec.set_peaks(([100.0, 200.0, 300.0], [10.0, 20.0, 30.0]))
    exp.addSpectrum(spec)

    meta_file = str(tmp_path / "cached_run.mzML")
    CachedmzML.store(meta_file, exp)

    reloaded = MSExperiment()
    MzMLFile().load(meta_file, reloaded)
    assert reloaded._contains_cached_data_marker() is True

    access = SimpleOpenMSSpectraFactory.getSpectrumAccessOpenMSPtr(reloaded)
    assert isinstance(access, SpectrumAccessOpenMSCached)
    assert access.getNrSpectra() == 1

    decoded = access.getSpectrumById(0)
    assert list(decoded.get_mz_array()) == [100.0, 200.0, 300.0]
    assert list(decoded.get_intensity_array()) == [10.0, 20.0, 30.0]
    assert access.getSpectraByRT(41.5, 1.0) == [0]
