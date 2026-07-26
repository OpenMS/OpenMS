"""Regression tests for borrowed-reference lifetime bugs (issue #9792, part V.4).

Two distinct defects are covered here.

**Wrong keep-alive target.**  ``PeakIndex.getFeature/getPeak/getSpectrum`` return a
reference *into the map argument*, but were annotated ``rv_policy::reference_internal``,
which keeps argument 1 -- the ``PeakIndex`` -- alive instead.  The map could therefore be
collected while the returned wrapper still pointed into it.  The edge has to target
argument 2.

**Retained raw pointers.**  Several setters hand OpenMS a raw pointer that it stores
verbatim and keeps using after the call.  ``nb::keep_alive`` cannot make those safe: the
receivers are copyable, and a C++ copy duplicates the stored pointer *without* carrying
the Python-side lifetime edge, so the copy dangles as soon as the original wrapper dies::

    d1 = EnzymaticDigestion()
    d1.setEnzyme(e)                 # even with keep_alive, the edge lives on d1
    d2 = EnzymaticDigestion(d1)     # C++ copy duplicates enzyme_ but not the edge
    del d1, e                       # edge released -> d2.enzyme_ dangles

The only sound fix is to guarantee the pointee has process lifetime.  Two mechanisms are
used, both verified below:

* **validation** -- reject a pointer that is not owned by the corresponding database
  (``RibonucleotideDB``, ``ProteaseDB``/``RNaseDB``, ``ResidueDB``, ``ModificationsDB``);
* **interning** -- route modification setters through the by-reference overloads, which
  copy the modification into ``ModificationsDB`` and store the database-owned copy.

``test_copy_outlives_original_*`` are the tests that specifically pin down the property a
``keep_alive``-based fix would *not* have provided.
"""

import copy
import gc
import sys

import pytest

import pyopenms


def _churn(n=20000):
    """Force collection and aggressively reuse freed heap blocks."""
    gc.collect()
    junk = [bytearray(b"\xa5" * n) for n in (16, 64, 256, 1024, 4096) for _ in range(200)]
    mods = []
    for i in range(n // 4):
        m = pyopenms.ResidueModification()
        m.setId("ZZZ%d" % i)
        m.setFullId("ZZZZZZZZZZZZ%d" % i)
        mods.append(m)
    del junk, mods
    gc.collect()


# ---------------------------------------------------------------------------
# PeakIndex: the keep-alive edge must target the map, not the index
# ---------------------------------------------------------------------------


def _feature_map(n=1):
    fm = pyopenms.FeatureMap()
    for i in range(n):
        f = pyopenms.Feature()
        f.setRT(11.0 + i)
        f.setMZ(222.0 + i)
        f.setIntensity(333.0 + i)
        fm.push_back(f)
    return fm


def _experiment(n_spec=1, n_peaks=4):
    exp = pyopenms.MSExperiment()
    for s in range(n_spec):
        spec = pyopenms.MSSpectrum()
        spec.setRT(10.0 + s)
        for i in range(n_peaks):
            pk = pyopenms.Peak1D()
            pk.setMZ(100.0 + i)
            pk.setIntensity(50.0 + i)
            spec.push_back(pk)
        exp.addSpectrum(spec)
    return exp


@pytest.mark.parametrize(
    "build_map, call",
    [
        (_feature_map, lambda idx, m: idx.getFeature(m)),
        (_experiment, lambda idx, m: idx.getPeak(m)),
        (_experiment, lambda idx, m: idx.getSpectrum(m)),
    ],
    ids=["getFeature", "getPeak", "getSpectrum"],
)
def test_peakindex_keeps_the_map_alive_not_the_index(build_map, call):
    """The returned wrapper must hold a reference to the *map*, not to the PeakIndex.

    This is a deterministic check: nanobind's keep_alive increments the patient's
    refcount, so a correctly targeted edge is directly observable.
    """
    container = build_map()
    idx = pyopenms.PeakIndex(0, 0)

    rc_map_before = sys.getrefcount(container)
    rc_idx_before = sys.getrefcount(idx)
    got = call(idx, container)

    assert sys.getrefcount(container) == rc_map_before + 1, (
        "returned wrapper does not keep the map alive; it would dangle once the "
        "map is collected"
    )
    # The index is not the owner of the referent and must not be the nurse.
    assert sys.getrefcount(idx) == rc_idx_before
    assert got is not None


@pytest.mark.parametrize(
    "build_map, call",
    [
        (_feature_map, lambda idx, m: idx.getFeature(m)),
        (_experiment, lambda idx, m: idx.getPeak(m)),
        (_experiment, lambda idx, m: idx.getSpectrum(m)),
    ],
    ids=["getFeature", "getPeak", "getSpectrum"],
)
def test_peakindex_rejects_out_of_range_index(build_map, call):
    """The accessors range-check only via OPENMS_PRECONDITION, which expands to nothing in a
    release build, so an out-of-range PeakIndex read past the end of the container."""
    container = build_map()
    with pytest.raises(IndexError):
        call(pyopenms.PeakIndex(5, 5), container)


@pytest.mark.parametrize(
    "build_map, call",
    [
        (pyopenms.FeatureMap, lambda idx, m: idx.getFeature(m)),
        (pyopenms.MSExperiment, lambda idx, m: idx.getPeak(m)),
        (pyopenms.MSExperiment, lambda idx, m: idx.getSpectrum(m)),
    ],
    ids=["getFeature", "getPeak", "getSpectrum"],
)
def test_peakindex_on_empty_container_raises(build_map, call):
    with pytest.raises(IndexError):
        call(pyopenms.PeakIndex(0, 0), build_map())


def test_featuremappinginfo_does_not_expose_the_removed_kd_tree_type():
    """`kd_tree` is a KDTreeFeatureMaps member; leaving it bound after removing that type made
    reading the attribute raise a confusing TypeError about an unconvertible return value."""
    fm_info = pyopenms.FeatureMapping_FeatureMappingInfo()
    assert isinstance(fm_info.feature_maps, list)
    assert not hasattr(fm_info, "kd_tree")


def test_peakindex_feature_survives_map_deletion():
    fm = _feature_map()
    idx = pyopenms.PeakIndex(0, 0)
    got = idx.getFeature(fm)
    assert got.getRT() == 11.0

    del fm
    _churn()

    assert got.getRT() == 11.0
    assert got.getMZ() == 222.0


def test_peakindex_spectrum_survives_experiment_deletion():
    exp = _experiment()
    idx = pyopenms.PeakIndex(0, 0)
    spec = idx.getSpectrum(exp)
    assert spec.getRT() == 10.0

    del exp
    _churn()

    assert spec.getRT() == 10.0
    assert spec.size() == 4


# ---------------------------------------------------------------------------
# NASequence: ribonucleotide pointers must come from RibonucleotideDB
# ---------------------------------------------------------------------------


def _db_ribo(code=b"A"):
    return pyopenms.RibonucleotideDB().getRibonucleotide(code)


def _foreign_ribo():
    return pyopenms.Ribonucleotide(
        "GrillNuc",
        "G!",
        "",
        "",
        pyopenms.EmpiricalFormula("C5H9O6P"),
        "G",
        0.0,
        0.0,
        pyopenms.Ribonucleotide.TermSpecificityNuc.ANYWHERE,
        pyopenms.EmpiricalFormula("C5H9O6P"),
    )


def test_nasequence_accepts_database_ribonucleotides():
    seq = pyopenms.NASequence.fromString("AUGC")
    seq.set(0, _db_ribo(b"G"))
    assert seq.toString().startswith("G")

    del seq
    gc.collect()


def test_nasequence_set_rejects_foreign_ribonucleotide():
    seq = pyopenms.NASequence.fromString("AUGC")
    with pytest.raises(ValueError, match="RibonucleotideDB"):
        seq.set(0, _foreign_ribo())
    # the sequence must be left untouched
    assert seq.toString() == "AUGC"


@pytest.mark.parametrize("setter", ["setFivePrimeMod", "setThreePrimeMod"])
def test_nasequence_terminal_mod_rejects_foreign_ribonucleotide(setter):
    seq = pyopenms.NASequence.fromString("AUGC")
    with pytest.raises(ValueError, match="RibonucleotideDB"):
        getattr(seq, setter)(_foreign_ribo())


def test_nasequence_setsequence_rejects_foreign_element():
    seq = pyopenms.NASequence.fromString("AUGC")
    good = _db_ribo(b"A")
    with pytest.raises(ValueError, match="RibonucleotideDB"):
        seq.setSequence([good, _foreign_ribo(), good])
    # a fully valid list still works
    seq.setSequence([good, good])
    assert seq.size() == 2


def test_nasequence_setsequence_is_atomic_on_rejection():
    """Every element is validated before any is stored, so a rejected list must leave the
    sequence exactly as it was."""
    seq = pyopenms.NASequence.fromString("AUGC")
    before = seq.toString()
    good = _db_ribo(b"A")
    with pytest.raises(ValueError):
        seq.setSequence([good, _foreign_ribo(), good])
    assert seq.toString() == before
    assert seq.size() == 4


@pytest.mark.parametrize("index", [4, 99, 2**31])
def test_nasequence_set_rejects_out_of_range_index(index):
    """``NASequence::set`` is an unchecked ``seq_[index] = r``, so an out-of-range index
    was an out-of-bounds write rather than an error."""
    seq = pyopenms.NASequence.fromString("AUGC")
    with pytest.raises(IndexError):
        seq.set(index, _db_ribo(b"G"))
    assert seq.toString() == "AUGC"


@pytest.mark.parametrize("index", [4, 99, 2**31])
def test_nasequence_get_rejects_out_of_range_index(index):
    """``NASequence::get`` is an unchecked ``seq_[index]``; ``__getitem__`` already guarded."""
    seq = pyopenms.NASequence.fromString("AUGC")
    with pytest.raises(IndexError):
        seq.get(index)


@pytest.mark.parametrize("index", [7, 99, 2**31])
def test_aasequence_setmodification_rejects_out_of_range_index(index):
    """``AASequence::setModification(Size, const Residue*)`` is an unchecked
    ``peptide_[index] = ...``; the ResidueModification overloads throw IndexOverflow."""
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    with pytest.raises(IndexError):
        seq.setModification(index, pyopenms.ResidueDB().getResidue("Alanine"))
    assert seq.toString() == "PEPTIDE"


@pytest.mark.parametrize(
    "call",
    [
        lambda: pyopenms.NASequence.fromString("AUGC").set(0, None),
        lambda: pyopenms.NASequence.fromString("AUGC").setFivePrimeMod(None),
        lambda: pyopenms.NASequence.fromString("AUGC").setThreePrimeMod(None),
        lambda: pyopenms.EnzymaticDigestion().setEnzyme(None),
        lambda: pyopenms.RNaseDigestion().setEnzyme(None),
        lambda: pyopenms.AASequence.fromString("PEPTIDE").setNTerminalModification(None),
        lambda: pyopenms.AASequence.fromString("PEPTIDE").setCTerminalModification(None),
        lambda: pyopenms.Residue().setModification(None),
    ],
    ids=[
        "NASequence.set", "setFivePrimeMod", "setThreePrimeMod", "EnzymaticDigestion.setEnzyme",
        "RNaseDigestion.setEnzyme", "setNTerminalModification", "setCTerminalModification",
        "Residue.setModification",
    ],
)
def test_none_is_rejected_not_dereferenced(call):
    """The guarded setters dereference their argument, so None must never reach them."""
    with pytest.raises((TypeError, ValueError)):
        call()


def test_nasequence_survives_deletion_of_the_python_wrapper():
    """A DB pointer stays valid after the Python wrapper that produced it is gone."""
    seq = pyopenms.NASequence.fromString("AUGC")
    r = _db_ribo(b"G")
    seq.set(0, r)
    before = seq.toString()

    del r
    _churn()

    assert seq.toString() == before


def test_copy_outlives_original_nasequence():
    """The property a keep_alive-based fix would not have provided."""
    seq = pyopenms.NASequence.fromString("AUGC")
    r = _db_ribo(b"G")
    seq.set(0, r)
    clone = pyopenms.NASequence(seq)
    expected = seq.toString()

    del seq, r
    _churn()

    assert clone.toString() == expected


# ---------------------------------------------------------------------------
# Enzymes: must come from ProteaseDB / RNaseDB
# ---------------------------------------------------------------------------


def _foreign_enzyme():
    return pyopenms.DigestionEnzyme("GrillZyme", "K", set(), "")


def test_enzymatic_digestion_accepts_database_enzyme():
    d = pyopenms.EnzymaticDigestion()
    d.setEnzyme(pyopenms.ProteaseDB().getEnzyme("Trypsin"))
    assert d.getEnzymeName() == "Trypsin"


def test_enzymatic_digestion_rejects_foreign_enzyme():
    d = pyopenms.EnzymaticDigestion()
    with pytest.raises(ValueError, match="ProteaseDB or RNaseDB"):
        d.setEnzyme(_foreign_enzyme())


def test_enzyme_survives_deletion_of_the_python_wrapper():
    d = pyopenms.EnzymaticDigestion()
    e = pyopenms.ProteaseDB().getEnzyme("Trypsin")
    d.setEnzyme(e)

    del e
    _churn()

    assert d.getEnzymeName() == "Trypsin"


def test_copy_outlives_original_digestion():
    """The exact scenario that defeats keep_alive: a copy with no lifetime edge."""
    d1 = pyopenms.EnzymaticDigestion()
    e = pyopenms.ProteaseDB().getEnzyme("Trypsin")
    d1.setEnzyme(e)
    d2 = pyopenms.EnzymaticDigestion(d1)

    del d1, e
    _churn()

    assert d2.getEnzymeName() == "Trypsin"


def test_repeated_setenzyme_does_not_accumulate_state():
    d = pyopenms.EnzymaticDigestion()
    db = pyopenms.ProteaseDB()
    for name in ("Trypsin", "Lys-C", "Arg-C", "Chymotrypsin", "Trypsin"):
        if db.hasEnzyme(name):
            d.setEnzyme(db.getEnzyme(name))
    assert d.getEnzymeName() in ("Trypsin", "Chymotrypsin", "Arg-C", "Lys-C")


def test_rnase_digestion_accepts_rnase():
    d = pyopenms.RNaseDigestion()
    d.setEnzyme(pyopenms.RNaseDB().getEnzyme("RNase_T1"))
    assert d.getEnzymeName() == "RNase_T1"


def test_rnase_digestion_rejects_protease():
    """RNaseDigestion::setEnzyme dynamic_casts to DigestionEnzymeRNA and dereferences
    the result unchecked, so a protease used to be an immediate null-pointer crash."""
    d = pyopenms.RNaseDigestion()
    with pytest.raises(ValueError, match="RNaseDB"):
        d.setEnzyme(pyopenms.ProteaseDB().getEnzyme("Trypsin"))


def test_rnase_digestion_rejects_foreign_enzyme():
    d = pyopenms.RNaseDigestion()
    with pytest.raises(ValueError, match="RNaseDB"):
        d.setEnzyme(_foreign_enzyme())


@pytest.mark.parametrize(
    "enzyme_factory",
    [
        lambda: pyopenms.ProteaseDB().getEnzyme("Trypsin"),
        _foreign_enzyme,
    ],
    ids=["protease", "foreign"],
)
def test_base_setenzyme_descriptor_guards_on_the_receiver_type(enzyme_factory):
    """``setEnzyme`` is virtual, so calling the *base* descriptor explicitly with a derived
    receiver dispatches to ``RNaseDigestion::setEnzyme`` -- whose unchecked downcast to
    ``DigestionEnzymeRNA`` terminates the process for a non-RNase.  The guard must key off
    the receiver's dynamic type, not off which binding was entered."""
    d = pyopenms.RNaseDigestion()
    with pytest.raises(ValueError, match="RNaseDB"):
        pyopenms.EnzymaticDigestion.setEnzyme(d, enzyme_factory())


def test_base_setenzyme_descriptor_still_works_for_the_base_class():
    d = pyopenms.EnzymaticDigestion()
    pyopenms.EnzymaticDigestion.setEnzyme(d, pyopenms.ProteaseDB().getEnzyme("Trypsin"))
    assert d.getEnzymeName() == "Trypsin"


# ---------------------------------------------------------------------------
# Modifications: interned into ModificationsDB rather than stored by reference
# ---------------------------------------------------------------------------


def _well_formed_mod(tag="GrillOx"):
    m = pyopenms.ResidueModification()
    m.setId(tag)
    m.setFullId("%s (N-term)" % tag)
    m.setFullName("Grill oxidation %s" % tag)
    m.setTermSpecificity(pyopenms.ResidueModification.TermSpecificity.N_TERM)
    m.setMonoMass(15.9949)
    m.setDiffMonoMass(15.9949)
    return m


@pytest.mark.parametrize("setter", ["setNTerminalModification", "setCTerminalModification"])
def test_terminal_modification_from_database_survives_deletion(setter):
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    mod = pyopenms.ModificationsDB().getModification("Oxidation")
    getattr(seq, setter)(mod)
    expected = seq.toString()

    del mod
    _churn()

    assert seq.toString() == expected


@pytest.mark.parametrize("setter", ["setNTerminalModification", "setCTerminalModification"])
def test_terminal_modification_is_interned_not_aliased(setter):
    """A Python-built modification must be copied into ModificationsDB, not stored raw."""
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    mod = _well_formed_mod("Grill" + setter[3:8])
    getattr(seq, setter)(mod)

    stored = (
        seq.getNTerminalModification()
        if setter.startswith("setN")
        else seq.getCTerminalModification()
    )
    assert stored is not mod, "modification was stored by reference instead of interned"

    expected = seq.toString()
    del mod
    _churn()
    assert seq.toString() == expected


def test_terminal_modification_accepts_both_keyword_spellings():
    """Both overloads were reachable before the fix; neither spelling may regress."""
    mod = pyopenms.ModificationsDB().getModification("Oxidation")

    positional = pyopenms.AASequence.fromString("PEPTIDE")
    positional.setNTerminalModification(mod)

    by_modification = pyopenms.AASequence.fromString("PEPTIDE")
    by_modification.setNTerminalModification(modification=mod)

    by_mod = pyopenms.AASequence.fromString("PEPTIDE")
    by_mod.setNTerminalModification(mod=mod)

    assert positional.toString() == by_modification.toString() == by_mod.toString()


def _anywhere_mod(tag):
    mod = _well_formed_mod(tag)
    mod.setTermSpecificity(pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    mod.setOrigin("P")
    mod.setFullId(tag)
    return mod


def test_indexed_modification_is_interned_into_modificationsdb():
    """The indexed setter funnels the pointer into ResidueDB, which caches the modified
    Residue for the life of the process.  A raw pointer stored there outlives every
    possible keep-alive nurse and corrupts *later*, unrelated sequences.

    Interning is directly observable: the by-reference overload registers the modification
    in ModificationsDB, whereas the raw-pointer path never did.
    """
    tag = "GrillIdxIntern"
    assert not pyopenms.ModificationsDB().has(tag)

    seq = pyopenms.AASequence.fromString("PEPTIDE")
    seq.setModification(0, _anywhere_mod(tag))

    assert pyopenms.ModificationsDB().has(tag), (
        "modification was stored by raw pointer instead of being interned into "
        "ModificationsDB; ResidueDB would cache a pointer to a Python-owned object"
    )


def test_indexed_modification_does_not_poison_residuedb():
    """A second, independently built sequence must not see the first one's modification
    after both the first sequence and the modification object are gone."""
    tag = "GrillIdxPoison"
    first = pyopenms.AASequence.fromString("PEPTIDE")
    first.setModification(0, _anywhere_mod(tag))
    expected = first.toString()

    del first
    _churn()

    second = pyopenms.AASequence.fromString("PEPTIDE")
    second.setModification(0, _anywhere_mod(tag))

    assert second.toString() == expected


def test_copy_outlives_original_aasequence():
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    mod = _well_formed_mod("GrillCopy")
    seq.setNTerminalModification(mod)
    clone = pyopenms.AASequence(seq)
    expected = seq.toString()

    del seq, mod
    _churn()

    assert clone.toString() == expected


def test_residue_setmodification_is_interned():
    res = pyopenms.Residue()
    mod = _well_formed_mod("GrillRes")
    mod.setTermSpecificity(pyopenms.ResidueModification.TermSpecificity.ANYWHERE)
    res.setModification(mod)
    expected = res.getModificationName()

    del mod
    _churn()

    assert res.getModificationName() == expected


def test_aasequence_setmodification_rejects_foreign_residue():
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    with pytest.raises(ValueError, match="ResidueDB"):
        seq.setModification(0, pyopenms.Residue())


def test_aasequence_setmodification_accepts_database_residue():
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    seq.setModification(0, pyopenms.ResidueDB().getResidue("Alanine"))
    assert seq.toString().startswith("A")


def test_residuedb_getmodifiedresidue_rejects_unregistered_modification():
    db = pyopenms.ResidueDB()
    with pytest.raises(ValueError, match="ModificationsDB"):
        db.getModifiedResidue(db.getResidue("Alanine"), _well_formed_mod("GrillUnreg"))


def test_residuedb_getmodifiedresidue_accepts_database_modification():
    db = pyopenms.ResidueDB()
    mod = pyopenms.ModificationsDB().getModification("Oxidation")
    res = db.getModifiedResidue(db.getResidue("Methionine"), mod)
    assert res is not None
    assert res.isModified()


# ---------------------------------------------------------------------------
# Isobaric wrappers store the quantitation method by reference
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        lambda qm: pyopenms.IsobaricQuantifier(qm),
        lambda qm: pyopenms.IsobaricNormalizer(qm),
        lambda qm: pyopenms.IsobaricChannelExtractor(qm),
    ],
    ids=["IsobaricQuantifier", "IsobaricNormalizer", "IsobaricChannelExtractor"],
)
def test_isobaric_wrapper_keeps_quant_method_alive(factory):
    qm = pyopenms.ItraqFourPlexQuantitationMethod()
    rc_before = sys.getrefcount(qm)
    obj = factory(qm)

    assert sys.getrefcount(qm) == rc_before + 1, (
        "wrapper stores the quantitation method by reference but does not keep it alive"
    )
    del qm
    _churn()
    assert obj is not None


@pytest.mark.parametrize(
    "cls",
    ["IsobaricQuantifier", "IsobaricNormalizer", "IsobaricChannelExtractor"],
    ids=["Quantifier", "Normalizer", "ChannelExtractor"],
)
def test_isobaric_wrapper_has_no_pointer_duplicating_copy(cls):
    """No reachable copy may duplicate the stored quantitation-method pointer.

    Such a copy would carry the raw pointer without the keep-alive edge and dangle once
    the original wrapper died.  Two outcomes are acceptable: copying raises, or the
    inherited ``DefaultParamHandler.__copy__`` slices to a plain ``DefaultParamHandler``
    that holds no quantitation method at all.  What must never happen is a copy that is
    still an instance of ``cls``.
    """
    qm = pyopenms.ItraqFourPlexQuantitationMethod()
    obj = getattr(pyopenms, cls)(qm)
    try:
        clone = copy.copy(obj)
    except TypeError:
        return
    assert not isinstance(clone, getattr(pyopenms, cls)), (
        "copying produced a %s that duplicates the quantitation-method pointer "
        "without a keep-alive edge" % cls
    )


# ---------------------------------------------------------------------------
# The KDTree bindings were unusable and have been removed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", ["KDTreeFeatureMaps", "KDTreeFeatureNode", "MapAlignmentAlgorithmKD"]
)
def test_unusable_kdtree_bindings_are_gone(name):
    """``KDTreeFeatureMaps.addMaps`` received the maps through a *copying* STL caster
    while C++ stored pointers into that call-local copy, so every accessor read freed
    memory -- silently wrong from ~100 features and a hard segfault from ~5000.  The
    other two classes were reachable only through it (``MapAlignmentAlgorithmKD``) or
    not constructible at all (``KDTreeFeatureNode``).  ``FeatureGroupingAlgorithmKD``
    is the supported entry point and builds the tree internally."""
    assert not hasattr(pyopenms, name)


def test_kd_feature_grouping_still_works():
    """The supported replacement for the removed low-level bindings."""
    assert hasattr(pyopenms, "FeatureGroupingAlgorithmKD")
    algo = pyopenms.FeatureGroupingAlgorithmKD()

    maps = []
    for m in range(2):
        fm = pyopenms.FeatureMap()
        for i in range(5):
            f = pyopenms.Feature()
            f.setRT(100.0 + 10 * i)
            f.setMZ(500.0 + i)
            f.setIntensity(1000.0)
            f.setCharge(1)
            fm.push_back(f)
        maps.append(fm)

    out = pyopenms.ConsensusMap()
    algo.group(maps, out)
    assert out.size() > 0
