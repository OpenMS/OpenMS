"""
Deep-copy *independence* tests for pyopenms classes (OpenMS issue #9460, §2).

``test_copy_constructors.py`` already proves that ``copy.copy()``,
``copy.deepcopy()`` and ``ClassName(obj)`` return an object of the correct type
that compares equal to the (default-constructed) original. That sweep does
**not** catch the aliasing / object-slicing regression the nanobind rewrite
fixed: an old binding could hand back a "copy" that still shares the original's
underlying C++ storage, so mutating the copy would silently corrupt the
original. Comparing two freshly default-constructed empty objects for equality
passes regardless of whether they share state.

These tests pin the contract the issue asks for:

    copy -> mutate the copy -> assert the original is unchanged

For each class we build a *non-trivial* instance with a known observable value,
copy it three ways, mutate only the copy, and assert that

  * the copy starts out carrying the original's value,
  * mutating the copy leaves the ORIGINAL untouched (the real regression guard),
  * the mutation actually took effect on the copy.

The class list deliberately includes ``AASequence`` (the example called out in
the issue) and the ``MetaInfoInterface`` container types
(``MSExperiment``/``FeatureMap``/``ConsensusMap``), where shared internal
vectors/maps are exactly where an aliasing bug would hide.
"""

import copy

import pytest

import pyopenms


# --- builders: each returns a non-trivial instance with a known observable ---
def _b_peak1d():
    o = pyopenms.Peak1D(); o.setMZ(100.0); return o

def _b_peak2d():
    o = pyopenms.Peak2D(); o.setRT(10.0); return o

def _b_residue():
    o = pyopenms.Residue(); o.setName("Ala"); return o

def _b_peptide_hit():
    o = pyopenms.PeptideHit(); o.setScore(0.5); return o

def _b_peptide_id():
    o = pyopenms.PeptideIdentification(); o.setScoreType("q-value"); return o

def _b_protein_id():
    o = pyopenms.ProteinIdentification(); o.setIdentifier("run_1"); return o

def _b_spectrum():
    o = pyopenms.MSSpectrum(); o.setRT(5.0); return o

def _b_chromatogram():
    o = pyopenms.MSChromatogram(); o.setName("chrom_1"); return o

def _b_precursor():
    o = pyopenms.Precursor(); o.setMZ(500.0); return o

def _b_feature():
    o = pyopenms.Feature(); o.setRT(3.0); return o

def _b_consensus_feature():
    o = pyopenms.ConsensusFeature(); o.setCharge(2); return o

def _b_aasequence():
    return pyopenms.AASequence.fromString("PEPTIDE")

def _b_param():
    o = pyopenms.Param(); o.setValue("k", "x"); return o

def _b_datetime():
    o = pyopenms.DateTime(); o.set("2020-01-01 00:00:00"); return o

def _b_experiment():
    o = pyopenms.MSExperiment(); o.setMetaValue("k", "a"); return o

def _b_feature_map():
    o = pyopenms.FeatureMap(); o.setMetaValue("k", "a"); return o

def _b_consensus_map():
    o = pyopenms.ConsensusMap(); o.setMetaValue("k", "a"); return o


# (id, build, mutate, observe):
#   build()    -> non-trivial instance with a known observable value
#   mutate(o)  -> change that observable in place
#   observe(o) -> read the observable (hashable/comparable)
CASES = [
    ("Peak1D",                _b_peak1d,            lambda o: o.setMZ(200.0),                       lambda o: o.getMZ()),
    ("Peak2D",                _b_peak2d,            lambda o: o.setRT(20.0),                        lambda o: o.getRT()),
    ("Residue",               _b_residue,          lambda o: o.setName("Gly"),                     lambda o: o.getName()),
    ("PeptideHit",            _b_peptide_hit,       lambda o: o.setScore(0.9),                      lambda o: o.getScore()),
    ("PeptideIdentification", _b_peptide_id,        lambda o: o.setScoreType("FDR"),                lambda o: o.getScoreType()),
    ("ProteinIdentification", _b_protein_id,        lambda o: o.setIdentifier("run_2"),             lambda o: o.getIdentifier()),
    ("MSSpectrum",            _b_spectrum,          lambda o: o.setRT(15.0),                        lambda o: o.getRT()),
    ("MSChromatogram",        _b_chromatogram,      lambda o: o.setName("chrom_2"),                 lambda o: o.getName()),
    ("Precursor",             _b_precursor,         lambda o: o.setMZ(600.0),                       lambda o: o.getMZ()),
    ("Feature",               _b_feature,           lambda o: o.setRT(7.0),                         lambda o: o.getRT()),
    ("ConsensusFeature",      _b_consensus_feature, lambda o: o.setCharge(3),                       lambda o: o.getCharge()),
    ("AASequence",            _b_aasequence,        lambda o: o.setNTerminalModification("Acetyl"), lambda o: o.toString()),
    ("Param",                 _b_param,             lambda o: o.setValue("k", "y"),                 lambda o: str(o.getValue("k"))),
    ("DateTime",              _b_datetime,          lambda o: o.set("2021-06-15 12:00:00"),         lambda o: o.get()),
    ("MSExperiment",          _b_experiment,        lambda o: o.setMetaValue("k", "b"),             lambda o: str(o.getMetaValue("k"))),
    ("FeatureMap",            _b_feature_map,       lambda o: o.setMetaValue("k", "b"),             lambda o: str(o.getMetaValue("k"))),
    ("ConsensusMap",          _b_consensus_map,     lambda o: o.setMetaValue("k", "b"),             lambda o: str(o.getMetaValue("k"))),
]

CASE_IDS = [c[0] for c in CASES]


def _assert_independent(case, make_copy):
    name, build, mutate, observe = case

    original = build()
    baseline = observe(original)

    dup = make_copy(original)
    assert type(dup) is type(original), f"{name}: copy returned the wrong type"
    assert observe(dup) == baseline, f"{name}: copy did not preserve the original's value"

    mutate(dup)

    # The actual regression guard: a shared/aliased copy would change the original here.
    assert observe(original) == baseline, (
        f"{name}: mutating the copy changed the ORIGINAL "
        f"(copy-constructor aliasing/slicing regression)"
    )
    assert observe(dup) != baseline, f"{name}: mutation did not take effect on the copy"


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_copy_copy_is_independent(case):
    _assert_independent(case, copy.copy)


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_deepcopy_is_independent(case):
    _assert_independent(case, copy.deepcopy)


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_copy_constructor_is_independent(case):
    _assert_independent(case, lambda o: type(o)(o))
