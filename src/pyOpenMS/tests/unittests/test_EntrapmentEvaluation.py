"""Tests for the EntrapmentEvaluation bindings (entrapment FDP estimation).

Run directly:
    PYTHONPATH=OpenMS-build/pyOpenMS python3 -m pytest src/pyOpenMS/tests/unittests/test_EntrapmentEvaluation.py

Run from /tmp to avoid path shadowing between the source and the build pyopenms/.
"""

import math

import pyopenms as oms

EE = oms.EntrapmentEvaluation


def _entry(identifier, sequence):
    e = oms.FASTAEntry()
    e.identifier = identifier
    e.sequence = sequence
    return e


def _hit(sequence, accessions, decoy=False, cls=None):
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    evidences = []
    for acc in accessions:
        ev = oms.PeptideEvidence()
        ev.setProteinAccession(acc)
        evidences.append(ev)
    hit.setPeptideEvidences(evidences)
    hit.setMetaValue("target_decoy", "decoy" if decoy else "target")
    if cls is not None:
        hit.setMetaValue("mod_class", cls)
    return hit


def _id(sequence, q, accessions, decoy=False, cls=None):
    hit = _hit(sequence, accessions, decoy, cls)
    hit.setScore(q)
    pep_id = oms.PeptideIdentification()
    pep_id.setScoreType("q-value")
    pep_id.setHigherScoreBetter(False)
    pep_id.setHits([hit])
    return pep_id


def test_classify():
    assert EE.classify(_hit("PEPTIDEK", ["P1"]), "ENTRAP_") == EE.Kind.TARGET
    assert EE.classify(_hit("PEPTIDEK", ["ENTRAP_P1"]), "ENTRAP_") == EE.Kind.ENTRAPMENT
    # shared with a target protein: target
    assert EE.classify(_hit("PEPTIDEK", ["ENTRAP_P1", "P2"]), "ENTRAP_") == EE.Kind.TARGET
    # a decoy accession is ignored when the decoy string is given
    hit = _hit("PEPTIDEK", ["ENTRAP_P1", "DECOY_P2"])
    assert EE.classify(hit, "ENTRAP_") == EE.Kind.TARGET
    assert EE.classify(hit, "ENTRAP_", decoy_string="DECOY_") == EE.Kind.ENTRAPMENT
    assert EE.classify(_hit("PEPTIDEK", ["DECOY_P1"], decoy=True), "ENTRAP_") == EE.Kind.DECOY


def test_te_ratio_target_peptides_and_pairs():
    db = [
        _entry("P1", "AAAAAAAKCCCCICCKDDDDDDDK"),
        _entry("ENTRAP_P1", "AAAAAAAKCCCCLCCKEEEEEEEK"),  # AAAAAAAK shared, CCCCLCCK = CCCCICCK with I = L
        _entry("DECOY_P1", "GGGGGGGKHHHHHHHKIIIIIIIK"),
    ]
    db_pairs = db + [_entry("P2", "ACDEFGHK"), _entry("ENTRAP_P2", "HGFEDCAK")]
    # targets: 3 peptides; entrapment-only: EEEEEEEK
    assert math.isclose(EE.estimateTERatio(db, "ENTRAP_", "DECOY_", True, "Trypsin", 0, 7, 40), 3.0)

    targets = EE.targetPeptides(db, "ENTRAP_", "DECOY_", True, "Trypsin", 0, 7, 40)
    assert isinstance(targets, set)
    assert targets == {"AAAAAAAK", "CCCCLCCK", "DDDDDDDK"}  # I replaced by L

    pairs, mismatched = EE.pairPeptides(db_pairs, "ENTRAP_", "DECOY_", True, "Trypsin", 0, 7, 40)
    assert isinstance(pairs, dict)
    # HGFEDCAK is a permutation of ACDEFGHK: paired
    assert pairs == {"HGFEDCAK": "ACDEFGHK"}
    # AAAAAAAK equals its target: not paired. CCCCLCCK vs CCCCICCK and EEEEEEEK vs DDDDDDDK differ in composition
    # (compared before I/L folding): skipped and counted
    assert mismatched == 2


def test_qvalue():
    pep_id = _id("PEPTIDEK", 0.02, ["P1"])
    assert math.isclose(EE.qValue(pep_id, pep_id.getHits()[0]), 0.02)
    pep_id.setScoreType("hyperscore")
    try:
        EE.qValue(pep_id, pep_id.getHits()[0])
        assert False, "expected an exception without a q-value"
    except Exception:
        pass


def _ids():
    ids = oms.PeptideIdentificationList()
    for seq in ["AAAK", "CCCK", "DDDK", "EEEK", "FFFK", "GGGK", "HHHK"]:
        ids.append(_id(seq, 0.001, ["P1"], cls="unmod"))
    ids.append(_id("AAAK", 0.002, ["P1"], cls="unmod"))
    ids.append(_id("KLLLK", 0.005, ["ENTRAP_P9"], cls="mod"))
    ids.append(_id("KMMMK", 0.009, ["ENTRAP_P9"], cls="mod"))
    ids.append(_id("NNNK", 0.001, ["DECOY_P1"], decoy=True))
    ids.append(_id("PPPK", 0.5, ["ENTRAP_P2"]))
    return ids


def test_evaluate():
    res = EE.evaluate(_ids(), 0.01, "ENTRAP_", 1.0, class_meta_value="mod_class")
    assert isinstance(res, dict)
    overall = res["all"]
    assert isinstance(overall, EE.Result)
    assert overall.accepted == 10
    assert overall.entrapment == 2
    assert math.isclose(overall.lowerBoundFDP(), 0.2)
    assert math.isclose(overall.combinedFDP(), 0.4)
    assert not overall.paired
    assert math.isnan(overall.pairedFDP())
    assert res["mod"].accepted == 2
    assert res["unmod"].entrapment == 0

    # peptide level: 9 distinct peptides
    res = EE.evaluate(_ids(), 0.01, "ENTRAP_", 3.0, level=EE.Level.PEPTIDE)
    assert res["all"].accepted == 9
    assert math.isclose(res["all"].combinedFDP(), 2.0 / 9.0 * 4.0)


def test_evaluate_target_peptides_and_pairs():
    ids = _ids()
    ids.append(_id("KILLK", 0.003, ["ENTRAP_P8"]))
    # KILLK equals the target peptide KLLLK with I = L: counted as target
    res = EE.evaluate(ids, 0.01, "ENTRAP_", 1.0, target_peptides={"KLLLK", "AAAK"})
    assert res["all"].accepted == 11
    assert res["all"].entrapment == 1
    assert EE.evaluate(ids, 0.01, "ENTRAP_", 1.0, target_peptides=None)["all"].entrapment == 3

    # paired counts at peptide level: KMMMK's target MKMMK is not accepted
    res = EE.evaluate(ids, 0.01, "ENTRAP_", 1.0, level=EE.Level.PEPTIDE,
                      target_peptides={"KLLLK"}, pairs={"KMMMK": "MKMMK"})
    overall = res["all"]
    assert overall.paired
    assert overall.entrapment == 1
    assert overall.entrapment_target_not_accepted == 1
    assert overall.entrapment_above_target == 0
    assert math.isclose(overall.pairedFDP(), (1 + 1) / overall.accepted)


def test_results_are_independent():
    ids = _ids()
    first = EE.evaluate(ids, 0.01, "ENTRAP_", 1.0)["all"]
    first.accepted = 0
    # a new evaluation does not see the change, and a copy is independent of its source
    assert EE.evaluate(ids, 0.01, "ENTRAP_", 1.0)["all"].accepted == 10
    copied = EE.Result(first)
    copied.entrapment = 99
    assert first.entrapment == 2
