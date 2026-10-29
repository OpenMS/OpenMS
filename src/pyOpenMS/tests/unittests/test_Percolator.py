"""Smoke tests for OpenMS.Percolator bindings.

Run directly:
    PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/unittests/test_Percolator.py

Note: running this via pytest from the repo root may fail due to path-shadowing
between the source pyopenms/ (no compiled extensions) and the build pyopenms/.
Run from /tmp OR invoke via the command above.
"""

import pyopenms as oms


def test_import_percolator():
    assert hasattr(oms, "Percolator")
    assert hasattr(oms, "RescoreInput")
    assert hasattr(oms, "RescoreOutput")


def test_low_level_rescore_on_separable_data():
    """Low-level rescore on a trivially separable feature matrix."""
    ri = oms.RescoreInput()
    ri.feature_names = ["separator", "noise"]
    features = []
    labels = []
    for i in range(200):
        features.append([+1.0 + 0.01 * i, float(i % 7)])
        labels.append(False)
    for i in range(200):
        features.append([-1.0 + 0.01 * i, float(i % 11)])
        labels.append(True)
    ri.features = features
    ri.is_decoy = labels

    p = oms.Percolator()
    par = p.getDefaults()
    par.setValue("seed", 42)
    p.setParameters(par)

    out = p.rescore(ri)
    assert len(out.scores) == len(features)

    target_scores = sorted(s for s, d in zip(out.scores, labels) if not d)
    decoy_scores = sorted(s for s, d in zip(out.scores, labels) if d)
    assert target_scores[len(target_scores) // 2] > decoy_scores[len(decoy_scores) // 2]


def test_reproducible_with_same_seed():
    ri = oms.RescoreInput()
    ri.feature_names = ["f"]
    features = []
    labels = []
    for i in range(400):
        is_decoy = (i % 2 == 1)
        features.append([(-1.0 if is_decoy else +1.0) + 0.01 * i])
        labels.append(is_decoy)
    ri.features = features
    ri.is_decoy = labels

    p1, p2 = oms.Percolator(), oms.Percolator()
    for p in (p1, p2):
        par = p.getDefaults()
        par.setValue("seed", 17)
        p.setParameters(par)

    o1 = p1.rescore(ri)
    o2 = p2.rescore(ri)
    assert len(o1.scores) == len(o2.scores)
    for a, b in zip(o1.scores, o2.scores):
        assert abs(a - b) < 1e-9



def test_dense_feature_matrix_and_peptide_level():
    """The dense row-major feature_matrix gives the scores of the same rows as features; with peptide_level and
    peptide_keys the peptide-level q-values and PEPs are filled, one per row."""
    rows, labels, keys = [], [], []
    for i in range(400):
        is_decoy = (i % 2 == 1)
        rows.append([(-1.0 if is_decoy else +1.0) + 0.01 * i, float(i % 5)])
        labels.append(is_decoy)
        keys.append("PEP" + str(i // 2))  # two PSMs per peptide and label

    def make_percolator(peptide_level):
        p = oms.Percolator()
        par = p.getDefaults()
        par.setValue("seed", 23)
        par.setValue("peptide_level", "true" if peptide_level else "false")
        p.setParameters(par)
        return p

    nested = oms.RescoreInput()
    nested.feature_names = ["f", "g"]
    nested.features = rows
    nested.is_decoy = labels
    dense = oms.RescoreInput()
    dense.feature_names = ["f", "g"]
    dense.feature_matrix = [v for row in rows for v in row]
    dense.is_decoy = labels
    assert len(dense.feature_matrix) == 2 * len(rows)

    o_nested = make_percolator(False).rescore(nested)
    o_dense = make_percolator(False).rescore(dense)
    assert list(o_dense.scores) == list(o_nested.scores)
    assert len(o_dense.peptide_q_values) == 0

    dense.peptide_keys = keys
    assert list(dense.peptide_keys) == keys
    o_peptide = make_percolator(True).rescore(dense)
    assert len(o_peptide.peptide_q_values) == len(rows)
    assert len(o_peptide.peptide_peps) == len(rows)
    assert all(0.0 <= q <= 1.0 for q in o_peptide.peptide_q_values)

def _synthetic_psms(n=600):
    import random
    rng = random.Random(7)
    ids = oms.PeptideIdentificationList()
    aas = "ACDEFGHILMNPQRSTVWY"
    for i in range(n):
        decoy = (i % 2 == 1)
        seq = "".join(aas[(i * 7 + j * 3) % len(aas)] for j in range(8)) + "K"
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setCharge(2)
        hit.setScore((0.0 if decoy else 2.0) + rng.gauss(0.0, 1.0))
        hit.setMetaValue("target_decoy", "decoy" if decoy else "target")
        ev = oms.PeptideEvidence()
        ev.setProteinAccession(("DECOY_P" if decoy else "P") + str(i))
        hit.setPeptideEvidences([ev])
        pid = oms.PeptideIdentification()
        pid.setScoreType("hyperscore")
        pid.setHigherScoreBetter(True)
        pid.setMZ(hit.getSequence().getMZ(2))
        pid.setRT(10.0 * i)
        pid.setMetaValue("spectrum_reference", "scan=" + str(i + 1))
        pid.setHits([hit])
        ids.append(pid)
    return ids


def test_rescore_psms():
    ids = _synthetic_psms()
    p = oms.Percolator()
    features = list(oms.PercolatorInfile.getStandardFeatureSet(2, 2)) + ["score", "Peptide", "Proteins"]
    n = p.rescorePSMs(ids, features, "trypsin", 2, 2)
    assert n == ids.size()
    pid = ids[0]
    assert pid.getScoreType() == "q-value"
    assert not pid.isHigherScoreBetter()
    hit = pid.getHits()[0]
    assert hit.metaValueExists("MS:1001491")
    assert hit.metaValueExists("MS:1001493")
    assert hit.metaValueExists("hyperscore")


if __name__ == "__main__":
    test_import_percolator()
    print("test_import_percolator: passed")
    test_low_level_rescore_on_separable_data()
    print("test_low_level_rescore_on_separable_data: passed")
    test_reproducible_with_same_seed()
    print("test_reproducible_with_same_seed: passed")
    test_rescore_psms()
    print("test_rescore_psms: passed")
    print("\nAll Percolator pyOpenMS tests passed.")
