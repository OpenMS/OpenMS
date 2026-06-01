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


if __name__ == "__main__":
    test_import_percolator()
    print("test_import_percolator: passed")
    test_low_level_rescore_on_separable_data()
    print("test_low_level_rescore_on_separable_data: passed")
    test_reproducible_with_same_seed()
    print("test_reproducible_with_same_seed: passed")
    print("\nAll Percolator pyOpenMS tests passed.")
