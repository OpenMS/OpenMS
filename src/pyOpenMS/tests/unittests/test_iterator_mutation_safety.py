"""Mutation during iteration is defined, list-like behavior (issue #9792).

The __iter__ bindings are index-based value iterators: each step
bounds-checks against the container's *current* size and returns an owned
copy of the element. Mutating the container while iterating therefore
behaves like mutating a Python list during iteration -- appended elements
are visited, shrinking ends the loop early -- and can never dereference
invalidated C++ iterators. (The previous nb::make_iterator bindings held raw
begin()/end() iterators, so a push_back reallocation during iteration was a
use-after-free.) Param has no index access; its __iter__ snapshots owned
entry copies instead, so mutation during iteration changes nothing there.
"""
import pytest
import pyopenms

np = pytest.importorskip("numpy")


def _spectrum(rt):
    s = pyopenms.MSSpectrum()
    s.setRT(rt)
    s.set_peaks((np.linspace(100.0, 200.0, 50), np.ones(50, dtype=np.float32)))
    return s


def test_msexperiment_growth_during_iteration_is_visited_and_safe():
    exp = pyopenms.MSExperiment()
    exp.addSpectrum(_spectrum(1.0))
    exp.addSpectrum(_spectrum(2.0))

    visited = []
    grown = False
    for spec in exp:
        visited.append(spec.getRT())
        if not grown:
            grown = True
            # force at least one reallocation while the iterator is live
            for i in range(64):
                exp.addSpectrum(_spectrum(100.0 + i))

    assert visited[:2] == [1.0, 2.0]
    assert visited[2:] == [100.0 + i for i in range(64)]


def test_msspectrum_shrink_during_iteration_stops_early():
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 1.0]))

    visited = []
    for peak in spec:
        visited.append(peak.getMZ())
        spec.set_peaks(([9.0], [9.0]))  # shrink to one peak mid-iteration

    # step 0 sees the original first peak; afterwards the size is 1, so the
    # next bounds check ends the loop -- like iterating a shrinking list
    assert visited == [1.0]


def test_featuremap_growth_during_iteration():
    fm = pyopenms.FeatureMap()
    f = pyopenms.Feature()
    f.setRT(1.0)
    fm.push_back(f)

    visited = 0
    for feat in fm:
        visited += 1
        if visited == 1:
            f2 = pyopenms.Feature()
            f2.setRT(2.0)
            fm.push_back(f2)

    assert visited == 2


def test_iterated_elements_are_owned_copies():
    """Elements yielded during mutation are intact owned copies, never views
    of moved-from storage."""
    exp = pyopenms.MSExperiment()
    exp.addSpectrum(_spectrum(1.0))
    exp.addSpectrum(_spectrum(2.0))

    kept = []
    for spec in exp:
        kept.append(spec)
        exp.addSpectrum(_spectrum(50.0))
        if len(kept) == 10:
            break

    for spec in kept:
        mz, _ = spec.get_peaks()
        assert len(mz) == 50
        assert mz[0] == 100.0


def test_param_iteration_is_a_snapshot():
    p = pyopenms.Param()
    p.setValue("a", 1)
    p.setValue("b", 2)

    seen = 0
    for _entry in p:
        seen += 1
        p.setValue("c", 3)  # added mid-iteration: not visited (snapshot)

    assert seen == 2
    assert p.exists("c")


def test_aasequence_and_isotopedistribution_iteration_still_work():
    seq = pyopenms.AASequence.fromString("PEPTIDE")
    assert len(list(seq)) == 7

    iso = pyopenms.EmpiricalFormula("C100").getIsotopeDistribution(
        pyopenms.CoarseIsotopePatternGenerator(3))
    assert len(list(iso)) == 3
