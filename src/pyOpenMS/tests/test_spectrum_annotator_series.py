"""Ion-series evidence shared by peptide annotation and RNA graph scoring."""

import pyopenms


def test_longest_consecutive_ion_series():
    count = pyopenms.SpectrumAnnotator.compute_longest_ion_series
    assert count([]) == 0
    assert count([False, False]) == 0
    assert count([False, True, True, False, True]) == 2
    assert count([True, False, True, False, True]) == 1


def test_missing_rungs_are_connected_but_not_counted():
    count = pyopenms.SpectrumAnnotator.compute_longest_ion_series
    assert count([True, False, True, False, True], max_missing=1) == 3
    assert count([True, False, False, True, True], max_missing=1) == 2
    assert count([True, False, False, True, True], max_missing=2) == 3
    assert count([False, False, True, False, False], max_missing=2) == 1
