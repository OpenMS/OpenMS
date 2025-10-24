import unittest
import numpy as np
from pyopenms import RankData

def test_rankdata_statics_match_docs():
    a = [0.0, 2.0, 3.0, 2.0]
    r_avg = RankData.rankdata_double(a, RankData.Average, RankData.Propagate)
    r_min = RankData.rankdata_double(a, RankData.Min, RankData.Propagate)
    r_max = RankData.rankdata_double(a, RankData.Max, RankData.Propagate)
    r_den = RankData.rankdata_double(a, RankData.Dense, RankData.Propagate)
    r_ord = RankData.rankdata_double(a, RankData.Ordinal, RankData.Propagate)
    assert list(r_avg) == [1.0, 2.5, 4.0, 2.5]
    assert list(r_min) == [1.0, 2.0, 4.0, 2.0]
    assert list(r_max) == [1.0, 3.0, 4.0, 3.0]
    assert list(r_den) == [1.0, 2.0, 3.0, 2.0]
    assert list(r_ord) == [1.0, 2.0, 4.0, 3.0]

def test_rankdata_addon_convenience():
    import pyopenms
    a = [0, 2, 3, 2]
    r = pyopenms.rankdata(a, method="average", nan_policy="propagate")
    assert np.allclose(r, [1.0, 2.5, 4.0, 2.5])
