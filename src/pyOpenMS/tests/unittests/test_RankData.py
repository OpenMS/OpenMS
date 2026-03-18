import math
import numpy as np
import pyopenms as poms
import pytest

def test_tie_methods():
    a = np.array([0, 2, 3, 2], dtype=np.float64)
    m = poms.RankData.Method
    p = poms.RankData.NaNPolicy.Propagate

    np.testing.assert_allclose(poms.RankData.rankdata_double(a, m.Average, p), [1.0, 2.5, 4.0, 2.5])
    np.testing.assert_allclose(poms.RankData.rankdata_double(a, m.Min,     p), [1.0, 2.0, 4.0, 2.0])
    np.testing.assert_allclose(poms.RankData.rankdata_double(a, m.Max,     p), [1.0, 3.0, 4.0, 3.0])
    np.testing.assert_allclose(poms.RankData.rankdata_double(a, m.Dense,   p), [1.0, 2.0, 3.0, 2.0])
    np.testing.assert_allclose(poms.RankData.rankdata_double(a, m.Ordinal, p), [1.0, 2.0, 4.0, 3.0])

def test_nan_policies():
    a = np.array([1.0, np.nan, 1.0])
    m = poms.RankData.Method.Average

    # propagate -> all NaN
    out = poms.RankData.rankdata_double(a, m, poms.RankData.NaNPolicy.Propagate)
    assert out.shape == (3,)
    assert all(math.isnan(x) for x in out)

    # omit -> NaN stays at NaN position; others ranked among non-NaN
    out = poms.RankData.rankdata_double(a, m, poms.RankData.NaNPolicy.Omit)
    assert math.isclose(out[0], 1.5)
    assert math.isnan(out[1])
    assert math.isclose(out[2], 1.5)

    # raise -> ValueError from C++ (std::invalid_argument)
    with pytest.raises(ValueError):
        poms.RankData.rankdata_double(a, m, poms.RankData.NaNPolicy.Raise)

def test_empty_and_dtypes():
    out = poms.RankData.rankdata_double(np.array([], dtype=np.float64),
                                        poms.RankData.Method.Average,
                                        poms.RankData.NaNPolicy.Propagate)
    assert isinstance(out, np.ndarray)
    assert out.size == 0

    m = poms.RankData.Method.Min
    p = poms.RankData.NaNPolicy.Propagate

    np.testing.assert_allclose(
        poms.RankData.rankdata_double(np.array([0, 2, 2], dtype=np.float64), m, p),
        [1.0, 2.0, 2.0]
    )
    np.testing.assert_allclose(
        poms.RankData.rankdata_float(np.array([0, 2, 2], dtype=np.float32), m, p),
        [1.0, 2.0, 2.0]
    )
    np.testing.assert_allclose(
        poms.RankData.rankdata_int(np.array([0, 2, 2], dtype=np.int32), m, p),
        [1.0, 2.0, 2.0]
    )

def test_ordinal_stability():
    # two equal blocks; check stable order within each tie
    a = np.array([5, 5, 1, 5], dtype=float)
    out = poms.RankData.rankdata_double(a, poms.RankData.Method.Ordinal, poms.RankData.NaNPolicy.Propagate)
    # sorted stable indices by value: [2, 0, 1, 3] -> ranks: [2,3,1,4]
    np.testing.assert_allclose(out, [2.0, 3.0, 1.0, 4.0])
