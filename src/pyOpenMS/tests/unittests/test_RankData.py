import math
import numpy as np
import pyopenms as poms
import pytest

def test_tie_methods():
    a = [0, 2, 3, 2]

    np.testing.assert_allclose(poms.rankdata(a, "average"), [1.0, 2.5, 4.0, 2.5])
    np.testing.assert_allclose(poms.rankdata(a, "min"),     [1.0, 2.0, 4.0, 2.0])
    np.testing.assert_allclose(poms.rankdata(a, "max"),     [1.0, 3.0, 4.0, 3.0])
    np.testing.assert_allclose(poms.rankdata(a, "dense"),   [1.0, 2.0, 3.0, 2.0])
    np.testing.assert_allclose(poms.rankdata(a, "ordinal"), [1.0, 2.0, 4.0, 3.0])

def test_nan_policies():
    a = np.array([1.0, np.nan, 1.0])

    # propagate -> all NaN
    out = poms.rankdata(a, "average", "propagate")
    assert out.shape == (3,)
    assert all(math.isnan(x) for x in out)

    # omit -> NaN stays at NaN position; others ranked among non-NaN
    out = poms.rankdata(a, "average", "omit")
    assert math.isclose(out[0], 1.5)
    assert math.isnan(out[1])
    assert math.isclose(out[2], 1.5)

    # raise -> ValueError from the wrapper
    with pytest.raises(ValueError):
        poms.rankdata(a, "average", "raise")

def test_empty_and_dtypes():
    out = poms.rankdata([], "average")
    assert isinstance(out, np.ndarray)
    assert out.size == 0

    # dtype routing
    np.testing.assert_allclose(
        poms.rankdata(np.array([0, 2, 2], dtype=np.float64), "min", dtype="float64"),
        [1.0, 2.0, 2.0]
    )
    np.testing.assert_allclose(
        poms.rankdata(np.array([0, 2, 2], dtype=np.float32), "min", dtype="float32"),
        [1.0, 2.0, 2.0]
    )
    np.testing.assert_allclose(
        poms.rankdata(np.array([0, 2, 2], dtype=np.int32), "min", dtype="int32"),
        [1.0, 2.0, 2.0]
    )

def test_ordinal_stability():
    # two equal blocks; check stable order within each tie
    a = np.array([5, 5, 1, 5], dtype=float)
    out = poms.rankdata(a, "ordinal")
    # sorted stable indices by value: [2, 0, 1, 3] -> ranks: [2,3,1,4]
    np.testing.assert_allclose(out, [2.0, 3.0, 1.0, 4.0])
