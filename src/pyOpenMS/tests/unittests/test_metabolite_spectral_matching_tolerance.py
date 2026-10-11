"""HyperScore validates tolerances through both existing Python entry points."""

import math

import numpy as np
import pyopenms as oms
import pytest


def spectrum():
    result = oms.MSSpectrum()
    result.set_peaks((np.array([100.0, 200.0, 300.0]), np.array([10.0, 10.0, 10.0])))
    return result


@pytest.mark.parametrize("tolerance", [-1.0, math.nan, math.inf, -math.inf])
@pytest.mark.parametrize("ppm", [False, True])
@pytest.mark.parametrize("empty", [False, True])
@pytest.mark.parametrize("annotations", [False, True])
def test_invalid_hyperscore_tolerance(tolerance, ppm, empty, annotations):
    experimental = oms.MSSpectrum() if empty else spectrum()
    method = (
        oms.MetaboliteSpectralMatching.computeHyperScoreWithAnnotations
        if annotations
        else oms.MetaboliteSpectralMatching.computeHyperScore
    )
    with pytest.raises(RuntimeError, match="finite nonnegative fragment tolerance"):
        method(tolerance, ppm, experimental, spectrum(), 0.0)


@pytest.mark.parametrize("tolerance", [0.0, 0.5])
@pytest.mark.parametrize("ppm", [False, True])
def test_valid_hyperscore_tolerance(tolerance, ppm):
    experimental = spectrum()
    expected = math.log(300.0) + math.log(6.0)
    score = oms.MetaboliteSpectralMatching.computeHyperScore(
        tolerance, ppm, experimental, experimental, 0.0
    )
    annotated_score, _ = oms.MetaboliteSpectralMatching.computeHyperScoreWithAnnotations(
        tolerance, ppm, experimental, experimental, 0.0
    )
    assert score == pytest.approx(expected)
    assert annotated_score == pytest.approx(expected)
    assert oms.MetaboliteSpectralMatching.computeHyperScore(
        tolerance, ppm, oms.MSSpectrum(), experimental, 0.0
    ) == 0.0
