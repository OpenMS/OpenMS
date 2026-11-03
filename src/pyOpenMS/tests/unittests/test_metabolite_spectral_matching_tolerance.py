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


@pytest.mark.parametrize("annotations", [False, True])
@pytest.mark.parametrize("tolerance", [0.0, -1.0, math.nan, math.inf])
def test_mvh_tolerance(annotations, tolerance):
    """Both MVH overloads reject invalid tolerances before scoring."""
    method = (oms.MetaboliteSpectralMatching.computeMVHScoreWithAnnotations
              if annotations else oms.MetaboliteSpectralMatching.computeMVHScore)
    with pytest.raises(RuntimeError, match="finite positive"):
        method(tolerance, False, spectrum(), spectrum())


def test_mvh_peak_credit_and_annotations():
    """Shared peaks retain both annotations and allow historical per-ion counting."""
    experimental = oms.MSSpectrum()
    experimental.set_peaks((np.arange(100.0, 107.0), np.array([10., 11., 12., 13., 14., 1000., 16.])))
    theoretical = oms.MSSpectrum()
    theoretical.set_peaks((np.array([104.98, 105.005]), np.ones(2)))
    names = oms.StringDataArray()
    names.setName("IonNames")
    names.push_back("c1")
    names.push_back("c2")
    charges = oms.IntegerDataArray()
    charges.setName("Charges")
    charges.push_back(-1)
    charges.push_back(-2)
    theoretical.setStringDataArrays([names])
    theoretical.setIntegerDataArrays([charges])
    score, annotations = oms.MetaboliteSpectralMatching.computeMVHScoreWithAnnotations(
        0.1, False, experimental, theoretical, 3, 1.0, 0.0, True, True)
    legacy = oms.MetaboliteSpectralMatching.computeMVHScore(
        0.1, False, experimental, theoretical, 3, 1.0, 0.0, True, False)
    assert score > 0
    assert legacy == 0
    assert len(annotations) == 2


def test_hyperscore_mass_weighting_is_explicit():
    """Existing callers keep unweighted scores and exact matching."""
    experimental = spectrum()
    theoretical = spectrum()
    theoretical.set_peaks((np.array([100.01, 200.01, 300.01]), np.full(3, 10.0)))
    plain = oms.MetaboliteSpectralMatching.computeHyperScore(0.1, False, experimental, theoretical)
    weighted = oms.MetaboliteSpectralMatching.computeHyperScore(0.1, False, experimental, theoretical, 0., True)
    assert weighted < plain
    with pytest.raises(RuntimeError, match="positive"):
        oms.MetaboliteSpectralMatching.computeHyperScore(0., False, experimental, theoretical, 0., True)
