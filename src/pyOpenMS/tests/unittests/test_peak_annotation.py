"""Tests for theoretical fragment m/z values on peak annotations."""

import pytest

import pyopenms as oms


def test_peak_annotation_mz_errors_are_optional_and_signed():
    annotation = oms.PeptideHit_PeakAnnotation()
    annotation.mz = 499.999

    assert annotation.theoretical_mz is None
    assert annotation.getMZError() is None
    assert annotation.getMZErrorPPM() is None

    annotation.theoretical_mz = 500.0

    assert annotation.theoretical_mz == 500.0
    assert annotation.getMZError() == pytest.approx(-0.001)
    assert annotation.getMZErrorPPM() == pytest.approx(-2.0)

    annotation.theoretical_mz = None

    assert annotation.theoretical_mz is None
    assert annotation.getMZError() is None
    assert annotation.getMZErrorPPM() is None
