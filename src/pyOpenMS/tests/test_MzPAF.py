"""Regression tests for mzPAF satellite-ion annotations (OpenMS #10175)."""

import copy

import pytest

import pyopenms as p


@pytest.mark.parametrize(
    "text,series,ordinal,subtype",
    [
        ("d5", p.MzPAFIonSeries.D, 5, None),
        ("v7", p.MzPAFIonSeries.V, 7, None),
        ("w3", p.MzPAFIonSeries.W, 3, None),
        ("da12", p.MzPAFIonSeries.D, 12, "a"),
        ("db4", p.MzPAFIonSeries.D, 4, "b"),
        ("wa12", p.MzPAFIonSeries.W, 12, "a"),
        ("wb4", p.MzPAFIonSeries.W, 4, "b"),
    ],
)
def test_satellite_roundtrip(text, series, ordinal, subtype):
    ann = p.MzPAF.parse(text)
    assert ann.ion_series == series
    assert ann.ordinal == ordinal
    assert ann.satellite_subtype == subtype
    assert ann.isValid()
    assert p.MzPAF.isMzPAFFormat(text)
    assert p.MzPAF.isStandardFragmentIon(series)
    assert p.MzPAF.ionSeriesToChar(series) == text[0]
    assert p.MzPAF.charToIonSeries(text[0]) == series
    assert p.MzPAF.toString(ann) == text
    assert p.MzPAF.parse(p.MzPAF.toString(ann)) == ann
    assert copy.copy(ann) == ann
    assert copy.deepcopy(ann) == ann

    constructed = p.MzPAFAnnotation()
    assert constructed.satellite_subtype is None
    constructed.ion_series = series
    constructed.ordinal = ordinal
    constructed.satellite_subtype = subtype
    assert constructed == ann
    assert p.MzPAF.toString(constructed) == text


@pytest.mark.parametrize("prefix", ["d", "v", "w", "da", "db", "wa", "wb"])
def test_satellite_modifiers_and_peak_annotation(prefix):
    ann = p.MzPAF.parse(f"1@{prefix}3{{LIR}}-H2O+2i^2/-1.4ppm*0.75")
    assert ann.analyte_index == 1
    assert ann.embedded_sequence == "LIR"
    assert len(ann.neutral_losses) == 1
    assert ann.isotope_offset == 2
    assert ann.charge == 2
    assert ann.mass_delta.value == pytest.approx(-1.4)
    assert ann.mass_delta.unit == p.MzPAFDeltaUnit.PPM
    assert ann.confidence == pytest.approx(0.75)
    assert p.MzPAF.parse(p.MzPAF.toString(ann)) == ann

    peak = p.MzPAF.toPeakAnnotation(ann, 500.123, 1000.0)
    assert peak.charge == 2
    assert peak.mz == pytest.approx(500.123)
    assert peak.intensity == pytest.approx(1000.0)
    restored = p.MzPAF.fromPeakAnnotation(peak)
    assert restored.size() == 1
    assert restored.annotations[0] == ann


def test_multiple_satellite_annotations():
    text = "d5,da5,db5,v7,w3,wa3,wb3,y4^2"
    anns = p.MzPAF.parseMultiple(text)
    assert anns.size() == 8
    assert p.MzPAF.toStringMultiple(anns) == text
    assert p.MzPAF.parseMultiple(p.MzPAF.toStringMultiple(anns)) == anns


def test_satellite_subtype_validation_and_equality():
    ann = p.MzPAF.parse("wa3")
    ann.satellite_subtype = "b"
    assert ann.isValid()
    assert ann == p.MzPAF.parse("wb3")
    assert not ann == p.MzPAF.parse("wa3")
    ann.satellite_subtype = None
    assert ann == p.MzPAF.parse("w3")

    ann.ordinal = None
    assert not ann.isValid()
    ann.ordinal = 3
    ann.satellite_subtype = "c"
    assert not ann.isValid()
    # toString() stays total: the non-conformant subtype is dropped, not raised on.
    assert p.MzPAF.toString(ann) == "w3"

    ann.satellite_subtype = "a"
    ann.ion_series = p.MzPAFIonSeries.V
    assert not ann.isValid()
    assert p.MzPAF.toString(ann) == "v3"


@pytest.mark.parametrize(
    "text",
    [
        "d", "v", "w", "da", "db", "wa", "wb", "va3", "vb3",
        "aa3", "ba3", "ca3", "xa3", "ya3", "za3", "dc3", "wc3",
        "daa3", "wba3", "wA3", "d3a", "da^2", "wb-H2O",
        "da99999999999999999999", "va 3", "dc 3", "wafoo 3",
    ],
)
def test_reject_malformed_satellite_annotations(text):
    assert p.MzPAF.tryParse(text) is None
    assert not p.MzPAF.isMzPAFFormat(text)
    assert p.MzPAF.tryParseMultiple(f"y4,{text}") is None
    with pytest.raises(RuntimeError):
        p.MzPAF.parse(text)
