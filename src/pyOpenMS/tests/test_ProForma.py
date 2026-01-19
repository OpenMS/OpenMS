"""
Tests for ProForma v2 parser/writer in pyOpenMS.

These tests cover:
- Basic parsing and roundtrip serialization
- WriteMode (LOSSLESS vs CANONICAL)
- AASequence conversion
- Error handling
"""

import pytest


def test_parse_simple_sequence():
    """Test parsing simple peptide sequences."""
    import pyopenms as p

    # Simple sequence without modifications
    pf = p.parse("PEPTIDE")
    assert pf is not None

    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "PEPTIDE"


def test_parse_with_unimod_modification():
    """Test parsing sequences with UNIMOD accessions."""
    import pyopenms as p

    # Oxidation on M
    pf = p.parse("EM[UNIMOD:35]K")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "EM[UNIMOD:35]K"

    # Phosphorylation on S
    pf = p.parse("PES[UNIMOD:21]TIDE")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "PES[UNIMOD:21]TIDE"


def test_parse_with_named_modification():
    """Test parsing sequences with named modifications."""
    import pyopenms as p

    pf = p.parse("EM[Oxidation]K")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "EM[Oxidation]K"


def test_parse_with_mass_delta():
    """Test parsing sequences with mass delta modifications."""
    import pyopenms as p

    pf = p.parse("EM[+15.9949]K")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    # In lossless mode, the original formatting is preserved
    assert "+15.9949" in result


def test_roundtrip_lossless():
    """Test that lossless mode preserves original formatting."""
    import pyopenms as p

    test_cases = [
        "PEPTIDE",
        "EM[UNIMOD:35]K",
        "EM[Oxidation]K",
        "[Acetyl]-PEPTIDE",
        "PEPTIDE-[Amidated]",
    ]

    for original in test_cases:
        pf = p.parse(original)
        result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
        assert result == original, f"Roundtrip failed for {original}: got {result}"


def test_write_mode_canonical():
    """Test that canonical mode normalizes output."""
    import pyopenms as p

    # Parse with limited precision
    pf = p.parse("EM[+15.99]K")

    # Canonical mode should output 4 decimal places
    canonical = p.toString(pf, p.ProFormaWriteMode.CANONICAL)
    assert "+15.9900" in canonical

    # Lossless mode should preserve original
    lossless = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert "+15.99" in lossless


def test_terminal_modifications():
    """Test parsing N-terminal and C-terminal modifications."""
    import pyopenms as p

    # N-terminal modification
    pf = p.parse("[Acetyl]-PEPTIDE")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "[Acetyl]-PEPTIDE"

    # C-terminal modification
    pf = p.parse("PEPTIDE-[Amidated]")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "PEPTIDE-[Amidated]"

    # Both terminal modifications
    pf = p.parse("[Acetyl]-PEPTIDE-[Amidated]")
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "[Acetyl]-PEPTIDE-[Amidated]"


def test_parse_ion_with_charge():
    """Test parsing peptidoform ions with charge states."""
    import pyopenms as p

    pfi = p.parseIon("PEPTIDE/2")
    # Verify parsing succeeded and chain was extracted
    assert len(pfi.chains) == 1
    # The charge should be accessible through the charge attribute
    assert pfi.charge == 2


def test_parse_ion_multiple_chains():
    """Test parsing multi-chain peptidoforms."""
    import pyopenms as p

    pfi = p.parseIon("PEPTIDE//SEQUENCE")
    # Verify multiple chains were parsed
    assert len(pfi.chains) == 2


def test_is_representable_as_aasequence():
    """Test checking if a peptidoform can be converted to AASequence."""
    import pyopenms as p

    # Simple sequence should be representable
    pf = p.parse("PEPTIDE")
    assert p.isRepresentableAsAASequence(pf) == True

    # Unlocalised modification is not representable
    pf = p.parse("[Phospho]?PEPTIDE")
    assert p.isRepresentableAsAASequence(pf) == False


def test_to_aasequence_simple():
    """Test converting simple peptidoform to AASequence."""
    import pyopenms as p

    pf = p.parse("PEPTIDE")
    aas = p.toAASequence(pf, p.AASequenceConversionPolicy.STRICT)
    assert aas.toString() == "PEPTIDE"


def test_to_aasequence_with_modification():
    """Test converting peptidoform with modification to AASequence."""
    import pyopenms as p

    pf = p.parse("EM[UNIMOD:35]K")
    p.resolveModifications(pf)
    aas = p.toAASequence(pf, p.AASequenceConversionPolicy.STRICT)
    # The AASequence should contain the modification
    assert "Oxidation" in aas.toString() or "(Oxidation)" in aas.toString()


def test_to_aasequence_strict_fails_on_unsupported():
    """Test that STRICT policy raises error for unsupported features."""
    import pyopenms as p

    # Unlocalised modification cannot be converted
    pf = p.parse("[Phospho]?PEPTIDE")

    with pytest.raises(Exception):
        p.toAASequence(pf, p.AASequenceConversionPolicy.STRICT)


def test_to_aasequence_best_effort():
    """Test that BEST_EFFORT policy skips unsupported features."""
    import pyopenms as p

    # Unlocalised modification should be skipped
    pf = p.parse("[Phospho]?PEPTIDE")

    # Should not raise
    aas = p.toAASequence(pf, p.AASequenceConversionPolicy.BEST_EFFORT)
    assert aas.toString() == "PEPTIDE"


def test_from_aasequence():
    """Test creating peptidoform from AASequence."""
    import pyopenms as p

    aas = p.AASequence.fromString("PEPTIDE")
    pf = p.fromAASequence(aas)
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    assert result == "PEPTIDE"


def test_from_aasequence_with_modification():
    """Test creating peptidoform from AASequence with modification."""
    import pyopenms as p

    aas = p.AASequence.fromString("PEM(Oxidation)TIDE")
    pf = p.fromAASequence(aas)
    result = p.toString(pf, p.ProFormaWriteMode.LOSSLESS)
    # Should contain modification info
    assert "M" in result
    # Could be UNIMOD:35 or Oxidation depending on implementation
    assert "Oxidation" in result or "UNIMOD:35" in result


def test_aasequence_roundtrip():
    """Test AASequence -> Peptidoform -> AASequence roundtrip."""
    import pyopenms as p

    original = p.AASequence.fromString("PEPTIDE")
    pf = p.fromAASequence(original)
    result = p.toAASequence(pf, p.AASequenceConversionPolicy.STRICT)
    assert original.toString() == result.toString()


def test_get_conversion_issues():
    """Test getting list of conversion issues."""
    import pyopenms as p

    # Simple sequence should have no issues
    pf = p.parse("PEPTIDE")
    issues = p.getAASequenceConversionIssues(pf)
    assert len(issues) == 0

    # Unlocalised modification should have issues
    pf = p.parse("[Phospho]?PEPTIDE")
    issues = p.getAASequenceConversionIssues(pf)
    assert len(issues) > 0


def test_parse_error_handling():
    """Test that parse errors are properly raised."""
    import pyopenms as p

    # Unclosed bracket should fail
    with pytest.raises(Exception):
        p.parse("[Invalid")

    # Empty sequence should fail
    with pytest.raises(Exception):
        p.parse("")


def test_complex_proforma():
    """Test parsing complex ProForma strings."""
    import pyopenms as p

    # Global isotope label
    pf = p.parse("<13C>PEPTIDE")
    assert "13C" in p.toString(pf, p.ProFormaWriteMode.LOSSLESS)

    # Labile modification
    pf = p.parse("{Glycan:Hex}PEPTIDE")
    assert "Glycan" in p.toString(pf, p.ProFormaWriteMode.LOSSLESS)


def test_peptidoform_class():
    """Test Peptidoform class functionality."""
    import pyopenms as p

    pf = p.parse("PEPTIDE")

    # Should be a Peptidoform instance
    assert isinstance(pf, p.Peptidoform)

    # __repr__ should work
    repr_str = repr(pf)
    assert "Peptidoform" in repr_str


def test_peptidoform_ion_class():
    """Test PeptidoformIon class functionality."""
    import pyopenms as p

    pfi = p.parseIon("PEPTIDE/2")

    # Should be a PeptidoformIon instance
    assert isinstance(pfi, p.PeptidoformIon)

    # Should have chains
    assert len(pfi.chains) > 0


def test_write_mode_enum():
    """Test ProFormaWriteMode enum values."""
    import pyopenms as p

    assert hasattr(p.ProFormaWriteMode, 'LOSSLESS')
    assert hasattr(p.ProFormaWriteMode, 'CANONICAL')


def test_conversion_policy_enum():
    """Test AASequenceConversionPolicy enum values."""
    import pyopenms as p

    assert hasattr(p.AASequenceConversionPolicy, 'STRICT')
    assert hasattr(p.AASequenceConversionPolicy, 'DROP_UNLOCALISED')
    assert hasattr(p.AASequenceConversionPolicy, 'BEST_EFFORT')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
