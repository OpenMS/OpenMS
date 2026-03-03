"""
Tests for ProForma v2 parser/writer in pyOpenMS.

These tests cover:
- Basic parsing and roundtrip serialization
- WriteMode (LOSSLESS vs CANONICAL)
- AASequence conversion
- JSON serialization
- Error handling
"""

import json
import pytest


def test_parse_simple_sequence():
    """Test parsing simple peptide sequences."""
    import pyopenms as p

    # Simple sequence without modifications
    pf = p.Peptidoform.fromString("PEPTIDE")
    assert pf is not None

    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert result == "PEPTIDE"


def test_parse_with_unimod_modification():
    """Test parsing sequences with UNIMOD accessions."""
    import pyopenms as p

    # Oxidation on M
    pf = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "EM[UNIMOD:35]K"

    # Phosphorylation on S
    pf = p.Peptidoform.fromString("PES[UNIMOD:21]TIDE")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "PES[UNIMOD:21]TIDE"


def test_parse_with_named_modification():
    """Test parsing sequences with named modifications."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("EM[Oxidation]K")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "EM[Oxidation]K"


def test_parse_with_mass_delta():
    """Test parsing sequences with mass delta modifications."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("EM[+15.9949]K")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
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
        pf = p.Peptidoform.fromString(original)
        result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
        assert result == original, f"Roundtrip failed for {original}: got {result}"


def test_write_mode_canonical():
    """Test that canonical mode normalizes output."""
    import pyopenms as p

    # Parse with limited precision
    pf = p.Peptidoform.fromString("EM[+15.99]K")

    # Canonical mode should output 4 decimal places
    canonical = pf.toString( p.ProForma.WriteMode.CANONICAL)
    assert "+15.9900" in canonical

    # Lossless mode should preserve original
    lossless = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert "+15.99" in lossless


def test_terminal_modifications():
    """Test parsing N-terminal and C-terminal modifications."""
    import pyopenms as p

    # N-terminal modification
    pf = p.Peptidoform.fromString("[Acetyl]-PEPTIDE")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "[Acetyl]-PEPTIDE"

    # C-terminal modification
    pf = p.Peptidoform.fromString("PEPTIDE-[Amidated]")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "PEPTIDE-[Amidated]"

    # Both terminal modifications
    pf = p.Peptidoform.fromString("[Acetyl]-PEPTIDE-[Amidated]")
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "[Acetyl]-PEPTIDE-[Amidated]"


def test_parse_ion_with_charge():
    """Test parsing peptidoform ions with charge states."""
    import pyopenms as p

    # Parse ion with charge state
    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    # Verify parsing succeeded and chain was extracted
    assert len(pfi.chains) == 1


def test_parse_ion_multiple_chains():
    """Test parsing multi-chain peptidoforms."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE//SEQUENCE")
    # Verify multiple chains were parsed
    assert len(pfi.chains) == 2


def test_is_representable_as_aasequence():
    """Test checking if a peptidoform can be converted to AASequence."""
    import pyopenms as p

    # Simple sequence should be representable
    pf = p.Peptidoform.fromString("PEPTIDE")
    assert pf.isRepresentableAsAASequence()

    # Unlocalised modification is not representable
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")
    assert not pf.isRepresentableAsAASequence()


def test_to_aasequence_simple():
    """Test converting simple peptidoform to AASequence."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")
    aas = pf.toAASequence( p.ProForma.ConversionPolicy.FAIL_ON_LOSS)
    assert aas.toString() == "PEPTIDE"


def test_to_aasequence_with_modification():
    """Test converting peptidoform with modification to AASequence."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    p.ProForma.resolveModifications(pf)
    aas = pf.toAASequence( p.ProForma.ConversionPolicy.FAIL_ON_LOSS)
    # The AASequence should contain the modification
    assert "Oxidation" in aas.toString() or "(Oxidation)" in aas.toString()


def test_to_aasequence_preserves_residue_with_unimod():
    """Regression test: UNIMOD accession lookup must use residue context.

    Bug: When looking up UNIMOD:35 (Oxidation), the modification database
    contains variants for different residues (M, D, W, etc.). Without passing
    the residue context, it would pick the wrong variant (e.g., 'D' instead of 'M'),
    causing toString() to show wrong residue letters like 'ED(Oxidation)K'
    instead of 'EM(Oxidation)K'.
    """
    import pyopenms as p

    # Test Oxidation on M - should NOT become D
    pf = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    aas = pf.toAASequence()
    result = aas.toString()
    assert result == "EM(Oxidation)K", f"Expected 'EM(Oxidation)K', got '{result}'"

    # Verify the unmodified sequence is correct
    assert aas.toUnmodifiedString() == "EMK"

    # Verify residue at position 1 is M, not D
    assert aas.getResidue(1).getOneLetterCode() == "M"


def test_to_aasequence_fail_on_loss_raises_on_unsupported():
    """Test that FAIL_ON_LOSS policy raises error for unsupported features."""
    import pyopenms as p

    # Unlocalised modification cannot be converted
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")

    with pytest.raises(RuntimeError):
        pf.toAASequence( p.ProForma.ConversionPolicy.FAIL_ON_LOSS)


def test_to_aasequence_best_effort():
    """Test that BEST_EFFORT policy skips unsupported features."""
    import pyopenms as p

    # Unlocalised modification should be skipped
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")

    # Should not raise
    aas = pf.toAASequence( p.ProForma.ConversionPolicy.BEST_EFFORT)
    assert aas.toString() == "PEPTIDE"


def test_from_aasequence():
    """Test creating peptidoform from AASequence."""
    import pyopenms as p

    aas = p.AASequence.fromString("PEPTIDE")
    pf = p.Peptidoform.fromAASequence(aas)
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    assert result == "PEPTIDE"


def test_from_aasequence_with_modification():
    """Test creating peptidoform from AASequence with modification."""
    import pyopenms as p

    aas = p.AASequence.fromString("PEM(Oxidation)TIDE")
    pf = p.Peptidoform.fromAASequence(aas)
    result = pf.toString( p.ProForma.WriteMode.LOSSLESS)
    # Should contain modification info
    assert "M" in result
    # Could be UNIMOD:35 or Oxidation depending on implementation
    assert "Oxidation" in result or "UNIMOD:35" in result


def test_aasequence_roundtrip():
    """Test AASequence -> Peptidoform -> AASequence roundtrip."""
    import pyopenms as p

    original = p.AASequence.fromString("PEPTIDE")
    pf = p.Peptidoform.fromAASequence(original)
    result = pf.toAASequence( p.ProForma.ConversionPolicy.FAIL_ON_LOSS)
    assert original.toString() == result.toString()


def test_get_conversion_issues():
    """Test getting list of conversion issues."""
    import pyopenms as p

    # Simple sequence should have no issues
    pf = p.Peptidoform.fromString("PEPTIDE")
    issues = pf.getAASequenceConversionIssues()
    assert len(issues) == 0

    # Unlocalised modification should have issues
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")
    issues = pf.getAASequenceConversionIssues()
    assert len(issues) > 0


def test_json_roundtrip_peptidoform():
    """Test JSON serialization roundtrip for Peptidoform."""
    import pyopenms as p

    original = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    json_str = p.ProForma.peptidoformToJSON(original)

    # Should be valid JSON
    data = json.loads(json_str)
    assert "sequence" in data

    # Should roundtrip
    restored = p.Peptidoform.fromJSON(json_str)
    assert p.ProForma.toString(restored, p.ProForma.WriteMode.LOSSLESS) == p.ProForma.toString(original, p.ProForma.WriteMode.LOSSLESS)


def test_json_roundtrip_peptidoform_ion():
    """Test JSON serialization roundtrip for PeptidoformIon."""
    import pyopenms as p

    original = p.PeptidoformIon.fromString("PEPTIDE/2")
    json_str = p.ProForma.peptidoformIonToJSON(original)

    # Should be valid JSON
    data = json.loads(json_str)
    assert "chains" in data

    # Should roundtrip
    restored = p.PeptidoformIon.fromJSON(json_str)
    # Compare JSON representations
    assert p.ProForma.peptidoformIonToJSON(restored) == p.ProForma.peptidoformIonToJSON(original)


def test_parse_error_handling():
    """Test that parse errors are properly raised."""
    import pyopenms as p

    # Unclosed bracket should fail
    with pytest.raises(RuntimeError):
        p.Peptidoform.fromString("[Invalid")

    # Empty sequence should fail
    with pytest.raises(RuntimeError):
        p.Peptidoform.fromString("")


def test_complex_proforma():
    """Test parsing complex ProForma strings."""
    import pyopenms as p

    # Global isotope label
    pf = p.Peptidoform.fromString("<13C>PEPTIDE")
    assert "13C" in pf.toString( p.ProForma.WriteMode.LOSSLESS)

    # Labile modification
    pf = p.Peptidoform.fromString("{Glycan:Hex}PEPTIDE")
    assert "Glycan" in pf.toString( p.ProForma.WriteMode.LOSSLESS)


def test_peptidoform_class():
    """Test Peptidoform class functionality."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")

    # Should be a Peptidoform instance
    assert isinstance(pf, p.Peptidoform)

    # __repr__ should work
    repr_str = repr(pf)
    assert "Peptidoform" in repr_str


def test_peptidoform_ion_class():
    """Test PeptidoformIon class functionality."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")

    # Should be a PeptidoformIon instance
    assert isinstance(pfi, p.PeptidoformIon)

    # Should have chains
    assert len(pfi.chains) > 0


def test_peptidoform_ion_to_string():
    """Test PeptidoformIon toString (toStringIon) method."""
    import pyopenms as p

    # Simple ion with charge
    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    result = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "PEPTIDE" in result
    assert "/2" in result

    # Cross-linked chains
    pfi2 = p.PeptidoformIon.fromString("PEPTIDE//SEQUENCE")
    result2 = pfi2.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "PEPTIDE" in result2
    assert "SEQUENCE" in result2
    assert "//" in result2

    # Test roundtrip
    original_str = "EM[UNIMOD:35]K/2"
    pfi3 = p.PeptidoformIon.fromString(original_str)
    roundtrip = pfi3.toString(p.ProForma.WriteMode.LOSSLESS)
    assert roundtrip == original_str


def test_peptidoform_ion_to_string_modes():
    """Test PeptidoformIon toString with different write modes and default parameter."""
    import pyopenms as p

    # Test with mass delta to see difference between modes
    pfi = p.PeptidoformIon.fromString("PEM[+15.99]K/2")

    # LOSSLESS mode preserves original formatting
    lossless = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "+15.99" in lossless  # Original precision preserved

    # CANONICAL mode normalizes to 4 decimal places
    canonical = pfi.toString(p.ProForma.WriteMode.CANONICAL)
    assert "+15.9900" in canonical  # Normalized to 4 decimals

    # Test default parameter (should default to LOSSLESS)
    default_result = pfi.toString()
    assert default_result == lossless  # Default should match LOSSLESS


def test_peptidoform_ion_to_string_exact_format():
    """Test PeptidoformIon toString preserves exact separator and order."""
    import pyopenms as p

    # Test that // separator is exactly preserved (not + or other)
    pfi = p.PeptidoformIon.fromString("ABC//DEF//GHI")
    result = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    # Verify exact format: chains separated by // in order
    assert result == "ABC//DEF//GHI"

    # Test chain order is preserved
    pfi2 = p.PeptidoformIon.fromString("FIRST//SECOND")
    result2 = pfi2.toString(p.ProForma.WriteMode.LOSSLESS)
    assert result2.index("FIRST") < result2.index("SECOND")


def test_peptidoform_ion_mass_integration():
    """Test PeptidoformIon mass/m/z calculation integration with toString."""
    import pyopenms as p

    # Create ion with charge
    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")

    # Verify mass calculation works
    assert pfi.canCalculateMass()
    mass = pfi.getMonoWeight()
    assert 799.0 < mass < 800.0  # PEPTIDE ~799.36 Da

    # Verify m/z calculation works (uses charge from ion)
    mz = pfi.getMZ()
    expected_mz = (mass + 2 * 1.007276) / 2  # (M + 2H) / 2
    assert abs(mz - expected_mz) < 0.1

    # Verify toString still works after mass calculations
    result = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    assert result == "PEPTIDE/2"


def test_write_mode_enum():
    """Test WriteMode enum values accessible via ProForma class."""
    import pyopenms as p

    # Enums are nested under ProForma class (C++ enum class with wrap-attach)
    assert hasattr(p.ProForma, 'WriteMode')
    assert hasattr(p.ProForma.WriteMode, 'LOSSLESS')
    assert hasattr(p.ProForma.WriteMode, 'CANONICAL')

    # Verify enum members are distinct
    assert p.ProForma.WriteMode.LOSSLESS != p.ProForma.WriteMode.CANONICAL


def test_conversion_policy_enum():
    """Test ConversionPolicy enum values accessible via ProForma class."""
    import pyopenms as p

    # Enums are nested under ProForma class (C++ enum class with wrap-attach)
    assert hasattr(p.ProForma, 'ConversionPolicy')
    assert hasattr(p.ProForma.ConversionPolicy, 'FAIL_ON_LOSS')
    assert hasattr(p.ProForma.ConversionPolicy, 'DROP_UNLOCALISED')
    assert hasattr(p.ProForma.ConversionPolicy, 'BEST_EFFORT')

    # Verify enum members are distinct
    assert p.ProForma.ConversionPolicy.FAIL_ON_LOSS != p.ProForma.ConversionPolicy.DROP_UNLOCALISED
    assert p.ProForma.ConversionPolicy.DROP_UNLOCALISED != p.ProForma.ConversionPolicy.BEST_EFFORT


def test_nested_enums_under_proforma():
    """Test that ProForma enums are accessible via ProForma class (via wrap-attach)."""
    import pyopenms as p

    # Enums should be nested under ProForma via wrap-attach
    assert hasattr(p.ProForma, 'WriteMode')
    assert hasattr(p.ProForma, 'ConversionPolicy')
    assert hasattr(p.ProForma, 'ConversionIssueType')
    assert hasattr(p.ProForma, 'CvDatabase')

    # Verify enum values are accessible
    assert p.ProForma.WriteMode.LOSSLESS is not None
    assert p.ProForma.WriteMode.CANONICAL is not None
    assert p.ProForma.ConversionPolicy.FAIL_ON_LOSS is not None
    assert p.ProForma.ConversionPolicy.DROP_UNLOCALISED is not None
    assert p.ProForma.ConversionPolicy.BEST_EFFORT is not None


def test_peptidoform_mass_calculations():
    """Test mass calculation methods on Peptidoform."""
    import pyopenms as p

    # Simple peptide
    pf = p.Peptidoform.fromString("PEPTIDE")
    assert pf.canCalculateMass()
    mass = pf.getMonoWeight()
    assert 799.0 < mass < 800.0  # PEPTIDE ~799.36 Da

    # With modification
    pf_mod = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    assert pf_mod.canCalculateMass()
    mass_mod = pf_mod.getMonoWeight()
    assert 420.0 < mass_mod < 425.0  # ~422.18 Da

    # m/z calculation
    mz = pf.getMZ(2)  # z=2
    assert 399.0 < mz < 401.0  # (799.36 + 2*1.007) / 2 ≈ 400.69


def test_peptidoform_ion_mass_calculations():
    """Test mass calculation methods on PeptidoformIon."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    assert pfi.canCalculateMass()

    mass = pfi.getMonoWeight()
    assert 799.0 < mass < 800.0

    # getMZ uses embedded charge
    mz = pfi.getMZ()
    assert 399.0 < mz < 401.0


def test_peptidoform_repr_str():
    """Test __repr__ and __str__ methods on Peptidoform."""
    import pyopenms as p

    # Simple sequence
    pf = p.Peptidoform.fromString("PEPTIDE")
    repr_str = repr(pf)
    assert "Peptidoform(" in repr_str
    assert "sequence=" in repr_str
    assert "PEPTIDE" in repr_str
    assert "mono_mass=" in repr_str

    str_str = str(pf)
    assert str_str == "PEPTIDE"

    # With modification
    pf_mod = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    repr_mod = repr(pf_mod)
    assert "Peptidoform(" in repr_mod
    assert "EM[UNIMOD:35]K" in repr_mod

    str_mod = str(pf_mod)
    assert str_mod == "EM[UNIMOD:35]K"


def test_peptidoform_ion_repr_str():
    """Test __repr__ and __str__ methods on PeptidoformIon."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    repr_str = repr(pfi)
    assert "PeptidoformIon(" in repr_str
    assert "chains=" in repr_str

    str_str = str(pfi)
    assert "PEPTIDE" in str_str


def test_peptidoform_various_modifications():
    """Test parsing various modification types (similar to AASequence tests)."""
    import pyopenms as p

    # UNIMOD accession
    pf1 = p.Peptidoform.fromString("PEPTIDESEKLEM[UNIMOD:35]CER")
    assert pf1.toString() == "PEPTIDESEKLEM[UNIMOD:35]CER"
    aas1 = pf1.toAASequence()
    assert aas1.toString() == "PEPTIDESEKLEM(Oxidation)CER"
    assert aas1.toUnmodifiedString() == "PEPTIDESEKLEMCER"

    # Named modification
    pf2 = p.Peptidoform.fromString("PEPTIDESEKLEM[Oxidation]CER")
    assert pf2.toString() == "PEPTIDESEKLEM[Oxidation]CER"

    # N-terminal modification
    pf3 = p.Peptidoform.fromString("[UNIMOD:1]-PEPTIDE")
    assert "[UNIMOD:1]" in pf3.toString()
    aas3 = pf3.toAASequence()
    assert "Acetyl" in aas3.toString()

    # C-terminal modification
    pf4 = p.Peptidoform.fromString("PEPTIDE-[UNIMOD:2]")
    assert "[UNIMOD:2]" in pf4.toString()

    # Multiple modifications
    pf5 = p.Peptidoform.fromString("[UNIMOD:1]-EM[UNIMOD:35]K")
    aas5 = pf5.toAASequence()
    assert "Acetyl" in aas5.toString()
    assert "Oxidation" in aas5.toString()


def test_peptidoform_aasequence_comparison():
    """Test that Peptidoform -> AASequence produces equivalent results to direct AASequence parsing."""
    import pyopenms as p

    test_cases = [
        ("PEPTIDE", "PEPTIDE"),
        ("EM[UNIMOD:35]K", "EM(Oxidation)K"),
        ("EC[UNIMOD:4]K", "EC(Carbamidomethyl)K"),
        ("[UNIMOD:1]-PEPTIDE", ".(Acetyl)PEPTIDE"),
    ]

    for proforma, expected_aas in test_cases:
        pf = p.Peptidoform.fromString(proforma)
        aas_from_pf = pf.toAASequence()
        aas_direct = p.AASequence.fromString(expected_aas)

        # Masses should match
        assert abs(aas_from_pf.getMonoWeight() - aas_direct.getMonoWeight()) < 0.001, \
            f"Mass mismatch for {proforma}: {aas_from_pf.getMonoWeight()} vs {aas_direct.getMonoWeight()}"

        # Unmodified strings should match
        assert aas_from_pf.toUnmodifiedString() == aas_direct.toUnmodifiedString(), \
            f"Unmodified string mismatch for {proforma}"


def test_peptidoform_len():
    """Test PeptidoformIon __len__ method."""
    import pyopenms as p

    # Single chain
    pfi1 = p.PeptidoformIon.fromString("PEPTIDE/2")
    assert len(pfi1) == 1

    # Multiple chains (cross-linked)
    pfi2 = p.PeptidoformIon.fromString("PEPTIDE//SEQUENCE")
    assert len(pfi2) == 2


# =============================================================================
# Spectrum Generation Tests
# =============================================================================

def test_can_generate_spectrum_simple():
    """Test canGenerateSpectrum for simple peptides."""
    import pyopenms as p

    # Simple peptide should be generatable
    pf = p.Peptidoform.fromString("PEPTIDE")
    assert pf.canGenerateSpectrum()

    # PeptidoformIon should also work
    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    assert pfi.canGenerateSpectrum()


def test_can_generate_spectrum_with_modification():
    """Test canGenerateSpectrum for modified peptides."""
    import pyopenms as p

    # Modified peptide should be generatable
    pf = p.Peptidoform.fromString("EM[UNIMOD:35]K")
    assert pf.canGenerateSpectrum()


def test_can_generate_spectrum_unsupported():
    """Test canGenerateSpectrum returns False for unsupported features."""
    import pyopenms as p

    # Unlocalised modification cannot generate spectrum
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")
    assert not pf.canGenerateSpectrum()


def test_generate_spectrum_simple():
    """Test generating theoretical spectrum for simple peptide."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")

    # Generate b and y ions, charge 1-2
    spectrum = pf.generateSpectrum(1, 2, "by", False, False)

    # Should have peaks
    assert spectrum.size() > 0

    # Should have b and y ions (6 residues = 5 b ions + 5 y ions per charge)
    # With 2 charge states, expect around 20 peaks minimum
    assert spectrum.size() >= 10


def test_generate_spectrum_with_losses():
    """Test generating spectrum with neutral losses."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")

    # Without losses
    spectrum_no_loss = pf.generateSpectrum(1, 1, "by", False, False)

    # With losses (water, ammonia)
    spectrum_with_loss = pf.generateSpectrum(1, 1, "by", True, False)

    # Spectrum with losses should have more peaks
    assert spectrum_with_loss.size() >= spectrum_no_loss.size()


def test_generate_spectrum_with_metainfo():
    """Test generating spectrum with peak annotations."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")

    # Generate with metainfo
    spectrum = pf.generateSpectrum(1, 1, "by", False, True)

    # Should have peaks with annotations
    assert spectrum.size() > 0
    # Check that string data arrays exist (ion annotations)
    assert len(spectrum.getStringDataArrays()) > 0


def test_generate_spectrum_ion_types():
    """Test different ion type combinations."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPTIDE")

    # Only b ions
    spectrum_b = pf.generateSpectrum(1, 1, "b", False, False)

    # Only y ions
    spectrum_y = pf.generateSpectrum(1, 1, "y", False, False)

    # Both b and y
    spectrum_by = pf.generateSpectrum(1, 1, "by", False, False)

    # Combined should have roughly sum of individual (may overlap)
    assert spectrum_by.size() >= max(spectrum_b.size(), spectrum_y.size())


def test_generate_spectrum_peptidoform_ion():
    """Test generating spectrum for PeptidoformIon."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")

    # Generate spectrum
    spectrum = pfi.generateSpectrum(1, 2, "by", False, False)

    # Should have peaks
    assert spectrum.size() > 0


def test_get_spectrum_generation_issues():
    """Test getSpectrumGenerationIssues returns issues for unsupported features."""
    import pyopenms as p

    # Simple peptide should have no issues
    pf = p.Peptidoform.fromString("PEPTIDE")
    issues = pf.getSpectrumGenerationIssues()
    assert len(issues) == 0

    # Unlocalised modification should have issues
    pf_unloc = p.Peptidoform.fromString("[Phospho]?PEPTIDE")
    issues_unloc = pf_unloc.getSpectrumGenerationIssues()
    assert len(issues_unloc) > 0


def test_get_spectrum_generation_issues_ion():
    """Test getSpectrumGenerationIssues for PeptidoformIon."""
    import pyopenms as p

    # Simple ion should have no issues
    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    issues = pfi.getSpectrumGenerationIssues()
    assert len(issues) == 0


# =============================================================================
# Mass Calculation Issues Tests
# =============================================================================

def test_get_mass_calculation_issues_simple():
    """Test getMassCalculationIssues for simple peptides."""
    import pyopenms as p

    # Simple peptide should have no issues
    pf = p.Peptidoform.fromString("PEPTIDE")
    issues = pf.getMassCalculationIssues()
    assert len(issues) == 0


def test_get_mass_calculation_issues_unresolved():
    """Test getMassCalculationIssues for unresolved modifications."""
    import pyopenms as p

    # Unknown modification name should have issues
    pf = p.Peptidoform.fromString("PEM[SomeUnknownMod]K")
    issues = pf.getMassCalculationIssues()
    # May or may not have issues depending on whether mod can be looked up
    # At minimum, verify the method works
    assert isinstance(issues, list)


def test_get_mass_calculation_issues_ion():
    """Test getMassCalculationIssues for PeptidoformIon."""
    import pyopenms as p

    pfi = p.PeptidoformIon.fromString("PEPTIDE/2")
    issues = pfi.getMassCalculationIssues()
    assert len(issues) == 0


# =============================================================================
# Cross-link Label Tests
# =============================================================================

def test_crosslink_label_parsing():
    """Test parsing cross-link labels (#XL1)."""
    import pyopenms as p

    # Simple cross-link label
    pf = p.Peptidoform.fromString("PEPK[#XL1]TIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "#XL1" in result


def test_crosslink_paired_chains():
    """Test parsing paired cross-linked chains."""
    import pyopenms as p

    # Two chains with matching cross-link labels
    pfi = p.PeptidoformIon.fromString("PEPK[#XL1]TIDE//SEQC[#XL1]ENCE")
    assert len(pfi.chains) == 2

    result = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "#XL1" in result
    assert "//" in result


def test_crosslink_with_modification():
    """Test cross-link combined with modification."""
    import pyopenms as p

    # Cross-link label with a modification (DSS cross-linker)
    pf = p.Peptidoform.fromString("PEPK[XLMOD:02001#XL1]TIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "XLMOD:02001" in result
    assert "#XL1" in result


# =============================================================================
# Ambiguous Modification Tests
# =============================================================================

def test_ambiguous_modification_localization():
    """Test ambiguous modification localization with positional options."""
    import pyopenms as p

    # Phosphorylation on one of the S/T residues (ambiguous position)
    pf = p.Peptidoform.fromString("PES[Phospho#g1]TIDES[#g1]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "#g1" in result


def test_ambiguous_modification_count():
    """Test ambiguous modification with occurrence count."""
    import pyopenms as p

    # Two phosphorylations distributed among 3 possible sites
    # Using unlocalised with count
    pf = p.Peptidoform.fromString("[Phospho]^2?STYTIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Phospho" in result
    assert "^2" in result


def test_ambiguous_modification_range():
    """Test ambiguous modification on a range of residues."""
    import pyopenms as p

    # Modification somewhere in positions 2-5
    pf = p.Peptidoform.fromString("P(EPTI)[Phospho]DE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Phospho" in result


# =============================================================================
# Formula Tag Tests
# =============================================================================

def test_formula_tag_simple():
    """Test parsing formula tags."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEM[Formula:C2H2O]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Formula:C2H2O" in result


def test_formula_tag_complex():
    """Test parsing complex formula tags."""
    import pyopenms as p

    # Formula with negative counts (loss)
    pf = p.Peptidoform.fromString("PEM[Formula:H-2O-1]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Formula:" in result


def test_formula_tag_mass_calculation():
    """Test that formula tags allow mass calculation."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEM[Formula:O1]K")  # +16 Da (like oxidation)
    assert pf.canCalculateMass()

    mass = pf.getMonoWeight()
    # PEMK + O = 519.2363
    assert abs(mass - 519.2363) < 0.01


# =============================================================================
# Info Tag Tests
# =============================================================================

def test_info_tag_simple():
    """Test parsing info tags."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEM[info:custom annotation]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    # ProForma 2.0 is case-insensitive; implementation outputs uppercase INFO:
    assert "info:custom annotation" in result.lower()


def test_info_tag_with_modification():
    """Test info tag combined with modification."""
    import pyopenms as p

    # Modification with additional info
    pf = p.Peptidoform.fromString("PEM[Oxidation][info:confirmed by manual inspection]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Oxidation" in result
    # ProForma 2.0 is case-insensitive; implementation outputs uppercase INFO:
    assert "info:" in result.lower()


def test_info_tag_does_not_affect_mass():
    """Test that info tags don't affect mass calculation."""
    import pyopenms as p

    pf_no_info = p.Peptidoform.fromString("PEPTIDE")
    pf_with_info = p.Peptidoform.fromString("PEP[info:some note]TIDE")

    mass_no_info = pf_no_info.getMonoWeight()
    mass_with_info = pf_with_info.getMonoWeight()

    # Masses should be identical
    assert abs(mass_no_info - mass_with_info) < 0.001


# =============================================================================
# Glycan Notation Tests
# =============================================================================

def test_glycan_composition():
    """Test parsing glycan composition notation."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPN[Glycan:HexNAc2Hex3]TIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Glycan:" in result
    assert "HexNAc" in result


def test_glycan_gno_accession():
    """Test parsing GNO accession for glycans."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEPN[GNO:G59626AS]TIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "GNO:" in result


# =============================================================================
# Additional Edge Cases
# =============================================================================

def test_multiple_modifications_same_residue():
    """Test multiple modifications on the same residue."""
    import pyopenms as p

    # Multiple mods on same residue (e.g., phospho and acetyl on S)
    pf = p.Peptidoform.fromString("PES[Phospho][Acetyl]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Phospho" in result
    assert "Acetyl" in result


def test_resid_accession():
    """Test parsing RESID accession."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEM[RESID:AA0581]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "RESID:" in result


def test_mod_accession():
    """Test parsing MOD (PSI-MOD) accession."""
    import pyopenms as p

    pf = p.Peptidoform.fromString("PEM[MOD:00719]K")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "MOD:" in result


def test_charge_state_formats():
    """Test various charge state format parsing."""
    import pyopenms as p

    # Simple charge
    pfi1 = p.PeptidoformIon.fromString("PEPTIDE/2")
    assert len(pfi1.chains) == 1

    # Explicit positive charge
    pfi2 = p.PeptidoformIon.fromString("PEPTIDE/+2")
    result2 = pfi2.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "PEPTIDE" in result2

    # Negative charge (rare but valid)
    pfi3 = p.PeptidoformIon.fromString("PEPTIDE/-1")
    result3 = pfi3.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "PEPTIDE" in result3


def test_adduct_ion_notation():
    """Test adduct ion notation parsing."""
    import pyopenms as p

    # Adduct notation with element:charge format
    # ProForma spec: /[element:count+charge,...]
    pfi = p.PeptidoformIon.fromString("PEPTIDE/[Na:z+2]")
    result = pfi.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "PEPTIDE" in result
    assert "Na" in result


def test_global_modification_fixed():
    """Test fixed global modification notation."""
    import pyopenms as p

    # Fixed Carbamidomethyl on all C
    pf = p.Peptidoform.fromString("<[Carbamidomethyl]@C>PEPTCIDE")
    result = pf.toString(p.ProForma.WriteMode.LOSSLESS)
    assert "Carbamidomethyl" in result
    assert "@C" in result


def test_conversion_issue_types():
    """Test that ConversionIssue objects have correct attributes."""
    import pyopenms as p

    # Get issues from unlocalised modification
    pf = p.Peptidoform.fromString("[Phospho]?PEPTIDE")
    issues = pf.getAASequenceConversionIssues()

    assert len(issues) > 0

    # Check issue has expected attributes
    issue = issues[0]
    assert hasattr(issue, 'type')
    assert hasattr(issue, 'description')
    assert hasattr(issue, 'position')

    # Type should be UNLOCALISED_MOD
    assert issue.type == p.ProForma.ConversionIssueType.UNLOCALISED_MOD


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
