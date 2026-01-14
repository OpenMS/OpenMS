#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tutorial: ProForma v2 Parser and Writer

This tutorial demonstrates how to use the ProForma v2 parser and writer
in pyOpenMS for representing peptides with modifications.

ProForma is a standardized notation for proteoforms developed by the
HUPO Proteomics Standards Initiative (PSI).
"""

import pyopenms as oms


def basic_parsing():
    """Basic parsing and serialization examples."""
    print("=== Basic Parsing ===")

    # Parse a simple peptide
    pf = oms.parse("PEPTIDE")
    print(f"Simple: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Parse with UNIMOD modification (Oxidation on M)
    pf = oms.parse("PEM[UNIMOD:35]TIDE")
    print(f"UNIMOD: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Parse with named modification
    pf = oms.parse("PES[Phospho]TIDE")
    print(f"Named:  {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Parse with mass delta
    pf = oms.parse("PEM[+15.9949]TIDE")
    print(f"Delta:  {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")
    print()


def terminal_modifications():
    """Terminal modification examples."""
    print("=== Terminal Modifications ===")

    # N-terminal acetylation
    pf = oms.parse("[Acetyl]-PEPTIDE")
    print(f"N-term: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # C-terminal amidation
    pf = oms.parse("PEPTIDE-[Amidated]")
    print(f"C-term: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Both terminals modified
    pf = oms.parse("[Acetyl]-PEPTIDE-[Amidated]")
    print(f"Both:   {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")
    print()


def write_modes():
    """Demonstrate LOSSLESS vs CANONICAL write modes."""
    print("=== Write Modes ===")

    # Parse with limited precision mass delta
    pf = oms.parse("EM[+15.99]K")

    # Lossless preserves original formatting
    lossless = oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)
    print(f"Lossless:  {lossless}")

    # Canonical normalizes to 4 decimal places
    canonical = oms.toString(pf, oms.ProFormaWriteMode.CANONICAL)
    print(f"Canonical: {canonical}")
    print()


def advanced_features():
    """Advanced ProForma features."""
    print("=== Advanced Features ===")

    # Global isotope labeling (SILAC)
    pf = oms.parse("<13C>PEPTIDE")
    print(f"Isotope label: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Global modification on specific residues (TMT)
    pf = oms.parse("<[TMT6plex]@K,N-term>PEPTIDEK")
    print(f"TMT labeling: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Unlocalised modification (position unknown)
    pf = oms.parse("[Phospho]?STYPEPTIDE")
    print(f"Unlocalised: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")

    # Labile modification (lost during fragmentation)
    pf = oms.parse("{Glycan:Hex}PEPTIDE")
    print(f"Labile: {oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)}")
    print()


def multi_chain_and_charge():
    """Multi-chain peptides and charge states."""
    print("=== Multi-chain and Charge States ===")

    # Parse with charge state
    pfi = oms.parseIon("PEPTIDE/2")
    print(f"Charged ion - chains: {len(pfi.chains)}")

    # Multi-chain (cross-linked peptides)
    pfi = oms.parseIon("PEPTIDECK[XLMOD:02001#XL1]//AC[#XL1]DEFGH")
    print(f"Cross-linked - chains: {len(pfi.chains)}")

    # Adduct charge notation
    pfi = oms.parseIon("PEPTIDE/[Na:z+1,H:z+1]")
    print(f"Adduct notation - chains: {len(pfi.chains)}")
    print()


def aasequence_conversion():
    """Converting between ProForma and AASequence."""
    print("=== AASequence Conversion ===")

    # Parse and check if convertible
    pf = oms.parse("PEM[UNIMOD:35]TIDE")

    if oms.isRepresentableAsAASequence(pf):
        # Resolve modifications using ModificationsDB
        oms.resolveModifications(pf)

        # Convert to AASequence
        aas = oms.toAASequence(pf, oms.AASequenceConversionPolicy.STRICT)
        print(f"AASequence: {aas.toString()}")
        print(f"Mono mass:  {aas.getMonoWeight():.4f} Da")

    # Check conversion issues for complex peptidoforms
    pf = oms.parse("[Phospho]?STYPEPTIDE")
    issues = oms.getAASequenceConversionIssues(pf)
    print(f"\nConversion issues for '[Phospho]?STYPEPTIDE':")
    for issue in issues:
        print(f"  - {issue.description}")

    # Use BEST_EFFORT to convert what's possible
    aas = oms.toAASequence(pf, oms.AASequenceConversionPolicy.BEST_EFFORT)
    print(f"Best effort: {aas.toString()}")
    print()


def aasequence_to_proforma():
    """Converting AASequence to ProForma."""
    print("=== AASequence to ProForma ===")

    # Create AASequence the traditional way
    aas = oms.AASequence.fromString("PEM(Oxidation)TIDE")
    print(f"AASequence: {aas.toString()}")

    # Convert to ProForma
    pf = oms.fromAASequence(aas)
    proforma_str = oms.toString(pf, oms.ProFormaWriteMode.CANONICAL)
    print(f"ProForma:   {proforma_str}")
    print()


def json_serialization():
    """JSON serialization for storage and interchange."""
    print("=== JSON Serialization ===")
    import json

    # Serialize Peptidoform to JSON
    pf = oms.parse("[Acetyl]-PEM[UNIMOD:35]TIDE")
    json_str = oms.peptidoformToJSON(pf)
    print("Peptidoform JSON (truncated):")
    print(f"  {json_str[:80]}...")

    # Deserialize from JSON
    restored = oms.peptidoformFromJSON(json_str)
    original = oms.toString(pf, oms.ProFormaWriteMode.LOSSLESS)
    roundtrip = oms.toString(restored, oms.ProFormaWriteMode.LOSSLESS)
    print(f"Roundtrip OK: {original == roundtrip}")

    # PeptidoformIon JSON
    pfi = oms.parseIon("PEPTIDE/2")
    json_str = oms.peptidoformIonToJSON(pfi)
    restored = oms.peptidoformIonFromJSON(json_str)
    print(f"Ion roundtrip OK: {len(pfi.chains) == len(restored.chains)}")
    print()


def batch_processing():
    """Process multiple peptides and calculate masses."""
    print("=== Batch Processing ===")

    peptides = [
        "PEPTIDE",
        "PEM[UNIMOD:35]TIDE",
        "[Acetyl]-PEPTIDEK",
        "PEPS[UNIMOD:21]TIDE",
    ]

    print(f"{'ProForma':<30} {'Mass (Da)':>12}")
    print("-" * 44)

    for proforma_str in peptides:
        pf = oms.parse(proforma_str)

        if oms.isRepresentableAsAASequence(pf):
            oms.resolveModifications(pf)
            aas = oms.toAASequence(pf, oms.AASequenceConversionPolicy.STRICT)
            mass = aas.getMonoWeight()
            print(f"{proforma_str:<30} {mass:>12.4f}")
        else:
            print(f"{proforma_str:<30} {'N/A':>12}")
    print()


def error_handling():
    """Demonstrate error handling."""
    print("=== Error Handling ===")

    invalid_inputs = [
        ("[Invalid", "Unclosed bracket"),
        ("", "Empty input"),
        ("PEP[", "Incomplete modification"),
    ]

    for proforma_str, description in invalid_inputs:
        try:
            pf = oms.parse(proforma_str)
            print(f"  '{proforma_str}': Unexpectedly succeeded")
        except Exception as e:
            print(f"  {description}: Caught exception (expected)")
    print()


if __name__ == "__main__":
    print("=" * 60)
    print("ProForma v2 Tutorial for pyOpenMS")
    print("=" * 60)
    print()

    basic_parsing()
    terminal_modifications()
    write_modes()
    advanced_features()
    multi_chain_and_charge()
    aasequence_conversion()
    aasequence_to_proforma()
    json_serialization()
    batch_processing()
    error_handling()

    print("Tutorial completed successfully!")
