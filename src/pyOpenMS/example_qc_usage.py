#!/usr/bin/env python
"""
Example usage of QC classes in pyOpenMS.

This script demonstrates how to use the QC metrics classes that are now
available in pyOpenMS for quality control of mass spectrometry data.
"""

import pyopenms as oms


def example_tic_computation():
    """
    Example: Compute Total Ion Count (TIC) from an MSExperiment
    """
    print("="*60)
    print("Example 1: Computing TIC from MSExperiment")
    print("="*60)
    
    # Create a TIC object
    tic = oms.TIC()
    
    # Load an mzML file (or create a dummy experiment)
    exp = oms.MSExperiment()
    
    # For demonstration, let's create a simple spectrum
    spec = oms.MSSpectrum()
    spec.setRT(100.0)
    spec.setMSLevel(1)
    
    # Add some peaks
    peak1 = oms.Peak1D()
    peak1.setMZ(500.0)
    peak1.setIntensity(1000.0)
    spec.push_back(peak1)
    
    peak2 = oms.Peak1D()
    peak2.setMZ(600.0)
    peak2.setIntensity(2000.0)
    spec.push_back(peak2)
    
    exp.addSpectrum(spec)
    
    # Compute TIC
    result = tic.compute(exp, 0.0, 1)  # bin_size=0, ms_level=1
    
    print(f"✓ TIC computed successfully")
    print(f"  Area under TIC: {result.area}")
    print(f"  Number of retention time points: {len(result.retention_times)}")
    print(f"  Number of intensity points: {len(result.intensities)}")


def example_peptide_mass():
    """
    Example: Annotate peptide masses in a FeatureMap
    """
    print("\n" + "="*60)
    print("Example 2: Computing Peptide Masses")
    print("="*60)
    
    # Create a PeptideMass object
    peptide_mass = oms.PeptideMass()
    
    # Create a simple FeatureMap
    feature_map = oms.FeatureMap()
    
    # Note: In real usage, you would load a featureXML file with:
    # oms.FeatureXMLFile().load("input.featureXML", feature_map)
    
    # The compute method will add 'mass' metavalues to PeptideHits
    # peptide_mass.compute(feature_map)
    
    print(f"✓ PeptideMass object created: {peptide_mass.getName()}")
    print("  Use compute(feature_map) to annotate theoretical masses")


def example_fwhm():
    """
    Example: Transfer FWHM values from features to PeptideIdentifications
    """
    print("\n" + "="*60)
    print("Example 3: Transferring FWHM values")
    print("="*60)
    
    # Create an FWHM object
    fwhm = oms.FWHM()
    
    # Create a FeatureMap
    feature_map = oms.FeatureMap()
    
    # The compute method transfers FWHM metavalues
    # fwhm.compute(feature_map)
    
    print(f"✓ FWHM object created: {fwhm.getName()}")
    print("  Use compute(feature_map) to transfer FWHM values")


def example_ms2_identification_rate():
    """
    Example: Calculate MS2 identification rate
    """
    print("\n" + "="*60)
    print("Example 4: Computing MS2 Identification Rate")
    print("="*60)
    
    # Create an Ms2IdentificationRate object
    ms2_rate = oms.Ms2IdentificationRate()
    
    print(f"✓ Ms2IdentificationRate object created: {ms2_rate.getName()}")
    print("  This metric computes the ratio of identified MS2 spectra")
    print("  Usage: ms2_rate.compute(feature_map, experiment, assume_all_target=False)")


def example_rt_alignment():
    """
    Example: Annotate retention times after alignment
    """
    print("\n" + "="*60)
    print("Example 5: RT Alignment Annotation")
    print("="*60)
    
    # Create an RTAlignment object
    rt_alignment = oms.RTAlignment()
    
    print(f"✓ RTAlignment object created: {rt_alignment.getName()}")
    print("  This metric annotates raw and aligned retention times")
    print("  Usage: rt_alignment.compute(feature_map, transformation_description)")


def example_contaminants():
    """
    Example: Check for contaminant peptides
    """
    print("\n" + "="*60)
    print("Example 6: Contaminant Detection")
    print("="*60)
    
    # Create a Contaminants object
    contaminants = oms.Contaminants()
    
    print(f"✓ Contaminants object created: {contaminants.getName()}")
    print("  This metric identifies contaminant peptides")
    print("  Usage: contaminants.compute(feature_map, fasta_entries)")


def main():
    """
    Main function demonstrating all QC class examples
    """
    print("\n" + "="*70)
    print(" QC Classes in pyOpenMS - Usage Examples")
    print("="*70)
    print("\nThe following QC metric classes are now available in pyOpenMS:\n")
    
    qc_classes = [
        ("TIC", "Total Ion Current"),
        ("FWHM", "Full Width at Half Maximum"),
        ("PeptideMass", "Theoretical peptide mass annotation"),
        ("Contaminants", "Contaminant peptide detection"),
        ("Ms2IdentificationRate", "MS2 identification rate calculation"),
        ("MzCalibration", "m/z calibration error annotation"),
        ("RTAlignment", "Retention time alignment annotation"),
        ("MissedCleavages", "Missed cleavages counting"),
        ("FragmentMassError", "Fragment mass error calculation"),
        ("PSMExplainedIonCurrent", "PSM explained ion current"),
        ("Ms2SpectrumStats", "MS2 spectrum statistics"),
    ]
    
    for class_name, description in qc_classes:
        print(f"  • {class_name:25s} - {description}")
    
    print("\n" + "="*70)
    
    # Run examples
    try:
        example_tic_computation()
        example_peptide_mass()
        example_fwhm()
        example_ms2_identification_rate()
        example_rt_alignment()
        example_contaminants()
        
        print("\n" + "="*70)
        print(" All examples completed successfully!")
        print("="*70)
        
    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
