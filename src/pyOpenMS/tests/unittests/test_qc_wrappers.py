#!/usr/bin/env python
"""
Simple test to verify that QC classes are wrapped and accessible in pyOpenMS.

This test imports the QC classes and creates instances to verify they're properly wrapped.
"""

def test_qc_classes_import():
    """Test that QC classes can be imported"""
    try:
        import pyopenms as oms
        
        # Test QCBase nested classes
        print("Testing QCBase...")
        spec_map = oms.SpectraMap()
        print("  ✓ SpectraMap created")
        
        # Test TIC
        print("\nTesting TIC...")
        tic = oms.TIC()
        print("  ✓ TIC created")
        print(f"  ✓ TIC name: {tic.getName()}")
        
        # Test FWHM
        print("\nTesting FWHM...")
        fwhm = oms.FWHM()
        print("  ✓ FWHM created")
        print(f"  ✓ FWHM name: {fwhm.getName()}")
        
        # Test PeptideMass
        print("\nTesting PeptideMass...")
        peptide_mass = oms.PeptideMass()
        print("  ✓ PeptideMass created")
        print(f"  ✓ PeptideMass name: {peptide_mass.getName()}")
        
        # Test Contaminants
        print("\nTesting Contaminants...")
        contaminants = oms.Contaminants()
        print("  ✓ Contaminants created")
        print(f"  ✓ Contaminants name: {contaminants.getName()}")
        
        # Test Ms2IdentificationRate
        print("\nTesting Ms2IdentificationRate...")
        ms2_id_rate = oms.Ms2IdentificationRate()
        print("  ✓ Ms2IdentificationRate created")
        print(f"  ✓ Ms2IdentificationRate name: {ms2_id_rate.getName()}")
        
        # Test MzCalibration
        print("\nTesting MzCalibration...")
        mz_calibration = oms.MzCalibration()
        print("  ✓ MzCalibration created")
        print(f"  ✓ MzCalibration name: {mz_calibration.getName()}")
        
        # Test RTAlignment
        print("\nTesting RTAlignment...")
        rt_alignment = oms.RTAlignment()
        print("  ✓ RTAlignment created")
        print(f"  ✓ RTAlignment name: {rt_alignment.getName()}")
        
        # Test MissedCleavages
        print("\nTesting MissedCleavages...")
        missed_cleavages = oms.MissedCleavages()
        print("  ✓ MissedCleavages created")
        print(f"  ✓ MissedCleavages name: {missed_cleavages.getName()}")
        
        # Test FragmentMassError
        print("\nTesting FragmentMassError...")
        fme = oms.FragmentMassError()
        print("  ✓ FragmentMassError created")
        print(f"  ✓ FragmentMassError name: {fme.getName()}")
        
        # Test PSMExplainedIonCurrent
        print("\nTesting PSMExplainedIonCurrent...")
        psm = oms.PSMExplainedIonCurrent()
        print("  ✓ PSMExplainedIonCurrent created")
        print(f"  ✓ PSMExplainedIonCurrent name: {psm.getName()}")
        
        # Test Ms2SpectrumStats
        print("\nTesting Ms2SpectrumStats...")
        ms2_stats = oms.Ms2SpectrumStats()
        print("  ✓ Ms2SpectrumStats created")
        print(f"  ✓ Ms2SpectrumStats name: {ms2_stats.getName()}")
        
        print("\n" + "="*60)
        print("All QC classes successfully imported and instantiated!")
        print("="*60)
        
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("\nNote: This test requires pyOpenMS to be built with the new QC wrappers.")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_qc_classes_import()
    exit(0 if success else 1)
