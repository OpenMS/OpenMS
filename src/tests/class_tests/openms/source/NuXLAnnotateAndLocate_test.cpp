// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotateAndLocate.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace std;
using namespace OpenMS;

START_TEST(NuXLAnnotateAndLocate, "$Id$")

/////////////////////////////////////////////////////////////

// Note: NuXLAnnotateAndLocate only contains static methods, so no constructor/destructor tests needed

START_SECTION((static void annotateAndLocate_(const PeakMap& exp, std::vector<std::vector<NuXLAnnotatedHit>>& annotated_hits, const NuXLModificationMassesResult& mm, const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, Size max_variable_mods_per_peptide, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const NuXLParameterParsing::PrecursorsToMS2Adducts& all_feasible_adducts)))
{
  // Create minimal test data for the annotateAndLocate_ method
  PeakMap exp;
  
  // Create a simple spectrum
  MSSpectrum spectrum;
  spectrum.setRT(100.0);
  spectrum.setMSLevel(2);
  
  // Add some peaks
  Peak1D peak;
  peak.setMZ(100.0);
  peak.setIntensity(1000.0);
  spectrum.push_back(peak);
  
  peak.setMZ(200.0);
  peak.setIntensity(500.0);
  spectrum.push_back(peak);
  
  // Add precursor information
  Precursor precursor;
  precursor.setMZ(300.0);
  precursor.setCharge(2);
  spectrum.getPrecursors().push_back(precursor);
  
  // Add required data arrays
  spectrum.getIntegerDataArrays().resize(1);
  spectrum.getIntegerDataArrays()[0].push_back(1); // charge for first peak
  spectrum.getIntegerDataArrays()[0].push_back(1); // charge for second peak
  
  exp.addSpectrum(spectrum);
  
  // Create empty annotated hits vector
  vector<vector<NuXLAnnotatedHit>> annotated_hits;
  annotated_hits.resize(1); // One spectrum
  
  // Create minimal NuXLModificationMassesResult
  NuXLModificationMassesResult mm;
  mm.formula2mass[""] = 0.0;
  mm.mod_combinations[""].insert("none");
  
  // Create empty modification maps
  ModifiedPeptideGenerator::MapToResidueType fixed_modifications;
  ModifiedPeptideGenerator::MapToResidueType variable_modifications;
  
  // Create empty feasible adducts
  NuXLParameterParsing::PrecursorsToMS2Adducts all_feasible_adducts;
  
  // Test parameters
  Size max_variable_mods_per_peptide = 2;
  double fragment_mass_tolerance = 0.1;
  bool fragment_mass_tolerance_unit_ppm = false;
  
  // Test with empty annotated hits (should not crash)
  TEST_EXCEPTION(Exception::BaseException, 
    NuXLAnnotateAndLocate::annotateAndLocate_(
      exp, 
      annotated_hits, 
      mm, 
      fixed_modifications, 
      variable_modifications, 
      max_variable_mods_per_peptide, 
      fragment_mass_tolerance, 
      fragment_mass_tolerance_unit_ppm, 
      all_feasible_adducts
    )
  )
  
  // Since this is a complex integration method that requires extensive setup,
  // we mainly test that it doesn't crash with minimal valid input.
  // More comprehensive tests would require setting up proper NuXLAnnotatedHit objects
  // with valid sequences, modification indices, etc.
  
  // Test that the method exists and can be called (basic smoke test)
  // The method should handle empty annotated_hits gracefully
  annotated_hits[0].clear(); // Ensure it's empty
  
  // This should not crash and should return without doing anything
  NuXLAnnotateAndLocate::annotateAndLocate_(
    exp, 
    annotated_hits, 
    mm, 
    fixed_modifications, 
    variable_modifications, 
    max_variable_mods_per_peptide, 
    fragment_mass_tolerance, 
    fragment_mass_tolerance_unit_ppm, 
    all_feasible_adducts
  );
  
  // Verify that empty input results in no changes
  TEST_EQUAL(annotated_hits[0].size(), 0)
}
END_SECTION

END_TEST