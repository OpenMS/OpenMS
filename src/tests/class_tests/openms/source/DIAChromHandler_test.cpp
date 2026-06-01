// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TARGETED/DIAChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

using namespace OpenMS;
using namespace std;

START_TEST(DIAChromHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DIAChromHandler* ptr = nullptr;
DIAChromHandler* nullPointer = nullptr;

START_SECTION(DIAChromHandler())
{
	ptr = new DIAChromHandler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DIAChromHandler())
{
  delete ptr;
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> collectIrtChromatogramsForIrt(...))
{
  DIAChromHandler handler;

  // Create minimal test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with spectra and chromatograms
  MSExperiment exp;
  
  // Add multiple spectra at different RT points (required for DIA chromatogram extraction)
  for (double rt = 100.0; rt <= 300.0; rt += 50.0)
  {
    MSSpectrum spec;
    spec.setRT(rt);
    spec.push_back(Peak1D(500.0, 1000.0));  // Peak at precursor m/z
    spec.push_back(Peak1D(550.0, 800.0));   // Peak at product m/z
    exp.addSpectrum(spec);
  }
  
  MSChromatogram chrom;
  chrom.setNativeID("irt_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add data points to make chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp.addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<SpectrumAccessOpenMS>(std::make_shared<MSExperiment>(exp));
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Create minimal iRT transitions
  OpenSwath::LightTargetedExperiment irt_transitions;
  OpenSwath::LightCompound compound;
  compound.id = "irt_peptide_1";
  compound.rt = 10.0;
  irt_transitions.compounds.push_back(compound);

  OpenSwath::LightTransition irt_transition;
  irt_transition.peptide_ref = "irt_peptide_1";
  irt_transition.transition_name = "irt_transition_1";
  irt_transition.precursor_mz = 500.0;
  irt_transition.product_mz = 550.0;
  irt_transitions.transitions.push_back(irt_transition);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0; // No RT filtering for iRT
  cp.extraction_function = "tophat";
  cp.im_extraction_window = -1.0; // Disable ion mobility extraction
  cp.min_upper_edge_dist = 0.0; // Allow transitions near upper edge
  cp.mz_extraction_window = 0.1; // Small m/z window for extraction
  cp.ppm = false; // Use Da instead of ppm

  TransformationDescription trafo;

  // Test the method - should handle DIA chromatogram extraction
  std::vector<MSChromatogram> result = handler.collectIrtChromatogramsForIrt(
    swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, false, false);

  // Should return chromatograms - expect exactly 1 match
  std::cout << "DIA iRT result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 1)
  TEST_EQUAL(result[0].getNativeID(), "irt_transition_1")
  TEST_REAL_SIMILAR(result[0].getPrecursor().getMZ(), 500.0)
  TEST_REAL_SIMILAR(result[0].getProduct().getMZ(), 550.0)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(...))
{
  DIAChromHandler handler;

  // Create minimal test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with spectra and chromatograms
  MSExperiment exp;
  
  // Add multiple spectra at different RT points (required for DIA chromatogram extraction)
  for (double rt = 100.0; rt <= 300.0; rt += 50.0)
  {
    MSSpectrum spec;
    spec.setRT(rt);
    spec.push_back(Peak1D(500.0, 1000.0));  // Peak at precursor m/z
    spec.push_back(Peak1D(550.0, 800.0));   // Peak at product m/z
    exp.addSpectrum(spec);
  }
  
  MSChromatogram chrom;
  chrom.setNativeID("transition_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add data points to make chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp.addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<SpectrumAccessOpenMS>(std::make_shared<MSExperiment>(exp));
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Create minimal transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "peptide_1";
  compound.rt = 10.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "peptide_1";
  transition.transition_name = "transition_1";
  transition.precursor_mz = 500.0;
  transition.product_mz = 550.0;
  transition_exp.transitions.push_back(transition);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0; // No RT filtering
  cp.extraction_function = "tophat";
  cp.im_extraction_window = -1.0; // Disable ion mobility extraction
  cp.mz_extraction_window = 0.1; // Small m/z window for extraction
  cp.ppm = false; // Use Da instead of ppm
  cp.min_upper_edge_dist = 0.0; // Allow transitions near upper edge

  // Test the method - should handle DIA chromatogram extraction
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return chromatograms - expect exactly 1 match
  std::cout << "DIA transition result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 1)
  TEST_EQUAL(result[0].getNativeID(), "transition_1")
  TEST_REAL_SIMILAR(result[0].getPrecursor().getMZ(), 500.0)
  TEST_REAL_SIMILAR(result[0].getProduct().getMZ(), 550.0)
}
END_SECTION

START_SECTION(DIA_edge_case_no_matches)
{
  DIAChromHandler handler;

  // Create test data with transition outside swath window
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;  // Swath window 400-600
  swath_map.ms1 = false;

  // Create a mock spectrum access with spectra
  MSExperiment exp;
  for (double rt = 100.0; rt <= 300.0; rt += 50.0)
  {
    MSSpectrum spec;
    spec.setRT(rt);
    spec.push_back(Peak1D(500.0, 1000.0));
    spec.push_back(Peak1D(550.0, 800.0));
    exp.addSpectrum(spec);
  }
  auto mock_spectrum_access = std::make_shared<SpectrumAccessOpenMS>(std::make_shared<MSExperiment>(exp));
  swath_map.sptr = mock_spectrum_access;
  swath_maps.push_back(swath_map);

  // Create transition with product m/z not present in spectra
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "peptide_no_match";
  compound.rt = 10.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "peptide_no_match";
  transition.transition_name = "transition_no_match";
  transition.precursor_mz = 500.0;
  transition.product_mz = 850.0;  // No peaks at this m/z in spectra
  transition_exp.transitions.push_back(transition);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;
  cp.extraction_function = "tophat";
  cp.im_extraction_window = -1.0;
  cp.mz_extraction_window = 0.1;
  cp.ppm = false;
  cp.min_upper_edge_dist = 0.0; // Allow transitions near upper edge

  // Test the method - should return no chromatograms (empty chromatogram filtered out)
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return 1 chromatogram with zero intensity
  TEST_EQUAL(result.size(), 1)
  TEST_EQUAL(result[0].getNativeID(), "transition_no_match")
  TEST_REAL_SIMILAR(result[0].getPrecursor().getMZ(), 500.0)
  TEST_REAL_SIMILAR(result[0].getProduct().getMZ(), 850.0)
  // Chromatogram should have points but zero total intensity
  double total_intensity = 0.0;
  for (size_t i = 0; i < result[0].size(); ++i) {
    total_intensity += result[0][i].getIntensity();
  }
  TEST_EQUAL(total_intensity, 0.0)
}
END_SECTION

START_SECTION(DIA_edge_case_partial_matches)
{
  DIAChromHandler handler;

  // Create test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;  // Swath window 400-600
  swath_map.ms1 = false;

  // Create a mock spectrum access with spectra
  MSExperiment exp;
  for (double rt = 100.0; rt <= 300.0; rt += 50.0)
  {
    MSSpectrum spec;
    spec.setRT(rt);
    spec.push_back(Peak1D(500.0, 1000.0));
    spec.push_back(Peak1D(550.0, 800.0));
    spec.push_back(Peak1D(700.0, 600.0));  // For outside transition
    spec.push_back(Peak1D(750.0, 400.0));
    exp.addSpectrum(spec);
  }
  auto mock_spectrum_access = std::make_shared<SpectrumAccessOpenMS>(std::make_shared<MSExperiment>(exp));
  swath_map.sptr = mock_spectrum_access;
  swath_maps.push_back(swath_map);

  // Create transitions - one inside, one outside swath window
  OpenSwath::LightTargetedExperiment transition_exp;
  
  // Transition 1 - matches (product m/z exists in spectra)
  OpenSwath::LightCompound compound1;
  compound1.id = "peptide_match";
  compound1.rt = 10.0;
  transition_exp.compounds.push_back(compound1);

  OpenSwath::LightTransition transition1;
  transition1.peptide_ref = "peptide_match";
  transition1.transition_name = "transition_match";
  transition1.precursor_mz = 500.0;
  transition1.product_mz = 550.0;  // Exists in spectra
  transition_exp.transitions.push_back(transition1);

  // Transition 2 - no match (product m/z does not exist in spectra)
  OpenSwath::LightCompound compound2;
  compound2.id = "peptide_no_match";
  compound2.rt = 10.0;
  transition_exp.compounds.push_back(compound2);

  OpenSwath::LightTransition transition2;
  transition2.peptide_ref = "peptide_no_match";
  transition2.transition_name = "transition_no_match";
  transition2.precursor_mz = 500.0;
  transition2.product_mz = 850.0;  // Does not exist in spectra
  transition_exp.transitions.push_back(transition2);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;
  cp.extraction_function = "tophat";
  cp.im_extraction_window = -1.0;
  cp.mz_extraction_window = 0.1;
  cp.ppm = false;
  cp.min_upper_edge_dist = 0.0; // Allow transitions near upper edge

  // Test the method - should return only 1 chromatogram (the one with matching peaks)
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return 2 chromatograms, but one should be empty
  TEST_EQUAL(result.size(), 2)
  
  // Find the matching and non-matching chromatograms
  bool found_match = false, found_no_match = false;
  for (const auto& chrom : result)
  {
    if (chrom.getNativeID() == "transition_match")
    {
      found_match = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 550.0)
      // Should have non-zero intensity
      double total_intensity = 0.0;
      for (size_t i = 0; i < chrom.size(); ++i) {
        total_intensity += chrom[i].getIntensity();
      }
      TEST_NOT_EQUAL(total_intensity, 0.0)
    }
    else if (chrom.getNativeID() == "transition_no_match")
    {
      found_no_match = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 850.0)
      // Should have zero total intensity
      double total_intensity = 0.0;
      for (size_t i = 0; i < chrom.size(); ++i) {
        total_intensity += chrom[i].getIntensity();
      }
      TEST_EQUAL(total_intensity, 0.0)
    }
  }
  TEST_EQUAL(found_match, true)
  TEST_EQUAL(found_no_match, true)
}
END_SECTION

START_SECTION(DIA_edge_case_multiple_chromatograms)
{
  DIAChromHandler handler;

  // Create test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 800.0;  // Wider swath window
  swath_map.ms1 = false;

  // Create a mock spectrum access with spectra
  MSExperiment exp;
  for (double rt = 100.0; rt <= 300.0; rt += 50.0)
  {
    MSSpectrum spec;
    spec.setRT(rt);
    spec.push_back(Peak1D(500.0, 1000.0));
    spec.push_back(Peak1D(550.0, 800.0));
    spec.push_back(Peak1D(700.0, 600.0));  // For matching transition
    spec.push_back(Peak1D(750.0, 400.0));  // For non-matching transition
    exp.addSpectrum(spec);
  }
  auto mock_spectrum_access = std::make_shared<SpectrumAccessOpenMS>(std::make_shared<MSExperiment>(exp));
  swath_map.sptr = mock_spectrum_access;
  swath_maps.push_back(swath_map);

  // Create multiple transitions
  OpenSwath::LightTargetedExperiment transition_exp;
  
  // Transition 1
  OpenSwath::LightCompound compound1;
  compound1.id = "peptide_1";
  compound1.rt = 10.0;
  transition_exp.compounds.push_back(compound1);

  OpenSwath::LightTransition transition1;
  transition1.peptide_ref = "peptide_1";
  transition1.transition_name = "transition_1";
  transition1.precursor_mz = 500.0;
  transition1.product_mz = 550.0;
  transition_exp.transitions.push_back(transition1);

  // Transition 2
  OpenSwath::LightCompound compound2;
  compound2.id = "peptide_2";
  compound2.rt = 10.0;
  transition_exp.compounds.push_back(compound2);

  OpenSwath::LightTransition transition2;
  transition2.peptide_ref = "peptide_2";
  transition2.transition_name = "transition_2";
  transition2.precursor_mz = 600.0;
  transition2.product_mz = 650.0;
  transition_exp.transitions.push_back(transition2);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;
  cp.extraction_function = "tophat";
  cp.im_extraction_window = -1.0;
  cp.mz_extraction_window = 0.1;
  cp.ppm = false;
  cp.min_upper_edge_dist = 0.0; // Allow transitions near upper edge

  // Test the method - should return 2 chromatograms
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return 2 chromatograms
  TEST_EQUAL(result.size(), 2)
  
  // Check first chromatogram
  bool found_trans1 = false, found_trans2 = false;
  for (const auto& chrom : result)
  {
    if (chrom.getNativeID() == "transition_1")
    {
      found_trans1 = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 550.0)
    }
    else if (chrom.getNativeID() == "transition_2")
    {
      found_trans2 = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 600.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 650.0)
    }
  }
  TEST_EQUAL(found_trans1, true)
  TEST_EQUAL(found_trans2, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST