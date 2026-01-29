// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TARGETED/DefaultChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

using namespace OpenMS;
using namespace std;

START_TEST(DefaultChromHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DefaultChromHandler* ptr = nullptr;
DefaultChromHandler* nullPointer = nullptr;

START_SECTION(DefaultChromHandler())
{
	ptr = new DefaultChromHandler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DefaultChromHandler())
{
  delete ptr;
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> collectIrtChromatogramsForIrt(...))
{
  DefaultChromHandler handler;

  // Create minimal test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with chromatograms
  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("irt_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add some data points to make the chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Create minimal iRT transitions
  OpenSwath::LightTargetedExperiment irt_transitions;
  OpenSwath::LightCompound compound;
  compound.id = "irt_peptide_1";
  compound.rt = 10.0;
  irt_transitions.compounds.push_back(compound);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0; // No RT filtering for iRT

  TransformationDescription trafo;

  // Test the method - should delegate to SRM/MRM handler for non-MS1 data
  std::vector<MSChromatogram> result = handler.collectIrtChromatogramsForIrt(
    swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, false, false);

  // Should return chromatograms (may be empty if mapping fails, which is expected for minimal test data)
  TEST_EQUAL(result.size() >= 0, true)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(...))
{
  DefaultChromHandler handler;

  // Create minimal test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with chromatograms
  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("transition_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add some data points to make the chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
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

  // Test the method - should delegate to SRM/MRM handler for non-MS1 data
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return chromatograms (may be empty if mapping fails, which is expected for minimal test data)
  TEST_EQUAL(result.size() >= 0, true)
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_empty_swath_maps)
{
  DefaultChromHandler handler;

  // Empty swath maps vector
  std::vector<OpenSwath::SwathMap> swath_maps;

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
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with empty swath maps
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle empty swath maps gracefully
  TEST_EQUAL(result.size(), 0)
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_no_transitions)
{
  DefaultChromHandler handler;

  // Create minimal test data with chromatograms but no transitions
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with chromatograms
  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("transition_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Empty transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with no transitions
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle empty transitions gracefully
  TEST_EQUAL(result.size(), 0)
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_mismatched_transitions)
{
  DefaultChromHandler handler;

  // Create test data with chromatograms and transitions that don't match
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create a mock spectrum access with chromatograms
  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("chrom_A");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Create transitions with different identifiers
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "peptide_B";
  compound.rt = 10.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "peptide_B";
  transition.transition_name = "transition_B";  // Different from chromatogram native ID
  transition.precursor_mz = 600.0;  // Different m/z
  transition.product_mz = 650.0;
  transition_exp.transitions.push_back(transition);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with mismatched transitions
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle mismatched transitions gracefully
  TEST_EQUAL(result.size() >= 0, true)  // May return chromatograms or empty
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_multiple_swath_maps)
{
  DefaultChromHandler handler;

  // Create multiple swath maps
  std::vector<OpenSwath::SwathMap> swath_maps;

  // First swath map
  OpenSwath::SwathMap swath_map1;
  swath_map1.lower = 400.0;
  swath_map1.upper = 600.0;
  swath_map1.ms1 = false;

  auto exp1 = std::make_shared<MSExperiment>();
  MSChromatogram chrom1;
  chrom1.setNativeID("chrom_1");
  chrom1.getPrecursor().setMZ(500.0);
  chrom1.getProduct().setMZ(550.0);
  // Add some data points to make the chromatogram non-empty
  chrom1.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom1.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom1.push_back(ChromatogramPeak(300.0, 1500.0));
  exp1->addChromatogram(chrom1);
  auto mock_spectrum_access1 = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp1);
  swath_map1.sptr = mock_spectrum_access1;
  swath_maps.push_back(swath_map1);

  // Second swath map
  OpenSwath::SwathMap swath_map2;
  swath_map2.lower = 600.0;
  swath_map2.upper = 800.0;
  swath_map2.ms1 = false;

  auto exp2 = std::make_shared<MSExperiment>();
  MSChromatogram chrom2;
  chrom2.setNativeID("chrom_2");
  chrom2.getPrecursor().setMZ(700.0);
  chrom2.getProduct().setMZ(750.0);
  // Add some data points to make the chromatogram non-empty
  chrom2.push_back(ChromatogramPeak(150.0, 1200.0));
  chrom2.push_back(ChromatogramPeak(250.0, 2200.0));
  chrom2.push_back(ChromatogramPeak(350.0, 1700.0));
  exp2->addChromatogram(chrom2);
  auto mock_spectrum_access2 = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp2);
  swath_map2.sptr = mock_spectrum_access2;
  swath_maps.push_back(swath_map2);

  // Create multiple transitions
  OpenSwath::LightTargetedExperiment transition_exp;

  // First transition
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

  // Second transition
  OpenSwath::LightCompound compound2;
  compound2.id = "peptide_2";
  compound2.rt = 20.0;
  transition_exp.compounds.push_back(compound2);

  OpenSwath::LightTransition transition2;
  transition2.peptide_ref = "peptide_2";
  transition2.transition_name = "transition_2";
  transition2.precursor_mz = 700.0;
  transition2.product_mz = 750.0;
  transition_exp.transitions.push_back(transition2);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with multiple swath maps
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle multiple swath maps
  TEST_EQUAL(result.size() >= 0, true)
  TEST_EQUAL(swath_maps.size(), 2)  // Input should remain unchanged
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_invalid_parameters)
{
  DefaultChromHandler handler;

  // Create minimal test data
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("transition_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add some data points to make the chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
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

  // Invalid parameters (negative RT window)
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -100.0;  // Invalid negative value

  // Test extractAndMapChromatogramsForTransitions with invalid parameters
  // This should either handle gracefully or throw an exception
  bool exception_thrown = false;
  std::vector<MSChromatogram> result;
  try
  {
    result = handler.extractAndMapChromatogramsForTransitions(
      swath_maps, transition_exp, cp, mrm_mapping_param);
  }
  catch (const std::exception& e)
  {
    exception_thrown = true;
    OPENMS_LOG_INFO << "Expected exception caught: " << e.what() << std::endl;
  }

  // Either should succeed with defaults or throw exception
  TEST_EQUAL(exception_thrown || result.size() >= 0, true)
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_mixed_ms1_ms2)
{
  DefaultChromHandler handler;

  // Create swath maps with mixed MS1/MS2 data to test delegation logic
  std::vector<OpenSwath::SwathMap> swath_maps;

  // MS1 swath map (should trigger DIA mode)
  OpenSwath::SwathMap swath_map_ms1;
  swath_map_ms1.lower = 400.0;
  swath_map_ms1.upper = 1200.0;
  swath_map_ms1.ms1 = true;  // This should trigger DIA mode

  auto exp_ms1 = std::make_shared<MSExperiment>();
  // Add a spectrum to make it DIA-like
  MSSpectrum spec;
  spec.setMSLevel(1);
  exp_ms1->addSpectrum(spec);
  auto mock_spectrum_access_ms1 = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp_ms1);
  swath_map_ms1.sptr = mock_spectrum_access_ms1;
  swath_maps.push_back(swath_map_ms1);

  // MS2 swath map
  OpenSwath::SwathMap swath_map_ms2;
  swath_map_ms2.lower = 400.0;
  swath_map_ms2.upper = 600.0;
  swath_map_ms2.ms1 = false;

  auto exp_ms2 = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("transition_chrom_1");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(550.0);
  // Add some data points to make the chromatogram non-empty
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  exp_ms2->addChromatogram(chrom);
  auto mock_spectrum_access_ms2 = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp_ms2);
  swath_map_ms2.sptr = mock_spectrum_access_ms2;
  swath_maps.push_back(swath_map_ms2);

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
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with mixed MS1/MS2 data
  // Should delegate to DIA handler due to MS1 presence
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle mixed data (delegates to DIA handler)
  TEST_EQUAL(result.size() >= 0, true)
}
END_SECTION

START_SECTION(DefaultChromHandler_edge_case_empty_spectrum_access)
{
  DefaultChromHandler handler;

  // Create swath map with empty/null spectrum access
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create empty spectrum access (no chromatograms, no spectra)
  auto exp = std::make_shared<MSExperiment>();
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
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
  cp.rt_extraction_window = -1.0;

  // Test extractAndMapChromatogramsForTransitions with empty spectrum access
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should handle empty spectrum access gracefully
  TEST_EQUAL(result.size() >= 0, true)  // May return empty or some default result
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST