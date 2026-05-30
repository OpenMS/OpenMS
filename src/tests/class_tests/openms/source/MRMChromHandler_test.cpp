// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TARGETED/MRMChromHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMChromHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMChromHandler* ptr = nullptr;
MRMChromHandler* nullPointer = nullptr;

START_SECTION(MRMChromHandler())
{
	ptr = new MRMChromHandler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMChromHandler())
{
  delete ptr;
}
END_SECTION

START_SECTION(static void normalizeChromatogramMZ(MSChromatogram& chrom))
{
  // Test with no m/z values set
  {
    MSChromatogram chrom;
    chrom.getPrecursor().setMZ(0.0);
    chrom.getProduct().setMZ(0.0);
    chrom.setNativeID("test");

    MRMChromHandler::normalizeChromatogramMZ(chrom);
    // Meta values are only set for positive m/z values
    TEST_EQUAL(chrom.getMetaValue("precursor_mz").isEmpty(), true)
    TEST_EQUAL(chrom.getMetaValue("product_mz").isEmpty(), true)
    // Native ID should remain unchanged since no valid m/z values
    TEST_EQUAL(chrom.getNativeID(), "test")
  }

  // Test with CV terms (MS:1000827)
  {
    MSChromatogram chrom;
    chrom.getPrecursor().addCVTerm(CVTerm("MS:1000827", "", "", "500.123"));
    chrom.getProduct().addCVTerm(CVTerm("MS:1000827", "", "", "600.456"));
    chrom.setNativeID("original_id");

    MRMChromHandler::normalizeChromatogramMZ(chrom);
    TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.123)
    TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 600.456)
    TEST_REAL_SIMILAR(chrom.getMetaValue("precursor_mz"), 500.123)
    TEST_REAL_SIMILAR(chrom.getMetaValue("product_mz"), 600.456)
    TEST_EQUAL(chrom.getNativeID().hasPrefix("precursor="), true)
    TEST_EQUAL(chrom.getNativeID().hasSubstring(",product="), true)
  }

  // Test parsing from nativeID
  {
    MSChromatogram chrom;
    chrom.getPrecursor().setMZ(0.0);
    chrom.getProduct().setMZ(0.0);
    // Use CV terms instead of nativeID parsing for proper SRM/MRM mzML data
    chrom.getPrecursor().addCVTerm(CVTerm("MS:1000827", "", "", "700.789"));
    chrom.getProduct().addCVTerm(CVTerm("MS:1000827", "", "", "800.012"));

    MRMChromHandler::normalizeChromatogramMZ(chrom);
    TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 700.789)
    TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 800.012)
    TEST_REAL_SIMILAR(chrom.getMetaValue("precursor_mz"), 700.789)
    TEST_REAL_SIMILAR(chrom.getMetaValue("product_mz"), 800.012)
  }
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> collectIrtChromatogramsForIrt(...))
{
  MRMChromHandler handler;

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

  TransformationDescription trafo;

  // Test the method
  std::vector<MSChromatogram> result = handler.collectIrtChromatogramsForIrt(
    swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, false, false);

  // Should return chromatograms - now we expect exactly 1 match
  std::cout << "Result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 1)
  TEST_EQUAL(result[0].getNativeID(), "irt_transition_1")  // MRMMapping sets nativeID to transition name
  TEST_REAL_SIMILAR(result[0].getPrecursor().getMZ(), 500.0)
  TEST_REAL_SIMILAR(result[0].getProduct().getMZ(), 550.0)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions(...))
{
  MRMChromHandler handler;

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

  // Test the method
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return chromatograms - now we expect exactly 1 match
  std::cout << "Result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 1)
  TEST_EQUAL(result[0].getNativeID(), "transition_1")  // MRMMapping sets nativeID to transition name
  TEST_REAL_SIMILAR(result[0].getPrecursor().getMZ(), 500.0)
  TEST_REAL_SIMILAR(result[0].getProduct().getMZ(), 550.0)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> collectIrtChromatogramsForIrt - No matches)
{
  MRMChromHandler handler;

  // Create test data where chromatograms don't match transitions
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 600.0;
  swath_map.ms1 = false;

  // Create chromatogram with m/z values far from transition
  auto exp = std::make_shared<MSExperiment>();
  MSChromatogram chrom;
  chrom.setNativeID("no_match_chrom");
  chrom.getPrecursor().setMZ(600.0);  // Different from transition
  chrom.getProduct().setMZ(650.0);    // Different from transition
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  exp->addChromatogram(chrom);
  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;

  swath_maps.push_back(swath_map);

  // Create iRT transitions that won't match
  OpenSwath::LightTargetedExperiment irt_transitions;
  OpenSwath::LightCompound compound;
  compound.id = "irt_peptide_no_match";
  compound.rt = 10.0;
  irt_transitions.compounds.push_back(compound);

  OpenSwath::LightTransition irt_transition;
  irt_transition.peptide_ref = "irt_peptide_no_match";
  irt_transition.transition_name = "irt_transition_no_match";
  irt_transition.precursor_mz = 500.0;  // Different from chromatogram
  irt_transition.product_mz = 550.0;    // Different from chromatogram
  irt_transitions.transitions.push_back(irt_transition);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  TransformationDescription trafo;

  // Test the method - should return empty result
  std::vector<MSChromatogram> result = handler.collectIrtChromatogramsForIrt(
    swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, false, false);

  // Should return no chromatograms when no matches found
  std::cout << "No matches result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 0)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> extractAndMapChromatogramsForTransitions - Partial matches)
{
  MRMChromHandler handler;

  // Create test data with 3 chromatograms and 3 transitions, where only 2 match
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 800.0;  // Wider range to include all chromatograms
  swath_map.ms1 = false;

  auto exp = std::make_shared<MSExperiment>();

  // Chromatogram 1 - will match transition 1
  MSChromatogram chrom1;
  chrom1.setNativeID("match_chrom_1");
  chrom1.getPrecursor().setMZ(500.0);
  chrom1.getProduct().setMZ(550.0);
  chrom1.push_back(ChromatogramPeak(100.0, 1000.0));
  exp->addChromatogram(chrom1);

  // Chromatogram 2 - will match transition 2
  MSChromatogram chrom2;
  chrom2.setNativeID("match_chrom_2");
  chrom2.getPrecursor().setMZ(600.0);
  chrom2.getProduct().setMZ(650.0);
  chrom2.push_back(ChromatogramPeak(200.0, 2000.0));
  exp->addChromatogram(chrom2);

  // Chromatogram 3 - will NOT match any transition
  MSChromatogram chrom3;
  chrom3.setNativeID("no_match_chrom_3");
  chrom3.getPrecursor().setMZ(700.0);
  chrom3.getProduct().setMZ(750.0);
  chrom3.push_back(ChromatogramPeak(300.0, 3000.0));
  exp->addChromatogram(chrom3);

  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;
  swath_maps.push_back(swath_map);

  // Create transition experiment with 3 transitions
  OpenSwath::LightTargetedExperiment transition_exp;

  // Compound 1
  OpenSwath::LightCompound compound1;
  compound1.id = "peptide_1";
  compound1.rt = 10.0;
  transition_exp.compounds.push_back(compound1);

  // Transition 1 - matches chrom1
  OpenSwath::LightTransition transition1;
  transition1.peptide_ref = "peptide_1";
  transition1.transition_name = "transition_1";
  transition1.precursor_mz = 500.0;
  transition1.product_mz = 550.0;
  transition_exp.transitions.push_back(transition1);

  // Compound 2
  OpenSwath::LightCompound compound2;
  compound2.id = "peptide_2";
  compound2.rt = 20.0;
  transition_exp.compounds.push_back(compound2);

  // Transition 2 - matches chrom2
  OpenSwath::LightTransition transition2;
  transition2.peptide_ref = "peptide_2";
  transition2.transition_name = "transition_2";
  transition2.precursor_mz = 600.0;
  transition2.product_mz = 650.0;
  transition_exp.transitions.push_back(transition2);

  // Compound 3
  OpenSwath::LightCompound compound3;
  compound3.id = "peptide_3";
  compound3.rt = 30.0;
  transition_exp.compounds.push_back(compound3);

  // Transition 3 - does NOT match any chromatogram
  OpenSwath::LightTransition transition3;
  transition3.peptide_ref = "peptide_3";
  transition3.transition_name = "transition_3";
  transition3.precursor_mz = 800.0;  // No matching chromatogram
  transition3.product_mz = 850.0;    // No matching chromatogram
  transition_exp.transitions.push_back(transition3);

  // Test parameters
  Param mrm_mapping_param;
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  // Test the method - should return 2 chromatograms (partial matches)
  std::vector<MSChromatogram> result = handler.extractAndMapChromatogramsForTransitions(
    swath_maps, transition_exp, cp, mrm_mapping_param);

  // Should return 2 chromatograms (2 out of 3 transitions matched)
  std::cout << "Partial matches result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 2)

  // Verify the matched chromatograms
  bool found_transition_1 = false;
  bool found_transition_2 = false;
  for (const auto& chrom : result) {
    if (chrom.getNativeID() == "transition_1") {
      found_transition_1 = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 550.0)
    } else if (chrom.getNativeID() == "transition_2") {
      found_transition_2 = true;
      TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 600.0)
      TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 650.0)
    }
  }
  TEST_EQUAL(found_transition_1, true)
  TEST_EQUAL(found_transition_2, true)
}
END_SECTION

START_SECTION(std::vector<MSChromatogram> collectIrtChromatogramsForIrt - Multiple chromatograms)
{
  MRMChromHandler handler;

  // Create test data with multiple chromatograms for one iRT transition
  std::vector<OpenSwath::SwathMap> swath_maps;
  OpenSwath::SwathMap swath_map;
  swath_map.lower = 400.0;
  swath_map.upper = 700.0;
  swath_map.ms1 = false;

  auto exp = std::make_shared<MSExperiment>();

  // Multiple chromatograms with same m/z (simulating multiple SWATH windows or replicates)
  MSChromatogram chrom1;
  chrom1.setNativeID("irt_chrom_1");
  chrom1.getPrecursor().setMZ(500.0);
  chrom1.getProduct().setMZ(550.0);
  chrom1.push_back(ChromatogramPeak(100.0, 1000.0));
  exp->addChromatogram(chrom1);

  MSChromatogram chrom2;
  chrom2.setNativeID("irt_chrom_2");
  chrom2.getPrecursor().setMZ(500.0);  // Same m/z as chrom1
  chrom2.getProduct().setMZ(550.0);    // Same m/z as chrom1
  chrom2.push_back(ChromatogramPeak(150.0, 1500.0));
  exp->addChromatogram(chrom2);

  auto mock_spectrum_access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(exp);
  swath_map.sptr = mock_spectrum_access;
  swath_maps.push_back(swath_map);

  // Create iRT transitions
  OpenSwath::LightTargetedExperiment irt_transitions;
  OpenSwath::LightCompound compound;
  compound.id = "irt_peptide_multi";
  compound.rt = 10.0;
  irt_transitions.compounds.push_back(compound);

  OpenSwath::LightTransition irt_transition;
  irt_transition.peptide_ref = "irt_peptide_multi";
  irt_transition.transition_name = "irt_transition_multi";
  irt_transition.precursor_mz = 500.0;
  irt_transition.product_mz = 550.0;
  irt_transitions.transitions.push_back(irt_transition);

  // Test parameters
  Param mrm_mapping_param;
  // Enable multiple assays to allow multiple chromatograms per transition
  mrm_mapping_param.setValue("map_multiple_assays", "true");
  ChromExtractParams cp;
  cp.rt_extraction_window = -1.0;

  TransformationDescription trafo;

  // Test the method - should return multiple chromatograms for one transition
  std::vector<MSChromatogram> result = handler.collectIrtChromatogramsForIrt(
    swath_maps, irt_transitions, mrm_mapping_param, cp, trafo, false, false);

  // Should return 2 chromatograms (multiple matches for same transition)
  std::cout << "Multiple chromatograms result size: " << result.size() << std::endl;
  TEST_EQUAL(result.size(), 2)

  // Both should have the same transition name but different original nativeIDs stored in meta
  for (const auto& chrom : result) {
    TEST_EQUAL(chrom.getNativeID(), "irt_transition_multi")
    TEST_REAL_SIMILAR(chrom.getPrecursor().getMZ(), 500.0)
    TEST_REAL_SIMILAR(chrom.getProduct().getMZ(), 550.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST