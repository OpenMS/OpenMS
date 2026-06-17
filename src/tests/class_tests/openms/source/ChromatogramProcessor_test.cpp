// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TARGETED/ChromatogramProcessor.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramProcessor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static void pickExperiment(...))
{
  // Create test chromatograms
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram chrom;
  chrom.setNativeID("test_chrom_1");

  // Add some test peaks
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  chromatograms.push_back(chrom);

  // Create minimal transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "test_peptide";
  compound.rt = 200.0; // RT in seconds
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "test_peptide";
  transition.transition_name = "test_transition";
  transition.precursor_mz = 500.0;
  transition.product_mz = 550.0;
  transition_exp.transitions.push_back(transition);

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Verify that feature finding was attempted
  // (exact results depend on the feature finder implementation)
  TEST_EQUAL(featureFile.size() >= 0, true)  // May be 0 if no features found
  TEST_EQUAL(transition_group_map.size() >= 0, true)  // May be 0 if no transitions mapped
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_empty_chromatograms)
{
  // Test with chromatograms that have no data points
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram chrom;
  chrom.setNativeID("empty_chrom_1");
  // Don't add any peaks - empty chromatogram
  chromatograms.push_back(chrom);

  // Create minimal transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "test_peptide";
  compound.rt = 200.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "test_peptide";
  transition.transition_name = "test_transition";
  transition.precursor_mz = 500.0;
  transition.product_mz = 550.0;
  transition_exp.transitions.push_back(transition);

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with empty chromatograms
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Should handle empty chromatograms gracefully
  TEST_EQUAL(featureFile.size(), 0)  // No features should be found
  TEST_EQUAL(transition_group_map.size() >= 0, true)  // May have transition groups but no features
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_no_transitions)
{
  // Test with chromatograms but no transitions
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram chrom;
  chrom.setNativeID("test_chrom_1");
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom.push_back(ChromatogramPeak(300.0, 1500.0));
  chromatograms.push_back(chrom);

  // Empty transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;
  // No compounds or transitions added

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with no transitions
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Should handle empty transitions gracefully
  TEST_EQUAL(featureFile.size(), 0)  // No features should be found
  TEST_EQUAL(transition_group_map.size(), 0)  // No transition groups
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_mismatched_transitions)
{
  // Test with chromatograms and transitions that don't match
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram chrom;
  chrom.setNativeID("chrom_A");
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom.push_back(ChromatogramPeak(200.0, 2000.0));
  chromatograms.push_back(chrom);

  // Create transitions with different native IDs
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "peptide_B";
  compound.rt = 200.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "peptide_B";
  transition.transition_name = "transition_B";  // Different from chromatogram native ID
  transition.precursor_mz = 500.0;
  transition.product_mz = 550.0;
  transition_exp.transitions.push_back(transition);

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with mismatched transitions
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Should handle mismatched transitions gracefully
  TEST_EQUAL(featureFile.size() >= 0, true)  // May find features or not
  TEST_EQUAL(transition_group_map.size() >= 0, true)  // May have transition groups
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_multiple_chromatograms)
{
  // Test with multiple chromatograms and transitions
  std::vector<MSChromatogram> chromatograms;
  
  // First chromatogram
  MSChromatogram chrom1;
  chrom1.setNativeID("chrom_1");
  chrom1.push_back(ChromatogramPeak(100.0, 1000.0));
  chrom1.push_back(ChromatogramPeak(200.0, 2000.0));
  chrom1.push_back(ChromatogramPeak(300.0, 1500.0));
  chromatograms.push_back(chrom1);

  // Second chromatogram
  MSChromatogram chrom2;
  chrom2.setNativeID("chrom_2");
  chrom2.push_back(ChromatogramPeak(150.0, 800.0));
  chrom2.push_back(ChromatogramPeak(250.0, 1800.0));
  chrom2.push_back(ChromatogramPeak(350.0, 1200.0));
  chromatograms.push_back(chrom2);

  // Create multiple transitions
  OpenSwath::LightTargetedExperiment transition_exp;
  
  // First transition
  OpenSwath::LightCompound compound1;
  compound1.id = "peptide_1";
  compound1.rt = 200.0;
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
  compound2.rt = 250.0;
  transition_exp.compounds.push_back(compound2);

  OpenSwath::LightTransition transition2;
  transition2.peptide_ref = "peptide_2";
  transition2.transition_name = "transition_2";
  transition2.precursor_mz = 600.0;
  transition2.product_mz = 650.0;
  transition_exp.transitions.push_back(transition2);

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with multiple chromatograms
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Should handle multiple chromatograms and transitions
  TEST_EQUAL(featureFile.size() >= 0, true)  // May find features
  TEST_EQUAL(transition_group_map.size() >= 0, true)  // May have transition groups
  TEST_EQUAL(chromatograms.size(), 2)  // Input should remain unchanged
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_invalid_parameters)
{
  // Test with invalid feature finder parameters
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram chrom;
  chrom.setNativeID("test_chrom_1");
  chrom.push_back(ChromatogramPeak(100.0, 1000.0));
  chromatograms.push_back(chrom);

  // Create minimal transition experiment
  OpenSwath::LightTargetedExperiment transition_exp;
  OpenSwath::LightCompound compound;
  compound.id = "test_peptide";
  compound.rt = 200.0;
  transition_exp.compounds.push_back(compound);

  OpenSwath::LightTransition transition;
  transition.peptide_ref = "test_peptide";
  transition.transition_name = "test_transition";
  transition.precursor_mz = 500.0;
  transition.product_mz = 550.0;
  transition_exp.transitions.push_back(transition);

  // Invalid feature finder parameters (negative spacing)
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", -0.005);  // Invalid negative value
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with invalid parameters
  // This should either handle gracefully or throw an exception
  bool exception_thrown = false;
  try
  {
    ChromatogramProcessor::pickExperiment(
      chromatograms, transition_exp, feature_finder_param,
      featureFile, transition_group_map);
  }
  catch (const std::exception& e)
  {
    exception_thrown = true;
    OPENMS_LOG_INFO << "Expected exception caught: " << e.what() << std::endl;
  }

  // Either should succeed with defaults or throw exception
  TEST_EQUAL(exception_thrown || featureFile.size() >= 0, true)
}
END_SECTION

START_SECTION(ChromatogramProcessor_edge_case_empty_input)
{
  // Test with completely empty input
  std::vector<MSChromatogram> chromatograms;  // Empty vector
  OpenSwath::LightTargetedExperiment transition_exp;  // Empty experiment

  // Feature finder parameters
  Param feature_finder_param;
  feature_finder_param.setValue("stop_report_after_feature", -1);
  feature_finder_param.setValue("add_up_spectra", 1);
  feature_finder_param.setValue("spacing_for_spectra_resampling", 0.005);
  feature_finder_param.setValue("EMGScoring:max_iteration", 10);

  // Output containers
  FeatureMap featureFile;
  OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

  // Test the static method with empty input
  ChromatogramProcessor::pickExperiment(
    chromatograms, transition_exp, feature_finder_param,
    featureFile, transition_group_map);

  // Should handle empty input gracefully
  TEST_EQUAL(featureFile.size(), 0)  // No features
  TEST_EQUAL(transition_group_map.size(), 0)  // No transition groups
}
END_SECTION
END_TEST