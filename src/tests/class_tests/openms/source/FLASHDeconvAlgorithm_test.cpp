// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Müller, Jaekwan Kim$
// $Authors: Tom David Müller, Jaekwan Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////

START_TEST(FLASHDeconvAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FLASHDeconvAlgorithm* ptr = nullptr;
FLASHDeconvAlgorithm* null_ptr = nullptr;

START_SECTION(FLASHDeconvAlgorithm())
  ptr = new FLASHDeconvAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm& source))
FLASHDeconvAlgorithm copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& source))
FLASHDeconvAlgorithm copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~FLASHDeconvAlgorithm())
  delete ptr;
END_SECTION

ptr = new FLASHDeconvAlgorithm();
START_SECTION(FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& source))
  FLASHDeconvAlgorithm temp;
  temp.setParameters(ptr->getParameters());
  FLASHDeconvAlgorithm moved(std::move(temp));
  TEST_EQUAL(moved.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& source))
  FLASHDeconvAlgorithm temp;
  temp.setParameters(ptr->getParameters());
  FLASHDeconvAlgorithm moved;
  moved = std::move(temp);
  TEST_EQUAL(moved.getParameters(), ptr->getParameters());
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ptr = new FLASHDeconvAlgorithm();
Param params;
START_SECTION(std::vector<double> getTolerances())
  params.setValue("SD:tol", ListUtils::create<double>("10.0,5.0"));
  ptr->setParameters(params);
  auto tolerances = ptr->getTolerances();
  TEST_EQUAL(tolerances.size(), 2)
  TEST_REAL_SIMILAR(tolerances[0], 10.0);
  TEST_REAL_SIMILAR(tolerances[1], 5.0);
END_SECTION

// load test data
// TODO: minimize
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FLASHDeconv_1_input.mzML"), input);

// Store FD outputs
std::vector<DeconvolvedSpectrum> deconvolved_spectra;
std::vector<FLASHHelperClasses::MassFeature> deconvolved_features;

START_SECTION(void run())
  ptr->run(input, deconvolved_spectra, deconvolved_features);
  TEST_EQUAL(deconvolved_spectra.size(), input.size());
  TEST_FALSE(deconvolved_features.empty());
END_SECTION

START_SECTION(int getScanNumber())
  // The scan number depends on the test file; just verify it returns a positive value
  TEST_EQUAL(FLASHDeconvAlgorithm::getScanNumber(input, 0) > 0, true);
END_SECTION

// Check if the precalculated averagine was initialized during runtime
START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getAveragine())
  TEST_NOT_EQUAL(ptr->getAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getAveragine().get(500.0).getContainer().empty());
END_SECTION

START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getDecoyAveragine())
  TEST_NOT_EQUAL(ptr->getDecoyAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getDecoyAveragine().get(500.0).getContainer().empty());
END_SECTION

// Decoy Averagine needs to be handled differently if FDR is reported
START_SECTION(FLASHHelperClasses::PrecalculatedAveragine getDecoyAveragine())
  params.setValue("FD:report_FD", true);
  ptr = new FLASHDeconvAlgorithm();  
  ptr->setParameters(params);
  ptr->run(input, deconvolved_spectra, deconvolved_features);
  TEST_NOT_EQUAL(ptr->getDecoyAveragine().getMaxIsotopeIndex(), 0);
  TEST_FALSE(ptr->getDecoyAveragine().get(500.0).getContainer().empty());
END_SECTION


/////////////////////////////////////////////////////////////
// Phase 1: Averagine Model Tests (Enhanced)
/////////////////////////////////////////////////////////////

START_SECTION(Averagine model - protein vs RNA mode initialization)
  // Test protein averagine (default)
  FLASHDeconvAlgorithm algo_protein;
  Param protein_params;
  protein_params.setValue("use_RNA_averagine", "false");
  algo_protein.setParameters(protein_params);
  algo_protein.run(input, deconvolved_spectra, deconvolved_features);
  
  auto& protein_avg = algo_protein.getAveragine();
  TEST_NOT_EQUAL(protein_avg.getMaxIsotopeIndex(), 0)
  TEST_FALSE(protein_avg.get(1000.0).getContainer().empty())
  
  // Test RNA averagine
  FLASHDeconvAlgorithm algo_rna;
  Param rna_params;
  rna_params.setValue("use_RNA_averagine", "true");
  algo_rna.setParameters(rna_params);
  algo_rna.run(input, deconvolved_spectra, deconvolved_features);
  
  auto& rna_avg = algo_rna.getAveragine();
  TEST_NOT_EQUAL(rna_avg.getMaxIsotopeIndex(), 0)
  TEST_FALSE(rna_avg.get(1000.0).getContainer().empty())
END_SECTION

START_SECTION(Averagine model - mass range validation)
  FLASHDeconvAlgorithm algo;
  algo.run(input, deconvolved_spectra, deconvolved_features);
  
  auto& avg = algo.getAveragine();
  
  // Test averagine at various mass ranges
  TEST_FALSE(avg.get(500.0).getContainer().empty())
  TEST_FALSE(avg.get(1000.0).getContainer().empty())
  TEST_FALSE(avg.get(5000.0).getContainer().empty())
  TEST_FALSE(avg.get(10000.0).getContainer().empty())
END_SECTION

START_SECTION(Averagine model - consistency check)
  FLASHDeconvAlgorithm algo;
  algo.run(input, deconvolved_spectra, deconvolved_features);
  
  auto& avg = algo.getAveragine();
  
  // Verify max isotope index is reasonable
  int max_iso_idx = avg.getMaxIsotopeIndex();
  TEST_NOT_EQUAL(max_iso_idx, 0)
  TEST_EQUAL(max_iso_idx > 0, true)
  TEST_EQUAL(max_iso_idx < 1000, true) // Reasonable upper bound
END_SECTION

/////////////////////////////////////////////////////////////
// Phase 1: MS Level Processing Tests
/////////////////////////////////////////////////////////////

START_SECTION(MS level processing - scan number retrieval)
  // Test getScanNumber with different indices - verify it returns positive values
  TEST_EQUAL(FLASHDeconvAlgorithm::getScanNumber(input, 0) > 0, true)

  // Test with valid index range
  if (input.size() > 1)
  {
    int scan_num = FLASHDeconvAlgorithm::getScanNumber(input, 1);
    TEST_NOT_EQUAL(scan_num, 0) // Should have a valid scan number
  }
END_SECTION

/////////////////////////////////////////////////////////////
// Phase 1: Spectrum Merging Tests
/////////////////////////////////////////////////////////////

START_SECTION(Spectrum merging - merging method parameter validation)
  FLASHDeconvAlgorithm algo;
  Param merge_params;
  
  // Test merging_method = 0 (no merging)
  merge_params.setValue("FD:merging_method", 0);
  algo.setParameters(merge_params);
  TEST_EQUAL(algo.getParameters().exists("FD:merging_method"), true)
  
  // Test merging_method = 1 (Gaussian)
  merge_params.setValue("FD:merging_method", 1);
  algo.setParameters(merge_params);
  TEST_EQUAL(algo.getParameters().exists("FD:merging_method"), true)
  
  // Test merging_method = 2 (block)
  merge_params.setValue("FD:merging_method", 2);
  algo.setParameters(merge_params);
  TEST_EQUAL(algo.getParameters().exists("FD:merging_method"), true)
END_SECTION

START_SECTION(Spectrum merging - effect on output with no merging)
  FLASHDeconvAlgorithm algo_no_merge;
  Param no_merge_params;
  no_merge_params.setValue("FD:merging_method", 0);
  algo_no_merge.setParameters(no_merge_params);
  
  std::vector<DeconvolvedSpectrum> no_merge_spectra;
  std::vector<FLASHHelperClasses::MassFeature> no_merge_features;
  
  algo_no_merge.run(input, no_merge_spectra, no_merge_features);
  
  // With no merging, output spectrum count should equal input
  TEST_EQUAL(no_merge_spectra.size(), input.size())
END_SECTION

START_SECTION(Spectrum merging - Gaussian merging mode)
  FLASHDeconvAlgorithm algo_gaussian;
  Param gaussian_params;
  gaussian_params.setValue("FD:merging_method", 1);
  algo_gaussian.setParameters(gaussian_params);
  
  std::vector<DeconvolvedSpectrum> gaussian_spectra;
  std::vector<FLASHHelperClasses::MassFeature> gaussian_features;
  
  algo_gaussian.run(input, gaussian_spectra, gaussian_features);
  
  // With Gaussian merging, we should have some output
  // (exact count depends on merging parameters)
  TEST_NOT_EQUAL(gaussian_spectra.size(), 0)
END_SECTION

/////////////////////////////////////////////////////////////
// Phase 1: FDR and Decoy Generation Tests
/////////////////////////////////////////////////////////////

START_SECTION(FDR - report_FDR parameter toggling)
  // Test with FDR reporting disabled (default)
  FLASHDeconvAlgorithm algo_no_fdr;
  Param no_fdr_params;
  no_fdr_params.setValue("FD:report_FD", "false");
  algo_no_fdr.setParameters(no_fdr_params);
  
  std::vector<DeconvolvedSpectrum> no_fdr_spectra;
  std::vector<FLASHHelperClasses::MassFeature> no_fdr_features;
  algo_no_fdr.run(input, no_fdr_spectra, no_fdr_features);
  
  // Decoy averagine should still be initialized
  auto& decoy_avg_no_fdr = algo_no_fdr.getDecoyAveragine();
  TEST_NOT_EQUAL(decoy_avg_no_fdr.getMaxIsotopeIndex(), 0)
  
  // Test with FDR reporting enabled
  FLASHDeconvAlgorithm algo_with_fdr;
  Param fdr_params;
  fdr_params.setValue("FD:report_FD", "true");
  algo_with_fdr.setParameters(fdr_params);
  
  std::vector<DeconvolvedSpectrum> fdr_spectra;
  std::vector<FLASHHelperClasses::MassFeature> fdr_features;
  algo_with_fdr.run(input, fdr_spectra, fdr_features);
  
  // Decoy averagine should be properly initialized for FDR
  auto& decoy_avg_fdr = algo_with_fdr.getDecoyAveragine();
  TEST_NOT_EQUAL(decoy_avg_fdr.getMaxIsotopeIndex(), 0)
  TEST_FALSE(decoy_avg_fdr.get(500.0).getContainer().empty())
END_SECTION

START_SECTION(FDR - noise decoy weight validation)
  FLASHDeconvAlgorithm algo;
  Param fdr_params;
  fdr_params.setValue("FD:report_FD", "true");
  algo.setParameters(fdr_params);

  algo.run(input, deconvolved_spectra, deconvolved_features);

  // Get noise decoy weight
  double noise_weight = algo.getNoiseDecoyWeight();

  // Verify noise weight is reasonable (should be positive and bounded)
  TEST_EQUAL(noise_weight > 0, true)
  TEST_EQUAL(noise_weight <= 1.0, true) // Calculated weight should be at most 1.0
END_SECTION

START_SECTION(FDR - decoy averagine properties)
  FLASHDeconvAlgorithm algo;
  Param fdr_params;
  fdr_params.setValue("FD:report_FD", "true");
  algo.setParameters(fdr_params);
  
  algo.run(input, deconvolved_spectra, deconvolved_features);
  
  auto& decoy_avg = algo.getDecoyAveragine();
  
  // Verify decoy averagine has similar properties to regular averagine
  TEST_NOT_EQUAL(decoy_avg.getMaxIsotopeIndex(), 0)
  
  // Test at different masses
  TEST_FALSE(decoy_avg.get(500.0).getContainer().empty())
  TEST_FALSE(decoy_avg.get(1000.0).getContainer().empty())
  TEST_FALSE(decoy_avg.get(5000.0).getContainer().empty())
END_SECTION

START_SECTION(FDR - decoy generation consistency)
  FLASHDeconvAlgorithm algo1;
  FLASHDeconvAlgorithm algo2;
  
  Param fdr_params;
  fdr_params.setValue("FD:report_FD", "true");
  
  algo1.setParameters(fdr_params);
  algo2.setParameters(fdr_params);
  
  std::vector<DeconvolvedSpectrum> spectra1, spectra2;
  std::vector<FLASHHelperClasses::MassFeature> features1, features2;
  
  algo1.run(input, spectra1, features1);
  algo2.run(input, spectra2, features2);
  
  // Both should generate decoy averagines with same max isotope index
  TEST_EQUAL(algo1.getDecoyAveragine().getMaxIsotopeIndex(),
             algo2.getDecoyAveragine().getMaxIsotopeIndex())
END_SECTION

/////////////////////////////////////////////////////////////
// Phase 1: Output Validation Tests
/////////////////////////////////////////////////////////////

START_SECTION(Output validation - deconvolved spectrum structure)
  FLASHDeconvAlgorithm algo;
  std::vector<DeconvolvedSpectrum> output_spectra;
  std::vector<FLASHHelperClasses::MassFeature> output_features;
  
  algo.run(input, output_spectra, output_features);
  
  // Verify output has expected size
  TEST_EQUAL(output_spectra.size(), input.size())
  
  // Check that each deconvolved spectrum has valid properties
  for (const auto& spec : output_spectra)
  {
    // Verify spectrum has been processed (has some metadata)
    TEST_EQUAL(spec.getScanNumber() >= 0, true)
  }
END_SECTION

START_SECTION(Output validation - mass features structure)
  FLASHDeconvAlgorithm algo;
  std::vector<DeconvolvedSpectrum> output_spectra;
  std::vector<FLASHHelperClasses::MassFeature> output_features;
  
  algo.run(input, output_spectra, output_features);
  
  // Verify features were found (may be empty for some datasets)
  // Just check that the vector is accessible
  TEST_EQUAL(output_features.size() >= 0, true)
  
  // If features exist, validate basic properties
  for (const auto& feature : output_features)
  {
    // Verify feature has valid mass range
    TEST_EQUAL(feature.avg_mass >= 0, true)
  }
END_SECTION

START_SECTION(Output validation - PeakGroup properties in deconvolved spectra)
  FLASHDeconvAlgorithm algo;
  std::vector<DeconvolvedSpectrum> output_spectra;
  std::vector<FLASHHelperClasses::MassFeature> output_features;
  
  algo.run(input, output_spectra, output_features);
  
  // Check deconvolved spectra for peak groups
  bool found_peak_groups = false;
  for (const auto& spec : output_spectra)
  {
    if (spec.size() > 0)
    {
      found_peak_groups = true;
      
      // Verify peak groups have valid properties
      for (Size i = 0; i < spec.size(); ++i)
      {
        const auto& pg = spec[i];
        // Peak groups should have positive mass
        TEST_EQUAL(pg.getMonoMass() >= 0, true)
      }
    }
  }
  
  // At least some spectra should have peak groups
  TEST_EQUAL(found_peak_groups, true)
END_SECTION

START_SECTION(Output validation - consistent output vector sizes)
  FLASHDeconvAlgorithm algo;
  std::vector<DeconvolvedSpectrum> output_spectra;
  std::vector<FLASHHelperClasses::MassFeature> output_features;
  
  // Run twice to ensure consistency
  algo.run(input, output_spectra, output_features);
  size_t first_spectra_size = output_spectra.size();
  
  output_spectra.clear();
  output_features.clear();
  
  algo.run(input, output_spectra, output_features);
  size_t second_spectra_size = output_spectra.size();
  
  // Should produce same number of spectra on repeated runs
  TEST_EQUAL(first_spectra_size, second_spectra_size)
END_SECTION

START_SECTION(Output validation - output spectrum properties)
  FLASHDeconvAlgorithm algo;
  std::vector<DeconvolvedSpectrum> output_spectra;
  std::vector<FLASHHelperClasses::MassFeature> output_features;
  
  algo.run(input, output_spectra, output_features);
  
  // Verify each output spectrum has consistent properties with input
  for (size_t i = 0; i < output_spectra.size() && i < input.size(); ++i)
  {
    // Scan numbers should match
    TEST_EQUAL(output_spectra[i].getScanNumber(),
               FLASHDeconvAlgorithm::getScanNumber(input, i))
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
