// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathPredictedLibraryPreparation, "$Id$")

START_SECTION((OpenSwathLibraryPreparation::LibraryStats preparePredictedLibraryToPQP(...)))
{
  std::string fasta_file;
  NEW_TMP_FILE(fasta_file)
  {
    std::ofstream fasta(fasta_file);
    fasta
      << ">ProteinA\nPEPTIDEK\n"
      << ">ProteinB\npeptidek\n"
      << ">ProteinC\nPEPCIDEK\n"
      << ">MetProtein\nMTESTPEPK\n"
      << ">Ambiguous\nPEPXIDEK\n"
      << ">DECOY_ProteinA\nQQQQQQQK\n"
      << ">DECOY_ProteinB\nQQQQQQQK\n"
      << ">DECOY_ProteinC\nQQQCQQQK\n"
      << ">DECOY_MetProtein\nQQQWQQQK\n";
  }

  // DecoyHelper intentionally requires the same decoy affix on at least 40%
  // of proteins. Keep the fixture above that threshold and assert the helper
  // contract explicitly before testing predicted-library preparation.
  FASTAContainer<TFI_File> decoy_scan(fasta_file);
  const DecoyHelper::Result detected_decoy =
    DecoyHelper::findDecoyString(decoy_scan, true);
  TEST_EQUAL(detected_decoy.success, true)
  TEST_EQUAL(detected_decoy.name, "DECOY_")
  TEST_EQUAL(detected_decoy.is_prefix, true)

  std::string output_pqp;
  NEW_TMP_FILE(output_pqp)
  File::remove(output_pqp);

  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  OpenSwathLibraryPreparation::PredictedLibraryParameters prediction;
  prediction.enzyme = "Trypsin";
  prediction.missed_cleavages = 2;
  prediction.min_peptide_length = 7;
  prediction.max_peptide_length = 30;
  prediction.precursor_charges = {2};
  prediction.fixed_modifications = {"Carbamidomethyl (C)"};
  prediction.variable_modifications = {};
  prediction.max_variable_modifications = 0;
  prediction.clip_nterm_methionine = true;
  prediction.prediction_batch_size = 1;
  prediction.inference_threads = 1;
  prediction.nce = 30.0;
  prediction.instrument_index = 0;
  prediction.predict_ccs = false;

  OpenSwathLibraryPreparation::AssayGeneratorParameters assay;
  assay.min_transitions = 1;
  assay.max_transitions = 6;
  assay.precursor_lower_mz_limit = 0.0;
  assay.precursor_upper_mz_limit = 2000.0;
  assay.product_lower_mz_limit = 0.0;
  assay.product_upper_mz_limit = 2000.0;

  OpenSwathLibraryPreparation::DecoyGeneratorParameters decoy;
  decoy.method = "reverse";
  decoy.min_decoy_fraction = 0.0;

  const auto stats = prep.preparePredictedLibraryToPQP(
    fasta_file, output_pqp, assay, decoy, prediction);

  TEST_EQUAL(File::exists(output_pqp), true)
  TEST_EQUAL(stats.compound_count > 0, true)
  TEST_EQUAL(stats.transition_count > 0, true)
  TEST_EQUAL(stats.decoy_transition_count > 0, true)

  OpenSwath::LightTargetedExperiment library;
  TransitionPQPFile reader;
  reader.convertPQPToTargetedExperiment(output_pqp.c_str(), library);

  bool saw_cam = false;
  bool saw_clipped_methionine = false;
  bool saw_ambiguous = false;
  bool saw_input_decoy_sequence = false;
  for (const auto& compound : library.compounds)
  {
    saw_cam = saw_cam || compound.sequence.find("UniMod:4") != std::string::npos;
    saw_clipped_methionine =
      saw_clipped_methionine || compound.sequence.find("TESTPEPK") != std::string::npos;
    saw_ambiguous = saw_ambiguous || compound.sequence.find('X') != std::string::npos;
    saw_input_decoy_sequence =
      saw_input_decoy_sequence ||
      compound.sequence.find("QQQQQQQK") != std::string::npos ||
      compound.sequence.find("QQQCQQQK") != std::string::npos ||
      compound.sequence.find("QQQWQQQK") != std::string::npos;
  }

  TEST_EQUAL(saw_cam, true)
  TEST_EQUAL(saw_clipped_methionine, true)
  TEST_EQUAL(saw_ambiguous, false)
  TEST_EQUAL(saw_input_decoy_sequence, false)

  bool saw_double_decoy_protein = false;
  for (const auto& protein : library.proteins)
  {
    saw_double_decoy_protein =
      saw_double_decoy_protein ||
      protein.id.find("DECOY_DECOY_") != std::string::npos;
  }
  TEST_EQUAL(saw_double_decoy_protein, false)

  // Entries carrying the configured decoy tag are skipped even when they are
  // too rare for DecoyHelper to detect a database-wide decoy affix.
  std::string sparse_decoy_fasta;
  NEW_TMP_FILE(sparse_decoy_fasta)
  {
    std::ofstream fasta(sparse_decoy_fasta);
    fasta
      << ">ProteinA\nPEPTIDEK\n"
      << ">ProteinB\nPEPCIDEK\n"
      << ">ProteinC\nSDVEEQK\n"
      << ">ProteinD\nYPELLANR\n"
      << ">DECOY_ProteinA\nQQQQQQQK\n";
  }
  FASTAContainer<TFI_File> sparse_decoy_scan(sparse_decoy_fasta);
  TEST_EQUAL(DecoyHelper::findDecoyString(sparse_decoy_scan, true).success, false)

  std::string sparse_decoy_pqp;
  NEW_TMP_FILE(sparse_decoy_pqp)
  File::remove(sparse_decoy_pqp);
  prep.preparePredictedLibraryToPQP(
    sparse_decoy_fasta, sparse_decoy_pqp, assay, decoy, prediction);

  OpenSwath::LightTargetedExperiment sparse_decoy_library;
  reader.convertPQPToTargetedExperiment(sparse_decoy_pqp.c_str(), sparse_decoy_library);
  bool saw_sparse_decoy_sequence = false;
  for (const auto& compound : sparse_decoy_library.compounds)
  {
    saw_sparse_decoy_sequence =
      saw_sparse_decoy_sequence ||
      compound.sequence.find("QQQQQQQK") != std::string::npos;
  }
  TEST_EQUAL(saw_sparse_decoy_sequence, false)
  bool saw_sparse_double_decoy_protein = false;
  for (const auto& protein : sparse_decoy_library.proteins)
  {
    saw_sparse_double_decoy_protein =
      saw_sparse_double_decoy_protein ||
      protein.id.find("DECOY_DECOY_") != std::string::npos;
  }
  TEST_EQUAL(saw_sparse_double_decoy_protein, false)

  // Fallback UIS/SWATH construction divides the precursor m/z range by the
  // precursor threshold. Reject a non-positive threshold before that division.
  // IPF requires a Unimod file; ModificationsDB is already initialized above,
  // so reuse it to get past that check and reach the threshold guard.
  OpenSwathLibraryPreparation::AssayGeneratorParameters invalid_uis_assay = assay;
  invalid_uis_assay.enable_ipf = true;
  invalid_uis_assay.unimod_file = File::find("CHEMISTRY/unimod.xml");
  invalid_uis_assay.reuse_existing_modifications_db = true;
  invalid_uis_assay.enable_swath_specifity = false;
  invalid_uis_assay.swathes.clear();
  invalid_uis_assay.precursor_mz_threshold = 0.0;

  std::string invalid_output_pqp;
  NEW_TMP_FILE(invalid_output_pqp)
  File::remove(invalid_output_pqp);

  TEST_EXCEPTION_WITH_MESSAGE(
    Exception::InvalidParameter,
    prep.preparePredictedLibraryToPQP(
      fasta_file, invalid_output_pqp, invalid_uis_assay, decoy, prediction),
    "AssayGeneratorParameters::precursor_mz_threshold must be greater than zero "
    "when constructing fallback UIS SWATH windows.")
}
END_SECTION

END_TEST
