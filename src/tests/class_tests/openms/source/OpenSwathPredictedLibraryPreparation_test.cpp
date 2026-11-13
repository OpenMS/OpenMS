// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <set>
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
      << ">Protein;A\nPEPTIDEK\n"
      << ">ProteinB\npeptidek\n"
      << ">ProteinC\nPEPCIDEK\n"
      << ">MetProtein\nMTESTPEPK\n"
      << ">Ambiguous\nPEPXIDEK\n"
      << ">DECOY_Protein;A\nQQQQQQQK\n"
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
  TEST_TRUE(detected_decoy.success)
  TEST_EQUAL(detected_decoy.name, "DECOY_")
  TEST_TRUE(detected_decoy.is_prefix)

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
  prediction.precursor_charges = {3, 2};
  prediction.fixed_modifications = {"Carbamidomethyl (C)"};
  prediction.variable_modifications = {"Oxidation (M)"};
  prediction.max_variable_modifications = 1;
  prediction.clip_nterm_methionine = true;
  prediction.prediction_batch_size = 1;
  prediction.inference_threads = 1;
  prediction.nce = 30.0;
  prediction.instrument_index = 0;
  prediction.predict_ccs = true;

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

  TempDir scratch_parent;
  const std::string scratch_directory = scratch_parent.getPath() + "/prediction";
  const auto stats = prep.preparePredictedLibraryToPQP(
    fasta_file, output_pqp, assay, decoy, prediction, scratch_directory);
  TEST_TRUE(std::filesystem::is_directory(scratch_directory))
  TEST_TRUE(std::filesystem::is_empty(scratch_directory))
  TEST_EQUAL(stats.protein_count, 8)
  {
    // Inspect SQL values directly: the generic reader treats ';' as a separator.
    SqliteConnector connection(output_pqp);
    connection.executeStatement("CREATE TEMP TABLE exact_accessions AS SELECT * FROM PROTEIN "
                                "WHERE PROTEIN_ACCESSION IN ('Protein;A', 'DECOY_Protein;A')");
    TEST_EQUAL(connection.countTableRows("exact_accessions"), 2)
    TEST_EQUAL(connection.countTableRows("PROTEIN"), 8)
  }

  TEST_TRUE(File::exists(output_pqp))
  TEST_TRUE(stats.compound_count > 0)
  TEST_TRUE(stats.transition_count > 0)
  TEST_TRUE(stats.decoy_transition_count > 0)

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

  TEST_TRUE(saw_cam)
  TEST_TRUE(saw_clipped_methionine)
  TEST_FALSE(saw_ambiguous)
  TEST_FALSE(saw_input_decoy_sequence)

  bool saw_double_decoy_protein = false;
  for (const auto& protein : library.proteins)
  {
    saw_double_decoy_protein =
      saw_double_decoy_protein ||
      protein.id.find("DECOY_DECOY_") != std::string::npos;
  }
  TEST_FALSE(saw_double_decoy_protein)

  // Several charges and a variable modification exercise append boundaries.
  // Compare by source ID: detectingTransitionsLight already orders transitions
  // within each batch, so canonical numbering can differ between batch sizes.
  auto single_batch_prediction = prediction;
  single_batch_prediction.prediction_batch_size = 64;
  std::string single_batch_pqp;
  NEW_TMP_FILE(single_batch_pqp)
  File::remove(single_batch_pqp);
  const auto single_batch_stats = prep.preparePredictedLibraryToPQP(
    fasta_file, single_batch_pqp, assay, decoy, single_batch_prediction);

  TEST_EQUAL(single_batch_stats.protein_count, stats.protein_count)
  TEST_EQUAL(single_batch_stats.compound_count, stats.compound_count)
  TEST_EQUAL(single_batch_stats.transition_count, stats.transition_count)
  TEST_EQUAL(single_batch_stats.decoy_transition_count, stats.decoy_transition_count)
  TEST_EQUAL(single_batch_stats.identifying_transition_count, stats.identifying_transition_count)

  const auto batched_precursor_sources =
    reader.getPQPCurrentIDToTraMLIDMap(output_pqp.c_str(), "PRECURSOR");
  const auto single_precursor_sources =
    reader.getPQPCurrentIDToTraMLIDMap(single_batch_pqp.c_str(), "PRECURSOR");
  std::set<std::string> batched_precursor_names, single_precursor_names;
  for (const auto& [id, source] : batched_precursor_sources) batched_precursor_names.insert(source);
  for (const auto& [id, source] : single_precursor_sources) single_precursor_names.insert(source);
  TEST_TRUE(batched_precursor_names == single_precursor_names)

  const auto batched_transition_sources =
    reader.getPQPCurrentIDToTraMLIDMap(output_pqp.c_str(), "TRANSITION");
  const auto single_transition_sources =
    reader.getPQPCurrentIDToTraMLIDMap(single_batch_pqp.c_str(), "TRANSITION");
  std::set<std::string> batched_transition_names, single_transition_names;
  for (const auto& [id, source] : batched_transition_sources) batched_transition_names.insert(source);
  for (const auto& [id, source] : single_transition_sources) single_transition_names.insert(source);
  TEST_TRUE(batched_transition_names == single_transition_names)

  OpenSwath::LightTargetedExperiment single_library;
  reader.convertPQPToTargetedExperiment(single_batch_pqp.c_str(), single_library);
  std::map<std::string, const OpenSwath::LightCompound*> single_compounds;
  for (const auto& compound : single_library.compounds)
  {
    single_compounds.emplace(single_precursor_sources.at(compound.id), &compound);
  }
  std::vector<PeptDeepLibraryPrecursor> reference_precursors;
  for (const auto& compound : library.compounds)
  {
    const auto& source = batched_precursor_sources.at(compound.id);
    const auto& other = *single_compounds.at(source);
    TEST_EQUAL(compound.sequence, other.sequence)
    TEST_EQUAL(compound.charge, other.charge)
    TEST_REAL_SIMILAR(compound.rt, other.rt)
    TEST_REAL_SIMILAR(compound.drift_time, other.drift_time)
    TEST_TRUE(compound.protein_refs == other.protein_refs)
    if (!source.starts_with(decoy.decoy_tag))
    {
      PeptDeepLibraryPrecursor precursor;
      precursor.id = source;
      precursor.peptide = AASequence::fromString(compound.sequence);
      precursor.charge = compound.charge;
      precursor.protein_refs = compound.protein_refs;
      precursor.nce = static_cast<float>(prediction.nce);
      precursor.instrument_index = prediction.instrument_index;
      reference_precursors.push_back(std::move(precursor));
    }
  }
  std::map<std::string, const OpenSwath::LightTransition*> single_transitions;
  for (const auto& transition : single_library.transitions)
  {
    single_transitions.emplace(single_transition_sources.at(transition.transition_name), &transition);
  }
  for (const auto& transition : library.transitions)
  {
    const auto& other = *single_transitions.at(batched_transition_sources.at(transition.transition_name));
    TEST_REAL_SIMILAR(transition.library_intensity, other.library_intensity)
    TEST_REAL_SIMILAR(transition.product_mz, other.product_mz)
    TEST_REAL_SIMILAR(transition.precursor_mz, other.precursor_mz)
    TEST_REAL_SIMILAR(transition.precursor_im, other.precursor_im)
    TEST_EQUAL(transition.getAnnotation(), other.getAnnotation())
    TEST_EQUAL(transition.getDecoy(), other.getDecoy())
    TEST_EQUAL(transition.isDetectingTransition(), other.isDetectingTransition())
    TEST_EQUAL(transition.isIdentifyingTransition(), other.isIdentifyingTransition())
    TEST_EQUAL(transition.isQuantifyingTransition(), other.isQuantifyingTransition())
  }

  // Check persisted RT, mobility and intensities against direct inference, not
  // just another run through the same spill reader.
  PeptDeepLibraryPredictor::Config reference_config;
  reference_config.batch_size = prediction.prediction_batch_size;
  reference_config.intra_op_threads = 1;
  reference_config.predict_ccs = true;
  const auto reference = PeptDeepLibraryPredictor(reference_config).predict(reference_precursors);
  for (const auto& compound : reference.compounds)
  {
    const auto& persisted = *single_compounds.at(compound.id);
    TEST_REAL_SIMILAR(compound.rt, persisted.rt)
    TEST_REAL_SIMILAR(compound.drift_time, persisted.drift_time)
  }
  Size checked_transitions = 0;
  for (const auto& transition : reference.transitions)
  {
    const auto it = single_transitions.find(transition.transition_name);
    if (it == single_transitions.end()) continue; // assay filtering keeps only the top fragments
    TEST_REAL_SIMILAR(transition.library_intensity, it->second->library_intensity)
    ++checked_transitions;
  }
  TEST_EQUAL(checked_transitions, stats.transition_count - stats.decoy_transition_count)

  // Force a failure after prediction/reload and verify exception-safe scratch cleanup.
  std::string missing_output_parent;
  NEW_TMP_FILE(missing_output_parent)
  File::remove(missing_output_parent);
  TEST_EXCEPTION(Exception::SqlOperationFailed,
    prep.preparePredictedLibraryToPQP(fasta_file, missing_output_parent + "/output.pqp",
                                     assay, decoy, prediction, scratch_directory))
  TEST_TRUE(std::filesystem::is_empty(scratch_directory))

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
  TEST_FALSE(DecoyHelper::findDecoyString(sparse_decoy_scan, true).success)

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
  TEST_FALSE(saw_sparse_decoy_sequence)
  bool saw_sparse_double_decoy_protein = false;
  for (const auto& protein : sparse_decoy_library.proteins)
  {
    saw_sparse_double_decoy_protein =
      saw_sparse_double_decoy_protein ||
      protein.id.find("DECOY_DECOY_") != std::string::npos;
  }
  TEST_FALSE(saw_sparse_double_decoy_protein)

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

  // Library callers do not get TOPP option checks, so invalid settings must be
  // rejected before any model or FASTA work starts.
  auto short_peptides = prediction;
  short_peptides.min_peptide_length = 1;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.preparePredictedLibraryToPQP(fasta_file, invalid_output_pqp, assay, decoy, short_peptides),
    "PredictedLibraryParameters::min_peptide_length must be at least 2.")
  auto no_max_length = prediction;
  no_max_length.max_peptide_length = 0;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.preparePredictedLibraryToPQP(fasta_file, invalid_output_pqp, assay, decoy, no_max_length),
    "PredictedLibraryParameters::max_peptide_length must be >= min_peptide_length.")
  auto bad_instrument = prediction;
  bad_instrument.instrument_index = 8;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.preparePredictedLibraryToPQP(fasta_file, invalid_output_pqp, assay, decoy, bad_instrument),
    "PredictedLibraryParameters::instrument_index must be between 0 and 7.")
  for (const Size batch_size : {Size{1}, Size{64}})
  {
    auto duplicate_charges = prediction;
    duplicate_charges.precursor_charges = {2, 2};
    duplicate_charges.prediction_batch_size = batch_size;
    TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
      prep.preparePredictedLibraryToPQP(fasta_file, invalid_output_pqp, assay, decoy, duplicate_charges),
      "PredictedLibraryParameters::precursor_charges must contain unique charges.")
  }
  auto bad_enzyme = prediction;
  bad_enzyme.enzyme = "NoSuchProtease";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter,
    prep.preparePredictedLibraryToPQP(fasta_file, invalid_output_pqp, assay, decoy, bad_enzyme),
    "PredictedLibraryParameters::enzyme 'NoSuchProtease' is not a known protease.")
  TEST_FALSE(File::exists(invalid_output_pqp))
}
END_SECTION

START_SECTION([EXTRA] generateDecoys with predict_decoys gives each decoy its own predicted assay)
{
  // A target library predicted at a known setting (NCE 27, Lumos, no CCS) without decoys.
  std::string fasta_file;
  NEW_TMP_FILE(fasta_file)
  {
    std::ofstream fasta(fasta_file);
    fasta << ">P1\nMSAEEVKLAGTPDDLVRAHSELLGQYK\n"
          << ">P2\nGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEK\n"
          << ">P3\nVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK\n"
          << ">P4\nTEFKEAFSLFDKDGDGTITTKELGTVMR\n"
          << ">P5\nSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMAR\n"
          << ">P6\nADQLTEEQIAEFKEAFSLFDKDGNGTISAAELR\n";
  }
  OpenSwathLibraryPreparation prep;
  prep.setLogType(ProgressLogger::NONE);

  OpenSwathLibraryPreparation::PredictedLibraryParameters prediction;
  prediction.enzyme = "Trypsin";
  prediction.missed_cleavages = 1;
  prediction.min_peptide_length = 7;
  prediction.max_peptide_length = 30;
  prediction.precursor_charges = {2, 3};
  prediction.fixed_modifications = {"Carbamidomethyl (C)"};
  prediction.max_variable_modifications = 0;
  prediction.clip_nterm_methionine = false;
  prediction.prediction_batch_size = 64;
  prediction.inference_threads = 1;
  prediction.nce = 27.0;
  prediction.instrument_index = ML::instrumentIndex("Lumos");
  prediction.predict_ccs = false;

  OpenSwathLibraryPreparation::AssayGeneratorParameters assay;
  assay.min_transitions = 6;
  assay.max_transitions = 6;
  assay.precursor_lower_mz_limit = 300.0;
  assay.precursor_upper_mz_limit = 1500.0;
  assay.product_lower_mz_limit = 150.0;
  assay.product_upper_mz_limit = 2000.0;

  OpenSwathLibraryPreparation::DecoyGeneratorParameters separate_decoys;
  separate_decoys.method = "pseudo-reverse";
  separate_decoys.min_decoy_fraction = 0.0;
  separate_decoys.separate = true;

  // Targets only: build with decoys, then keep the targets (the predicted path always adds decoys).
  std::string with_decoys;
  NEW_TMP_FILE_EXT(with_decoys, ".pqp")
  File::remove(with_decoys);
  OpenSwathLibraryPreparation::DecoyGeneratorParameters copy_decoys = separate_decoys;
  copy_decoys.separate = false;
  prep.preparePredictedLibraryToPQP(fasta_file, with_decoys, assay, copy_decoys, prediction);
  OpenSwath::LightTargetedExperiment library;
  TransitionPQPFile pqp;
  pqp.convertPQPToTargetedExperiment(with_decoys.c_str(), library);
  OpenSwath::LightTargetedExperiment targets;
  for (const auto& compound : library.compounds)
  {
    if (compound.id.rfind("DECOY_", 0) != 0) targets.compounds.push_back(compound);
  }
  for (const auto& transition : library.transitions)
  {
    if (!transition.getDecoy()) targets.transitions.push_back(transition);
  }
  for (const auto& protein : library.proteins)
  {
    if (protein.id.rfind("DECOY_", 0) != 0) targets.proteins.push_back(protein);
  }
  TEST_TRUE(targets.compounds.size() >= 20)
  std::string targets_pqp;
  NEW_TMP_FILE_EXT(targets_pqp, ".pqp")
  File::remove(targets_pqp);
  pqp.convertLightTargetedExperimentToPQP(targets_pqp.c_str(), targets);

  // Predicted decoys with explicit settings and with the settings estimated from the targets.
  OpenSwathLibraryPreparation::DecoyGeneratorParameters predicted = separate_decoys;
  predicted.predict_decoys = true;
  predicted.predict_threads = 1;
  std::map<std::string, std::string> outputs;
  for (const bool estimate : {false, true})
  {
    OpenSwathLibraryPreparation::DecoyGeneratorParameters parameters = predicted;
    if (!estimate)
    {
      parameters.predict_nce = 27.0;
      parameters.predict_instrument = "Lumos";
    }
    std::string out;
    NEW_TMP_FILE_EXT(out, ".pqp")
    File::remove(out);
    const auto stats = prep.generateDecoys(targets_pqp, FileTypes::PQP, out, FileTypes::PQP, parameters);
    TEST_TRUE(stats.decoy_transition_count > 0)
    outputs[estimate ? "estimated" : "explicit"] = out;
  }

  OpenSwath::LightTargetedExperiment decoys;
  pqp.convertPQPToTargetedExperiment(outputs["explicit"].c_str(), decoys);
  OpenSwath::LightTargetedExperiment estimated;
  pqp.convertPQPToTargetedExperiment(outputs["estimated"].c_str(), estimated);

  // Reference: the decoy sequences predicted directly at the library's setting.
  PeptDeepLibraryPredictor::Config config;
  config.intra_op_threads = 1;
  config.predict_ccs = false;
  PeptDeepLibraryPredictor predictor(config);
  std::vector<PeptDeepLibraryPrecursor> inputs;
  for (const auto& compound : decoys.compounds)
  {
    PeptDeepLibraryPrecursor precursor;
    precursor.peptide = AASequence::fromString(compound.sequence);
    precursor.id = compound.sequence + "/" + std::to_string(compound.charge);
    precursor.charge = compound.charge;
    precursor.nce = 27.0f;
    precursor.instrument_index = ML::instrumentIndex("Lumos");
    inputs.push_back(precursor);
  }
  const auto reference = predictor.predict(inputs);
  std::map<std::string, std::map<std::string, double>> reference_intensity;
  std::map<std::string, double> reference_rt;
  for (const auto& transition : reference.transitions)
  {
    reference_intensity[transition.peptide_ref][transition.getAnnotation()] = transition.library_intensity;
  }
  for (const auto& compound : reference.compounds) reference_rt[compound.id] = compound.rt;

  std::map<std::string, const OpenSwath::LightCompound*> decoy_by_id;
  for (const auto& compound : decoys.compounds) decoy_by_id[compound.id] = &compound;
  std::map<std::string, std::vector<const OpenSwath::LightTransition*>> decoy_transitions;
  for (const auto& transition : decoys.transitions) decoy_transitions[transition.peptide_ref].push_back(&transition);

  Size checked = 0;
  bool all_six = true, all_decoy = true, top_is_one = true, intensities_match = true, rt_match = true, top_k = true;
  for (const auto& [id, transitions] : decoy_transitions)
  {
    const auto* compound = decoy_by_id.at(id);
    const std::string key = compound->sequence + "/" + std::to_string(compound->charge);
    const auto& predicted_intensity = reference_intensity.at(key);
    all_six = all_six && transitions.size() == 6;
    double top = 0.0, max_kept = 0.0, min_kept = 1e300;
    for (const auto* transition : transitions)
    {
      max_kept = std::max(max_kept, predicted_intensity.at(transition->getAnnotation()));
    }
    for (const auto* transition : transitions)
    {
      all_decoy = all_decoy && transition->getDecoy();
      top = std::max(top, transition->library_intensity);
      // Kept intensities are the prediction normalised to the most intense fragment kept.
      const double expected = predicted_intensity.at(transition->getAnnotation()) / max_kept;
      intensities_match = intensities_match && std::abs(transition->library_intensity - expected) < 1e-4;
      min_kept = std::min(min_kept, predicted_intensity.at(transition->getAnnotation()));
    }
    // Every fragment not kept but inside the targets' window is no more intense than the weakest kept one.
    Size stronger_dropped = 0;
    for (const auto& [annotation, value] : predicted_intensity)
    {
      bool kept = false;
      for (const auto* transition : transitions) kept = kept || transition->getAnnotation() == annotation;
      if (!kept && value > min_kept + 1e-6) ++stronger_dropped;
    }
    // Stronger fragments can only be dropped for lying outside the targets' fragment m/z window.
    top_k = top_k && stronger_dropped <= 2;
    top_is_one = top_is_one && std::abs(top - 1.0) < 1e-6;
    rt_match = rt_match && std::abs(compound->rt - reference_rt.at(key)) < 1e-3;
    ++checked;
  }
  TEST_TRUE(checked >= 20)
  TEST_TRUE(all_six)
  TEST_TRUE(all_decoy)
  TEST_TRUE(top_is_one)
  TEST_TRUE(intensities_match)
  TEST_TRUE(rt_match)
  TEST_TRUE(top_k)

  // Estimating the setting from the targets recovers the one they were predicted at: the same decoy assays.
  TEST_EQUAL(estimated.transitions.size(), decoys.transitions.size())
  std::map<std::string, double> explicit_intensity;
  for (const auto& transition : decoys.transitions) explicit_intensity[transition.transition_name] = transition.library_intensity;
  bool same_assays = true;
  for (const auto& transition : estimated.transitions)
  {
    const auto it = explicit_intensity.find(transition.transition_name);
    same_assays = same_assays && it != explicit_intensity.end() && std::abs(it->second - transition.library_intensity) < 1e-6;
  }
  TEST_TRUE(same_assays)

  // TraML input has no predicted-decoy path.
  std::string traml_out;
  NEW_TMP_FILE_EXT(traml_out, ".TraML")
  TEST_EXCEPTION(Exception::InvalidParameter,
    prep.generateDecoys(targets_pqp, FileTypes::PQP, traml_out, FileTypes::TRAML, predicted))
}
END_SECTION

END_TEST
