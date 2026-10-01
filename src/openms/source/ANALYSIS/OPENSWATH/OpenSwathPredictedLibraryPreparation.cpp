// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Exception.h>

#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>
#endif

namespace OpenMS
{
#ifdef WITH_ONNX
  namespace
  {
    struct PredictedCandidate_
    {
      AASequence peptide;
      Int charge = 0;
      std::set<std::string> protein_refs;
    };

    void prepareDetectionTransitions_(
      OpenSwath::LightTargetedExperiment& experiment,
      const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters,
      const ProgressLogger::LogType log_type)
    {
      MRMAssay assays;
      assays.setLogType(log_type);
      assays.reannotateTransitionsLight(
        experiment,
        parameters.precursor_mz_threshold,
        parameters.product_mz_threshold,
        parameters.allowed_fragment_types,
        parameters.allowed_fragment_charges,
        parameters.enable_detection_specific_losses,
        parameters.enable_detection_unspecific_losses);
      assays.restrictTransitionsLight(
        experiment,
        parameters.product_lower_mz_limit,
        parameters.product_upper_mz_limit,
        parameters.swathes);
      assays.detectingTransitionsLight(
        experiment,
        parameters.min_transitions,
        parameters.max_transitions);
    }

    std::vector<std::pair<double, double>> buildPredictedUISSwathes_(
      const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters)
    {
      if (parameters.enable_swath_specifity && !parameters.swathes.empty())
      {
        return parameters.swathes;
      }

      if (parameters.precursor_mz_threshold <= 0.0)
      {
        throw Exception::InvalidParameter(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "AssayGeneratorParameters::precursor_mz_threshold must be greater than zero "
          "when constructing fallback UIS SWATH windows.");
      }

      std::vector<std::pair<double, double>> uis_swathes;
      const int num_precursor_windows = static_cast<int>(Math::round(
        (parameters.precursor_upper_mz_limit - parameters.precursor_lower_mz_limit) /
        parameters.precursor_mz_threshold));
      uis_swathes.reserve(num_precursor_windows);
      for (int i = 0; i < num_precursor_windows; ++i)
      {
        uis_swathes.emplace_back(
          parameters.precursor_lower_mz_limit + (i * parameters.precursor_mz_threshold),
          parameters.precursor_lower_mz_limit + ((i + 1) * parameters.precursor_mz_threshold));
      }
      return uis_swathes;
    }

    void prepareGlobalAssays_(
      OpenSwath::LightTargetedExperiment& experiment,
      const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters,
      const ProgressLogger::LogType log_type)
    {
      std::vector<std::pair<double, double>> uis_swathes;
      if (parameters.enable_ipf)
      {
        // Build/validate fallback UIS windows before transition processing so a
        // non-positive precursor threshold fails deterministically before use.
        uis_swathes = buildPredictedUISSwathes_(parameters);
      }

      prepareDetectionTransitions_(experiment, parameters, log_type);
      if (!parameters.enable_ipf)
      {
        return;
      }

      int uis_seed = parameters.ipf_decoy_seed;
      bool disable_decoy_transitions = parameters.disable_decoy_transitions;
      if (parameters.test_mode)
      {
        if (uis_seed == -1)
        {
          uis_seed = 42;
        }
        disable_decoy_transitions = true;
      }

      MRMAssay assays;
      assays.setLogType(log_type);
      assays.uisTransitionsLight(
        experiment,
        parameters.allowed_fragment_types,
        parameters.allowed_fragment_charges,
        parameters.enable_identification_specific_losses,
        parameters.enable_identification_unspecific_losses,
        parameters.enable_identification_ms2_precursors,
        parameters.product_mz_threshold,
        uis_swathes,
        -4,
        parameters.max_num_alternative_localizations,
        uis_seed,
        disable_decoy_transitions);
      assays.restrictTransitionsLight(
        experiment,
        parameters.product_lower_mz_limit,
        parameters.product_upper_mz_limit,
        {});
    }

    void validatePredictedDecoyCoverage_(
      const OpenSwath::LightTargetedExperiment& target,
      const OpenSwath::LightTargetedExperiment& decoy,
      const double min_decoy_fraction)
    {
      if (target.compounds.empty() || target.proteins.empty())
      {
        throw Exception::IllegalArgument(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Predicted library preparation produced no target compounds or proteins.");
      }

      const double compound_fraction =
        static_cast<double>(decoy.compounds.size()) / static_cast<double>(target.compounds.size());
      const double protein_fraction =
        static_cast<double>(decoy.proteins.size()) / static_cast<double>(target.proteins.size());
      if (compound_fraction < min_decoy_fraction || protein_fraction < min_decoy_fraction)
      {
        throw Exception::IllegalArgument(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Predicted library decoy generation did not reach the configured minimum decoy fraction.");
      }
    }

    OpenSwathLibraryPreparation::LibraryStats collectPredictedStats_(
      const OpenSwath::LightTargetedExperiment& experiment)
    {
      OpenSwathLibraryPreparation::LibraryStats stats;
      stats.protein_count = experiment.proteins.size();
      stats.compound_count = experiment.compounds.size();
      stats.transition_count = experiment.transitions.size();
      for (const auto& transition : experiment.transitions)
      {
        if (transition.getDecoy())
        {
          ++stats.decoy_transition_count;
        }
        if (transition.isIdentifyingTransition())
        {
          ++stats.identifying_transition_count;
        }
      }
      return stats;
    }
  } // namespace
#endif

  OpenSwathLibraryPreparation::LibraryStats OpenSwathLibraryPreparation::preparePredictedLibraryToPQP(
    const std::string& input_fasta,
    const std::string& output_pqp,
    const AssayGeneratorParameters& assay_parameters,
    const DecoyGeneratorParameters& decoy_parameters,
    const PredictedLibraryParameters& parameters) const
  {
#ifndef WITH_ONNX
    (void)input_fasta;
    (void)output_pqp;
    (void)assay_parameters;
    (void)decoy_parameters;
    (void)parameters;
    throw Exception::Precondition(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Predicted library preparation requires an OpenMS build configured with WITH_ONNX=ON.");
#else
    if (parameters.max_peptide_length < parameters.min_peptide_length)
    {
      throw Exception::InvalidParameter(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PredictedLibraryParameters::max_peptide_length must be >= min_peptide_length.");
    }
    if (parameters.precursor_charges.empty())
    {
      throw Exception::InvalidParameter(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PredictedLibraryParameters::precursor_charges cannot be empty.");
    }
    for (const Int charge : parameters.precursor_charges)
    {
      if (charge <= 0)
      {
        throw Exception::InvalidParameter(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "PredictedLibraryParameters::precursor_charges must contain only positive integers.");
      }
    }
    if (parameters.prediction_batch_size == 0)
    {
      throw Exception::InvalidParameter(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PredictedLibraryParameters::prediction_batch_size must be greater than zero.");
    }
    if (parameters.inference_threads <= 0)
    {
      throw Exception::InvalidParameter(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PredictedLibraryParameters::inference_threads must be greater than zero.");
    }

    ensureUnimodLoaded(assay_parameters);

    const auto fixed_modifications =
      ModifiedPeptideGenerator::getModifications(parameters.fixed_modifications);
    const auto variable_modifications =
      ModifiedPeptideGenerator::getModifications(parameters.variable_modifications);

    const auto reject_protein_terminal_modifications =
      [&](const ModifiedPeptideGenerator::MapToResidueType& modifications)
    {
      for (const auto& [modification, residue] : modifications.val)
      {
        (void)residue;
        const auto specificity = modification->getTermSpecificity();
        if (specificity == ResidueModification::PROTEIN_N_TERM ||
            specificity == ResidueModification::PROTEIN_C_TERM)
        {
          throw Exception::InvalidParameter(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Predicted library preparation does not support protein-terminal modifications because "
            "the predictor operates on digested peptide sequences without protein-terminal position metadata.");
        }
      }
    };
    reject_protein_terminal_modifications(fixed_modifications);
    reject_protein_terminal_modifications(variable_modifications);

    // Construct all ONNX sessions before touching the FASTA so model/configuration
    // problems fail before potentially expensive digestion and candidate enumeration.
    PeptDeepLibraryPredictor::Config predictor_config;
    predictor_config.rt_model_path = parameters.rt_model_path;
    predictor_config.ccs_model_path = parameters.ccs_model_path;
    predictor_config.ms2_model_path = parameters.ms2_model_path;
    predictor_config.intra_op_threads = parameters.inference_threads;
    predictor_config.batch_size = parameters.prediction_batch_size;
    predictor_config.predict_ccs = parameters.predict_ccs;
    PeptDeepLibraryPredictor predictor(predictor_config);

    ProteaseDigestion digestion;
    digestion.setEnzyme(parameters.enzyme);
    digestion.setMissedCleavages(parameters.missed_cleavages);

    FASTAContainer<TFI_File> decoy_scan(input_fasta);
    const DecoyHelper::Result decoy = DecoyHelper::findDecoyString(decoy_scan, true);
    const auto is_decoy_accession = [&](const std::string& accession)
    {
      if (!decoy.success) return false;
      if (decoy.is_prefix) return StringUtils::hasPrefix(accession, decoy.name);
      return StringUtils::hasSuffix(accession, decoy.name);
    };

    std::map<std::string, std::set<std::string>> peptide_proteins;
    Size fasta_proteins = 0;
    Size skipped_decoy_proteins = 0;
    Size skipped_ambiguous_peptides = 0;

    const auto collect_digested = [&](const AASequence& protein, const std::string& protein_id)
    {
      std::vector<AASequence> digested_peptides;
      digestion.digest(
        protein,
        digested_peptides,
        parameters.min_peptide_length,
        parameters.max_peptide_length);
      for (const auto& peptide : digested_peptides)
      {
        const std::string sequence = peptide.toUnmodifiedString();
        if (sequence.find_first_not_of("ACDEFGHIKLMNPQRSTVWY") != std::string::npos)
        {
          ++skipped_ambiguous_peptides;
          continue;
        }
        peptide_proteins[sequence].insert(protein_id);
      }
    };

    FASTAFile fasta;
    fasta.readStart(input_fasta);
    FASTAFile::FASTAEntry entry;
    while (fasta.readNext(entry))
    {
      ++fasta_proteins;
      if (entry.identifier.empty())
      {
        throw Exception::InvalidValue(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, entry.identifier,
          "FASTA entries used for predicted libraries require a non-empty protein identifier.");
      }
      if (is_decoy_accession(entry.identifier))
      {
        ++skipped_decoy_proteins;
        continue;
      }

      const std::string protein_sequence = StringUtils::toUppered(entry.sequence);
      const AASequence protein = AASequence::fromString(protein_sequence);
      collect_digested(protein, entry.identifier);

      if (parameters.clip_nterm_methionine &&
          protein.size() > 1 &&
          !protein_sequence.empty() &&
          protein_sequence.front() == 'M')
      {
        collect_digested(
          protein.getSubsequence(1, static_cast<UInt>(protein.size() - 1)),
          entry.identifier);
      }
    }

    std::map<std::string, PredictedCandidate_> candidates;
    for (const auto& [sequence, protein_refs] : peptide_proteins)
    {
      AASequence peptide = AASequence::fromString(sequence, false);
      ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, peptide);

      std::vector<AASequence> peptidoforms;
      if (parameters.variable_modifications.empty() ||
          parameters.max_variable_modifications == 0)
      {
        peptidoforms.push_back(std::move(peptide));
      }
      else
      {
        ModifiedPeptideGenerator::applyVariableModifications(
          variable_modifications,
          peptide,
          parameters.max_variable_modifications,
          peptidoforms,
          true);
      }

      for (const auto& peptidoform : peptidoforms)
      {
        const std::string peptidoform_id = peptidoform.toUniModString();
        for (const Int charge : parameters.precursor_charges)
        {
          const double precursor_mz = peptidoform.getMZ(charge);
          if (precursor_mz < assay_parameters.precursor_lower_mz_limit ||
              precursor_mz > assay_parameters.precursor_upper_mz_limit)
          {
            continue;
          }

          const std::string id = peptidoform_id + "/" + std::to_string(charge);
          auto [candidate_it, inserted] = candidates.try_emplace(id);
          if (inserted)
          {
            candidate_it->second.peptide = peptidoform;
            candidate_it->second.charge = charge;
          }
          candidate_it->second.protein_refs.insert(
            protein_refs.begin(), protein_refs.end());
        }
      }
    }

    if (fasta_proteins == 0)
    {
      throw Exception::InvalidValue(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_fasta,
        "The FASTA input does not contain any protein entries.");
    }
    if (candidates.empty())
    {
      throw Exception::InvalidValue(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_fasta,
        "FASTA digestion/modification/charge enumeration produced no supported precursor candidates "
        "inside the configured precursor m/z range.");
    }

    if (decoy.success)
    {
      OPENMS_LOG_INFO << "Detected FASTA decoy "
                      << (decoy.is_prefix ? "prefix" : "suffix")
                      << " '" << decoy.name << "'; skipped "
                      << skipped_decoy_proteins << " decoy protein entr"
                      << (skipped_decoy_proteins == 1 ? "y" : "ies") << "."
                      << std::endl;
    }
    if (skipped_ambiguous_peptides > 0)
    {
      OPENMS_LOG_WARN << "Skipped " << skipped_ambiguous_peptides
                      << " digested peptide occurrence(s) containing residues outside "
                      << "ACDEFGHIKLMNPQRSTVWY." << std::endl;
    }
    OPENMS_LOG_INFO << "Predicted library candidates: " << candidates.size()
                    << " unique precursors from " << fasta_proteins
                    << " FASTA proteins." << std::endl;

    OpenSwath::LightTargetedExperiment predicted_library;
    std::unordered_set<std::string> protein_ids;
    std::vector<PeptDeepLibraryPrecursor> batch;
    batch.reserve(parameters.prediction_batch_size);

    const auto flush_batch = [&]()
    {
      if (batch.empty()) return;

      OpenSwath::LightTargetedExperiment predicted_batch = predictor.predict(batch);

      // Bound the materialized transition set before retaining the batch. Global
      // UIS/IPF preparation is applied once below after all retained batches have
      // been merged.
      prepareDetectionTransitions_(predicted_batch, assay_parameters, log_type_);

      predicted_library.compounds.insert(
        predicted_library.compounds.end(),
        std::make_move_iterator(predicted_batch.compounds.begin()),
        std::make_move_iterator(predicted_batch.compounds.end()));
      predicted_library.transitions.insert(
        predicted_library.transitions.end(),
        std::make_move_iterator(predicted_batch.transitions.begin()),
        std::make_move_iterator(predicted_batch.transitions.end()));
      for (auto& protein : predicted_batch.proteins)
      {
        if (protein_ids.insert(protein.id).second)
        {
          predicted_library.proteins.push_back(std::move(protein));
        }
      }
      batch.clear();
    };

    for (const auto& [id, candidate] : candidates)
    {
      PeptDeepLibraryPrecursor precursor;
      precursor.peptide = candidate.peptide;
      precursor.id = id;
      precursor.charge = candidate.charge;
      precursor.nce = static_cast<float>(parameters.nce);
      precursor.instrument_index = parameters.instrument_index;
      precursor.protein_refs.assign(
        candidate.protein_refs.begin(), candidate.protein_refs.end());
      batch.push_back(std::move(precursor));
      if (batch.size() == parameters.prediction_batch_size)
      {
        flush_batch();
      }
    }
    flush_batch();

    // Perform the final shared assay-preparation semantics in memory rather than
    // writing and re-reading an intermediate raw PQP.
    prepareGlobalAssays_(predicted_library, assay_parameters, log_type_);
    if (predicted_library.transitions.empty())
    {
      throw Exception::Precondition(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Predicted assay preparation produced zero transitions.");
    }

    OpenSwath::LightTargetedExperiment predicted_decoys;
    MRMDecoy decoys;
    decoys.setLogType(log_type_);
    decoys.generateDecoysLight(
      predicted_library,
      predicted_decoys,
      decoy_parameters.method,
      decoy_parameters.aim_decoy_fraction,
      decoy_parameters.switch_kr,
      decoy_parameters.decoy_tag,
      decoy_parameters.shuffle_max_attempts,
      decoy_parameters.shuffle_sequence_identity_threshold,
      decoy_parameters.shift_precursor_mz_shift,
      decoy_parameters.shift_product_mz_shift,
      decoy_parameters.product_mz_threshold,
      decoy_parameters.allowed_fragment_types,
      decoy_parameters.allowed_fragment_charges,
      decoy_parameters.enable_detection_specific_losses,
      decoy_parameters.enable_detection_unspecific_losses);

    validatePredictedDecoyCoverage_(
      predicted_library,
      predicted_decoys,
      decoy_parameters.min_decoy_fraction);

    OpenSwath::LightTargetedExperiment prepared_library;
    if (decoy_parameters.separate)
    {
      prepared_library = std::move(predicted_decoys);
    }
    else
    {
      prepared_library = std::move(predicted_library);
      prepared_library.transitions.insert(
        prepared_library.transitions.end(),
        std::make_move_iterator(predicted_decoys.transitions.begin()),
        std::make_move_iterator(predicted_decoys.transitions.end()));
      prepared_library.compounds.insert(
        prepared_library.compounds.end(),
        std::make_move_iterator(predicted_decoys.compounds.begin()),
        std::make_move_iterator(predicted_decoys.compounds.end()));
      prepared_library.proteins.insert(
        prepared_library.proteins.end(),
        std::make_move_iterator(predicted_decoys.proteins.begin()),
        std::make_move_iterator(predicted_decoys.proteins.end()));
    }

    const LibraryStats stats = collectPredictedStats_(prepared_library);
    if (!stats.hasDecoys())
    {
      throw Exception::Precondition(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Predicted library decoy generation produced zero decoy transitions.");
    }

    TransitionPQPFile writer;
    writer.setLogType(log_type_);
    const auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(prepared_library);
    writer.convertLightTargetedExperimentToPQP(
      output_pqp.c_str(), prepared_library, &source_ids);
    return stats;
#endif
  }
} // namespace OpenMS
