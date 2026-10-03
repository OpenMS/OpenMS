// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <OpenMS/config.h>

#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <map>
#include <optional>
#include <set>
#include <string_view>
#include <unordered_set>
#endif

#include <cstdlib>
#include <iterator>
#include <mutex>
#include <utility>

namespace OpenMS
{
  namespace
  {
    bool useLightPath_(const FileTypes::Type input_type, const FileTypes::Type output_type)
    {
      const bool light_input = input_type == FileTypes::TSV || input_type == FileTypes::MRM ||
                               input_type == FileTypes::PQP || input_type == FileTypes::OSWPQ;
      const bool light_output = output_type == FileTypes::TSV || output_type == FileTypes::PQP ||
                                output_type == FileTypes::OSWPQ;
      return light_input && light_output;
    }

    OpenSwathLibraryPreparation::LibraryStats collectStats_(const OpenSwath::LightTargetedExperiment& exp)
    {
      OpenSwathLibraryPreparation::LibraryStats stats;
      stats.protein_count = exp.getProteins().size();
      stats.compound_count = exp.getCompounds().size();
      stats.transition_count = exp.getTransitions().size();
      for (const auto& transition : exp.getTransitions())
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

    OpenSwathLibraryPreparation::LibraryStats collectStats_(const TargetedExperiment& exp)
    {
      OpenSwathLibraryPreparation::LibraryStats stats;
      stats.protein_count = exp.getProteins().size();
      stats.compound_count = exp.getPeptides().size();
      stats.transition_count = exp.getTransitions().size();
      for (const auto& transition : exp.getTransitions())
      {
        if (transition.getDecoyTransitionType() == ReactionMonitoringTransition::DecoyTransitionType::DECOY)
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

    void loadLightLibrary_(const std::string& input_file,
                           const FileTypes::Type input_type,
                           const Param& reader_parameters,
                           const ProgressLogger::LogType log_type,
                           OpenSwath::LightTargetedExperiment& light_exp)
    {
      if (input_type == FileTypes::TSV || input_type == FileTypes::MRM)
      {
        TransitionTSVFile tsv_reader;
        tsv_reader.setLogType(log_type);
        tsv_reader.setParameters(reader_parameters);
        tsv_reader.convertTSVToTargetedExperiment(input_file.c_str(), input_type, light_exp);
      }
      else if (input_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_reader;
        pqp_reader.setLogType(log_type);
        pqp_reader.setParameters(reader_parameters);
        pqp_reader.convertPQPToTargetedExperiment(input_file.c_str(), light_exp);
      }
      else if (input_type == FileTypes::OSWPQ)
      {
        TransitionParquetFile parquet_reader;
        parquet_reader.convertParquetToTargetedExperiment(input_file, light_exp);
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported light-weight library input type '" +
                                          FileTypes::typeToName(input_type) + "'.");
      }
    }

    void loadHeavyLibrary_(const std::string& input_file,
                           const FileTypes::Type input_type,
                           const Param& reader_parameters,
                           const ProgressLogger::LogType log_type,
                           TargetedExperiment& targeted_exp)
    {
      if (input_type == FileTypes::TSV || input_type == FileTypes::MRM)
      {
        TransitionTSVFile tsv_reader;
        tsv_reader.setLogType(log_type);
        tsv_reader.setParameters(reader_parameters);
        tsv_reader.convertTSVToTargetedExperiment(input_file.c_str(), input_type, targeted_exp);
        tsv_reader.validateTargetedExperiment(targeted_exp);
      }
      else if (input_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_reader;
        pqp_reader.setLogType(log_type);
        pqp_reader.setParameters(reader_parameters);
        pqp_reader.convertPQPToTargetedExperiment(input_file.c_str(), targeted_exp);
        pqp_reader.validateTargetedExperiment(targeted_exp);
      }
      else if (input_type == FileTypes::TRAML)
      {
        FileHandler().loadTransitions(input_file, targeted_exp, {FileTypes::TRAML});
      }
      else if (input_type == FileTypes::OSWPQ)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Parquet input is only supported for light-weight library preparation paths.");
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported heavy library input type '" +
                                          FileTypes::typeToName(input_type) + "'.");
      }
    }

    void saveLightLibrary_(const std::string& output_file,
                           const FileTypes::Type output_type,
                           const ProgressLogger::LogType log_type,
                           const OpenSwath::LightTargetedExperiment& light_exp)
    {
      if (output_type == FileTypes::TSV)
      {
        TransitionTSVFile tsv_writer;
        tsv_writer.setLogType(log_type);
        tsv_writer.convertLightTargetedExperimentToTSV(output_file.c_str(), light_exp);
      }
      else if (output_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_writer;
        pqp_writer.setLogType(log_type);
        pqp_writer.convertLightTargetedExperimentToPQP(output_file.c_str(), light_exp);
      }
      else if (output_type == FileTypes::OSWPQ)
      {
        TransitionParquetFile parquet_writer;
        parquet_writer.convertLightTargetedExperimentToParquet(output_file, light_exp);
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported light-weight library output type '" +
                                          FileTypes::typeToName(output_type) + "'.");
      }
    }

    void saveHeavyLibrary_(const std::string& output_file,
                           const FileTypes::Type output_type,
                           const ProgressLogger::LogType log_type,
                           TargetedExperiment& targeted_exp)
    {
      if (output_type == FileTypes::TSV)
      {
        TransitionTSVFile tsv_writer;
        tsv_writer.setLogType(log_type);
        tsv_writer.convertTargetedExperimentToTSV(output_file.c_str(), targeted_exp);
      }
      else if (output_type == FileTypes::PQP)
      {
        TransitionPQPFile pqp_writer;
        pqp_writer.setLogType(log_type);
        pqp_writer.convertTargetedExperimentToPQP(output_file.c_str(), targeted_exp);
      }
      else if (output_type == FileTypes::TRAML)
      {
        FileHandler().storeTransitions(output_file, targeted_exp, {FileTypes::TRAML});
      }
      else if (output_type == FileTypes::OSWPQ)
      {
        OpenSwath::LightTargetedExperiment light_exp;
        OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
        TransitionParquetFile parquet_writer;
        parquet_writer.convertLightTargetedExperimentToParquet(output_file, light_exp);
      }
      else
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Unsupported heavy library output type '" +
                                          FileTypes::typeToName(output_type) + "'.");
      }
    }

    std::vector<std::pair<double, double>> buildUISSwathes_(const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters)
    {
      if (parameters.enable_swath_specifity && !parameters.swathes.empty())
      {
        return parameters.swathes;
      }

      if (parameters.precursor_mz_threshold <= 0.0)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "AssayGeneratorParameters::precursor_mz_threshold must be greater than zero "
                                          "when constructing fallback UIS SWATH windows.");
      }
      if (parameters.precursor_upper_mz_limit < parameters.precursor_lower_mz_limit)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "AssayGeneratorParameters::precursor_upper_mz_limit must not be below precursor_lower_mz_limit "
                                          "when constructing fallback UIS SWATH windows.");
      }

      std::vector<std::pair<double, double>> uis_swathes;
      const int num_precursor_windows = static_cast<int>(Math::round(
        (parameters.precursor_upper_mz_limit - parameters.precursor_lower_mz_limit) /
        parameters.precursor_mz_threshold));
      uis_swathes.reserve(num_precursor_windows);
      for (int i = 0; i < num_precursor_windows; ++i)
      {
        uis_swathes.emplace_back(parameters.precursor_lower_mz_limit + (i * parameters.precursor_mz_threshold),
                                 parameters.precursor_lower_mz_limit + ((i + 1) * parameters.precursor_mz_threshold));
      }
      return uis_swathes;
    }

    void prepareDetectionTransitionsLight_(OpenSwath::LightTargetedExperiment& exp,
                                           const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters,
                                           const ProgressLogger::LogType log_type)
    {
      MRMAssay assays;
      assays.setLogType(log_type);
      assays.reannotateTransitionsLight(exp, parameters.precursor_mz_threshold, parameters.product_mz_threshold,
                                        parameters.allowed_fragment_types, parameters.allowed_fragment_charges,
                                        parameters.enable_detection_specific_losses, parameters.enable_detection_unspecific_losses);
      assays.restrictTransitionsLight(exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit,
                                      parameters.swathes);
      assays.detectingTransitionsLight(exp, parameters.min_transitions, parameters.max_transitions);
    }

    std::pair<int, bool> resolveIPFDecoySettings_(const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters)
    {
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
      return {uis_seed, disable_decoy_transitions};
    }

    void addUISTransitionsLight_(OpenSwath::LightTargetedExperiment& exp,
                                 const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters,
                                 const std::vector<std::pair<double, double>>& uis_swathes,
                                 const ProgressLogger::LogType log_type)
    {
      const auto [uis_seed, disable_decoy_transitions] = resolveIPFDecoySettings_(parameters);
      MRMAssay assays;
      assays.setLogType(log_type);
      assays.uisTransitionsLight(exp,
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
      assays.restrictTransitionsLight(exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit, {});
    }

    void generateDecoysLight_(const OpenSwath::LightTargetedExperiment& targets,
                              OpenSwath::LightTargetedExperiment& decoys,
                              const OpenSwathLibraryPreparation::DecoyGeneratorParameters& parameters,
                              const ProgressLogger::LogType log_type)
    {
      MRMDecoy generator;
      generator.setLogType(log_type);
      generator.generateDecoysLight(targets, decoys, parameters.method,
                                    parameters.aim_decoy_fraction, parameters.switch_kr, parameters.decoy_tag,
                                    parameters.shuffle_max_attempts, parameters.shuffle_sequence_identity_threshold,
                                    parameters.shift_precursor_mz_shift, parameters.shift_product_mz_shift,
                                    parameters.product_mz_threshold, parameters.allowed_fragment_types,
                                    parameters.allowed_fragment_charges, parameters.enable_detection_specific_losses,
                                    parameters.enable_detection_unspecific_losses);
    }

    /// Target+decoy library, or only the decoys when @p separate is set.
    OpenSwath::LightTargetedExperiment mergeTargetDecoyLight_(OpenSwath::LightTargetedExperiment&& targets,
                                                              OpenSwath::LightTargetedExperiment&& decoys,
                                                              const bool separate)
    {
      if (separate)
      {
        // Only the decoys are written: release the targets now rather than when the caller returns.
        targets = OpenSwath::LightTargetedExperiment();
        return std::move(decoys);
      }
      OpenSwath::LightTargetedExperiment merged = std::move(targets);
      merged.transitions.insert(merged.transitions.end(),
                                std::make_move_iterator(decoys.transitions.begin()), std::make_move_iterator(decoys.transitions.end()));
      merged.compounds.insert(merged.compounds.end(),
                              std::make_move_iterator(decoys.compounds.begin()), std::make_move_iterator(decoys.compounds.end()));
      merged.proteins.insert(merged.proteins.end(),
                             std::make_move_iterator(decoys.proteins.begin()), std::make_move_iterator(decoys.proteins.end()));
      return merged;
    }

    void ensureUnimodLoaded_(const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters)
    {
      if (parameters.enable_ipf && parameters.unimod_file.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Please provide a valid Unimod XML file for IPF.");
      }

      // Preserve the historical TOPP behavior: an explicitly supplied Unimod file
      // is loaded even when IPF is disabled. With no requested file there is
      // nothing for this helper to initialize.
      if (parameters.unimod_file.empty())
      {
        return;
      }

      // ModificationsDB is a process-global singleton. Serialize both the
      // check/initialization sequence and the helper-owned provenance marker so
      // concurrent library preparation cannot race here.
      static std::mutex unimod_init_mutex;
      static std::string helper_initialized_unimod;
      const std::lock_guard<std::mutex> lock(unimod_init_mutex);

      const std::string requested_unimod = File::absolutePath(parameters.unimod_file);
      if (!ModificationsDB::isInstantiated())
      {
        const ModificationsDB* ptr = ModificationsDB::initializeModificationsDB(parameters.unimod_file, std::string(""), std::string(""));
        helper_initialized_unimod = requested_unimod;
        OPENMS_LOG_INFO << "Unimod XML: " << ptr->getNumberOfModifications()
                        << " modification types and residue specificities imported from file: "
                        << parameters.unimod_file << std::endl;
      }
      else if (!helper_initialized_unimod.empty() && helper_initialized_unimod == requested_unimod)
      {
        return;
      }
      else if (parameters.reuse_existing_modifications_db)
      {
        OPENMS_LOG_INFO << "ModificationsDB is already initialized; explicitly reusing the existing modification database for IPF.\n";
      }
      else
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "ModificationsDB was initialized before the configured Unimod XML could be loaded. "
                                      "Initialize the requested Unimod source through OpenSwathLibraryPreparation before parsing modified sequences, "
                                      "or explicitly opt in to reusing the existing database.");
      }
    }

    void validateDecoyCoverage_(const Size target_compounds,
                                const Size target_proteins,
                                const Size decoy_compounds,
                                const Size decoy_proteins,
                                const double min_decoy_fraction,
                                const bool require_compounds,
                                const bool require_proteins)
    {
      if ((require_compounds && target_compounds == 0) ||
          (require_proteins && target_proteins == 0))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "The input experiment has no compounds or proteins.");
      }

      // The heavy TraML path historically permits peptide-only/protein-less and
      // compound-only inputs. Only validate entity classes that are actually
      // populated; small-molecule TargetedExperiment::Compound entries are not
      // peptide-decoy targets and are therefore excluded by the heavy caller.
      const bool compounds_low = target_compounds > 0 &&
        static_cast<double>(decoy_compounds) / static_cast<double>(target_compounds) < min_decoy_fraction;
      const bool proteins_low = target_proteins > 0 &&
        static_cast<double>(decoy_proteins) / static_cast<double>(target_proteins) < min_decoy_fraction;
      if (compounds_low || proteins_low)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "The number of decoys for compounds/peptides or proteins is below the threshold of " +
                                         StringUtils::toStr(min_decoy_fraction * 100) + "% of the number of targets.");
      }
    }

    /// Removes an intermediate file when leaving scope, on success and on error.
    struct ScratchFileRemover_
    {
      std::string path;

      ~ScratchFileRemover_()
      {
        try
        {
          File::remove(path);
        }
        catch (...)
        {
          // Never throw from a destructor; a leftover scratch file is harmless.
        }
      }
    };
  } // namespace

  OpenSwathLibraryPreparation::OpenSwathLibraryPreparation() = default;

  void OpenSwathLibraryPreparation::setLogType(const ProgressLogger::LogType log_type)
  {
    log_type_ = log_type;
  }

  ProgressLogger::LogType OpenSwathLibraryPreparation::getLogType() const
  {
    return log_type_;
  }

  void OpenSwathLibraryPreparation::ensureUnimodLoaded(const AssayGeneratorParameters& parameters) const
  {
    ensureUnimodLoaded_(parameters);
  }

  OpenSwathLibraryPreparation::LibraryStats OpenSwathLibraryPreparation::normalizeLibraryToPQP(
    const std::string& input_file,
    FileTypes::Type input_type,
    const std::string& output_pqp,
    const Param& reader_parameters) const
  {
    if (input_type == FileTypes::UNKNOWN)
    {
      input_type = FileHandler().getType(input_file);
    }

    if (input_type == FileTypes::UNKNOWN)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Could not determine input library file type.");
    }

    if (input_type == FileTypes::TRAML)
    {
      TargetedExperiment targeted_exp;
      loadHeavyLibrary_(input_file, input_type, reader_parameters, log_type_, targeted_exp);
      const LibraryStats stats = collectStats_(targeted_exp);
      saveHeavyLibrary_(output_pqp, FileTypes::PQP, log_type_, targeted_exp);
      return stats;
    }

    OpenSwath::LightTargetedExperiment light_exp;
    loadLightLibrary_(input_file, input_type, reader_parameters, log_type_, light_exp);
    saveLightLibrary_(output_pqp, FileTypes::PQP, log_type_, light_exp);
    return collectStats_(light_exp);
  }

  OpenSwathLibraryPreparation::LibraryStats OpenSwathLibraryPreparation::prepareAssays(
    const std::string& input_file,
    FileTypes::Type input_type,
    const std::string& output_file,
    FileTypes::Type output_type,
    const AssayGeneratorParameters& parameters,
    const Param& reader_parameters) const
  {
    ensureUnimodLoaded(parameters);

    if (useLightPath_(input_type, output_type))
    {
      OpenSwath::LightTargetedExperiment light_exp;
      loadLightLibrary_(input_file, input_type, reader_parameters, log_type_, light_exp);

      prepareDetectionTransitionsLight_(light_exp, parameters, log_type_);
      if (parameters.enable_ipf)
      {
        addUISTransitionsLight_(light_exp, parameters, buildUISSwathes_(parameters), log_type_);
      }

      saveLightLibrary_(output_file, output_type, log_type_, light_exp);
      return collectStats_(light_exp);
    }

    TargetedExperiment targeted_exp;
    loadHeavyLibrary_(input_file, input_type, reader_parameters, log_type_, targeted_exp);

    MRMAssay assays;
    assays.setLogType(log_type_);
    assays.reannotateTransitions(targeted_exp, parameters.precursor_mz_threshold, parameters.product_mz_threshold,
                                 parameters.allowed_fragment_types, parameters.allowed_fragment_charges,
                                 parameters.enable_detection_specific_losses, parameters.enable_detection_unspecific_losses);
    assays.restrictTransitions(targeted_exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit,
                               parameters.swathes);
    assays.detectingTransitions(targeted_exp, parameters.min_transitions, parameters.max_transitions);

    if (parameters.enable_ipf)
    {
      const auto [uis_seed, disable_decoy_transitions] = resolveIPFDecoySettings_(parameters);
      std::vector<std::pair<double, double>> uis_swathes = buildUISSwathes_(parameters);
      assays.uisTransitions(targeted_exp,
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
      assays.restrictTransitions(targeted_exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit, {});
    }

    const LibraryStats stats = collectStats_(targeted_exp);
    saveHeavyLibrary_(output_file, output_type, log_type_, targeted_exp);
    return stats;
  }

  OpenSwathLibraryPreparation::LibraryStats OpenSwathLibraryPreparation::generateDecoys(
    const std::string& input_file,
    FileTypes::Type input_type,
    const std::string& output_file,
    FileTypes::Type output_type,
    const DecoyGeneratorParameters& parameters,
    const Param& reader_parameters) const
  {
    if (useLightPath_(input_type, output_type))
    {
      OpenSwath::LightTargetedExperiment light_exp;
      OpenSwath::LightTargetedExperiment light_decoy;
      loadLightLibrary_(input_file, input_type, reader_parameters, log_type_, light_exp);

      generateDecoysLight_(light_exp, light_decoy, parameters, log_type_);

      // Preserve the historical light-path contract: both compounds and
      // protein annotations are required for TSV/PQP/OpenSwath light libraries.
      validateDecoyCoverage_(light_exp.getCompounds().size(), light_exp.getProteins().size(),
                             light_decoy.getCompounds().size(), light_decoy.getProteins().size(),
                             parameters.min_decoy_fraction, true, true);

      const OpenSwath::LightTargetedExperiment light_merged =
        mergeTargetDecoyLight_(std::move(light_exp), std::move(light_decoy), parameters.separate);

      saveLightLibrary_(output_file, output_type, log_type_, light_merged);
      return collectStats_(light_merged);
    }

    TargetedExperiment targeted_exp;
    TargetedExperiment targeted_decoy;
    TargetedExperiment targeted_merged;
    loadHeavyLibrary_(input_file, input_type, reader_parameters, log_type_, targeted_exp);

    MRMDecoy decoys;
    decoys.setLogType(log_type_);
    decoys.generateDecoys(targeted_exp, targeted_decoy, parameters.method,
                          parameters.aim_decoy_fraction, parameters.switch_kr, parameters.decoy_tag,
                          parameters.shuffle_max_attempts, parameters.shuffle_sequence_identity_threshold,
                          parameters.shift_precursor_mz_shift, parameters.shift_product_mz_shift,
                          parameters.product_mz_threshold, parameters.allowed_fragment_types,
                          parameters.allowed_fragment_charges, parameters.enable_detection_specific_losses,
                          parameters.enable_detection_unspecific_losses);

    // MRMDecoy generates peptide/protein decoys, not small-molecule compound
    // decoys. Do not count TargetedExperiment::Compound entries in this ratio.
    validateDecoyCoverage_(targeted_exp.getPeptides().size(),
                           targeted_exp.getProteins().size(),
                           targeted_decoy.getPeptides().size(),
                           targeted_decoy.getProteins().size(),
                           parameters.min_decoy_fraction, false, false);

    if (parameters.separate)
    {
      targeted_merged = std::move(targeted_decoy);
    }
    else
    {
      targeted_merged = std::move(targeted_exp);
      targeted_merged += std::move(targeted_decoy);
    }

    const LibraryStats stats = collectStats_(targeted_merged);
    saveHeavyLibrary_(output_file, output_type, log_type_, targeted_merged);
    return stats;
  }

  OpenSwathLibraryPreparation::LibraryStats OpenSwathLibraryPreparation::prepareEmpiricalLibraryToPQP(
    const std::string& input_file,
    FileTypes::Type input_type,
    const std::string& output_pqp,
    const AssayGeneratorParameters& assay_parameters,
    const DecoyGeneratorParameters& decoy_parameters,
    const Param& reader_parameters,
    const std::string& scratch_directory) const
  {
    // Always isolate the assay-preparation intermediate per invocation. A supplied
    // scratch_directory is only a (possibly shared) parent, so write a uniquely named
    // file into it rather than nesting another TempDir: callers such as OpenDIA already
    // pass a unique TempDir here, and a second ~90-character TempDir segment pushes the
    // path past the Windows MAX_PATH limit.
    std::unique_ptr<TempDir> temp_dir;
    std::string scratch_base;
    if (scratch_directory.empty())
    {
      temp_dir = std::make_unique<TempDir>();
      scratch_base = temp_dir->getPath();
    }
    else
    {
      scratch_base = File::absolutePath(scratch_directory);
      File::makeDir(scratch_base);
    }
    StringUtils::ensureLastChar(scratch_base, '/');

    std::string assay_output;
    do
    {
      assay_output = scratch_base + "prepared_assays_" + File::getUniqueName(false) + ".pqp";
    } while (File::exists(assay_output));
    const ScratchFileRemover_ assay_output_remover{assay_output};

    const LibraryStats assay_stats = prepareAssays(input_file, input_type, assay_output, FileTypes::PQP, assay_parameters, reader_parameters);

    if (assay_stats.transition_count == 0)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Assay preparation produced zero transitions. Refusing to bypass the configured assay filters by falling back to the raw empirical library.");
    }

    const LibraryStats stats = generateDecoys(assay_output, FileTypes::PQP, output_pqp, FileTypes::PQP, decoy_parameters, reader_parameters);
    if (!stats.hasDecoys())
    {
      // Fail closed: generateDecoys has already written output_pqp, so do not leave a
      // target-only library at the caller's output path.
      File::remove(output_pqp);
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Decoy generation on the assay-prepared library produced zero decoy transitions.");
    }
    return stats;
  }

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
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Predicted library preparation requires an OpenMS build configured with WITH_ONNX=ON.");
#else
    const auto invalid = [](const std::string& message)
    {
      return Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "PredictedLibraryParameters::" + message);
    };
    if (parameters.min_peptide_length < 2) throw invalid("min_peptide_length must be at least 2.");
    if (parameters.max_peptide_length < parameters.min_peptide_length) throw invalid("max_peptide_length must be >= min_peptide_length.");
    if (parameters.precursor_charges.empty()) throw invalid("precursor_charges cannot be empty.");
    for (const Int charge : parameters.precursor_charges)
    {
      if (charge <= 0) throw invalid("precursor_charges must contain only positive integers.");
    }
    if (parameters.prediction_batch_size == 0) throw invalid("prediction_batch_size must be greater than zero.");
    if (parameters.inference_threads <= 0) throw invalid("inference_threads must be greater than zero.");
    // PeptDeep embeds the instrument in a table of 8 slots (the last one means "unknown").
    if (parameters.instrument_index < 0 || parameters.instrument_index > 7) throw invalid("instrument_index must be between 0 and 7.");
    if (!ProteaseDB::getInstance()->hasEnzyme(parameters.enzyme)) throw invalid("enzyme '" + parameters.enzyme + "' is not a known protease.");

    ensureUnimodLoaded(assay_parameters);
    std::vector<std::pair<double, double>> uis_swathes;
    if (assay_parameters.enable_ipf)
    {
      uis_swathes = buildUISSwathes_(assay_parameters);
    }

    const auto fixed_modifications = ModifiedPeptideGenerator::getModifications(parameters.fixed_modifications);
    const auto variable_modifications = ModifiedPeptideGenerator::getModifications(parameters.variable_modifications);
    for (const auto* modifications : {&fixed_modifications, &variable_modifications})
    {
      for (const auto& [modification, residue] : modifications->val)
      {
        (void)residue;
        const auto specificity = modification->getTermSpecificity();
        if (specificity == ResidueModification::PROTEIN_N_TERM || specificity == ResidueModification::PROTEIN_C_TERM)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Predicted library preparation does not support protein-terminal modifications because "
                                            "the predictor operates on digested peptide sequences without protein-terminal position metadata.");
        }
      }
    }

    // Construct all ONNX sessions before touching the FASTA so model/configuration
    // problems fail before potentially expensive digestion and candidate enumeration.
    PeptDeepLibraryPredictor::Config predictor_config;
    predictor_config.rt_model_path = parameters.rt_model_path;
    predictor_config.ccs_model_path = parameters.ccs_model_path;
    predictor_config.ms2_model_path = parameters.ms2_model_path;
    predictor_config.intra_op_threads = parameters.inference_threads;
    predictor_config.batch_size = parameters.prediction_batch_size;
    predictor_config.predict_ccs = parameters.predict_ccs;
    std::optional<PeptDeepLibraryPredictor> predictor(std::in_place, predictor_config);

    ProteaseDigestion digestion;
    digestion.setEnzyme(parameters.enzyme);
    digestion.setMissedCleavages(parameters.missed_cleavages);

    // Skip decoy entries of target-decoy FASTAs. DecoyHelper only reports an affix
    // when it marks a large share of the database, so entries carrying the
    // configured decoy tag are always skipped as well.
    DecoyHelper::Result decoy;
    {
      FASTAContainer<TFI_File> decoy_scan(input_fasta);
      decoy = DecoyHelper::findDecoyString(decoy_scan, true);
    }
    const auto is_decoy_accession = [&](const std::string& accession)
    {
      if (!decoy_parameters.decoy_tag.empty() && StringUtils::hasPrefix(accession, decoy_parameters.decoy_tag))
      {
        return true;
      }
      if (!decoy.success) return false;
      return decoy.is_prefix ? StringUtils::hasPrefix(accession, decoy.name) : StringUtils::hasSuffix(accession, decoy.name);
    };

    std::map<std::string, std::set<std::string>> peptide_proteins;
    Size fasta_proteins = 0;
    Size skipped_decoy_proteins = 0;
    Size skipped_ambiguous_peptides = 0;
    std::vector<std::string_view> digested_peptides;
    const auto collect_digested = [&](std::string_view protein, const std::string& protein_id, const bool n_terminal_only)
    {
      digestion.digestUnmodified(protein, digested_peptides, parameters.min_peptide_length, parameters.max_peptide_length);
      for (const auto& peptide : digested_peptides)
      {
        if (n_terminal_only && peptide.data() != protein.data()) continue;
        if (peptide.find_first_not_of("ACDEFGHIKLMNPQRSTVWY") != std::string_view::npos)
        {
          ++skipped_ambiguous_peptides;
          continue;
        }
        peptide_proteins[std::string(peptide)].insert(protein_id);
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
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, entry.identifier,
                                      "FASTA entries used for predicted libraries require a non-empty protein identifier.");
      }
      if (is_decoy_accession(entry.identifier))
      {
        ++skipped_decoy_proteins;
        continue;
      }

      const std::string protein = StringUtils::toUppered(entry.sequence);
      collect_digested(protein, entry.identifier, false);
      // Only peptides at the new N-terminus differ after initiator-Met removal.
      if (parameters.clip_nterm_methionine && protein.size() > 1 && protein.front() == 'M')
      {
        collect_digested(std::string_view(protein).substr(1), entry.identifier, true);
      }
    }
    if (fasta_proteins == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_fasta,
                                    "The FASTA input does not contain any protein entries.");
    }

    if (decoy.success)
    {
      OPENMS_LOG_INFO << "Detected FASTA decoy " << (decoy.is_prefix ? "prefix" : "suffix")
                      << " '" << decoy.name << "'; skipped " << skipped_decoy_proteins << " decoy protein entr"
                      << (skipped_decoy_proteins == 1 ? "y" : "ies") << "." << std::endl;
    }
    else if (skipped_decoy_proteins > 0)
    {
      OPENMS_LOG_WARN << "Skipped " << skipped_decoy_proteins << " FASTA protein entr"
                      << (skipped_decoy_proteins == 1 ? "y" : "ies")
                      << " carrying the decoy tag '" << decoy_parameters.decoy_tag
                      << "', although no database-wide decoy prefix/suffix was detected." << std::endl;
    }
    if (skipped_ambiguous_peptides > 0)
    {
      OPENMS_LOG_WARN << "Skipped " << skipped_ambiguous_peptides
                      << " digested peptide occurrence(s) containing residues outside ACDEFGHIKLMNPQRSTVWY." << std::endl;
    }

    // Each precursor ID derives from exactly one unmodified sequence, so precursors
    // can be streamed into prediction batches without a global candidate map.
    OpenSwath::LightTargetedExperiment predicted_library;
    std::unordered_set<std::string> protein_ids;
    std::vector<PeptDeepLibraryPrecursor> batch;
    batch.reserve(parameters.prediction_batch_size);
    Size precursor_count = 0;

    const auto flush_batch = [&]()
    {
      if (batch.empty()) return;
      OpenSwath::LightTargetedExperiment predicted_batch = predictor->predict(batch);
      batch.clear();

      // Bound the materialized transition set before retaining the batch. These
      // steps act per precursor, so running them per batch equals a global pass.
      prepareDetectionTransitionsLight_(predicted_batch, assay_parameters, ProgressLogger::NONE);

      predicted_library.compounds.insert(predicted_library.compounds.end(),
                                         std::make_move_iterator(predicted_batch.compounds.begin()),
                                         std::make_move_iterator(predicted_batch.compounds.end()));
      predicted_library.transitions.insert(predicted_library.transitions.end(),
                                           std::make_move_iterator(predicted_batch.transitions.begin()),
                                           std::make_move_iterator(predicted_batch.transitions.end()));
      for (auto& protein : predicted_batch.proteins)
      {
        if (protein_ids.insert(protein.id).second)
        {
          predicted_library.proteins.push_back(std::move(protein));
        }
      }
    };

    ProgressLogger progress;
    progress.setLogType(log_type_);
    progress.startProgress(0, peptide_proteins.size(), "Predicting library from FASTA");
    Size progress_index = 0;
    std::vector<AASequence> peptidoforms;
    for (auto& [sequence, protein_refs] : peptide_proteins)
    {
      progress.setProgress(progress_index++);
      AASequence peptide = AASequence::fromString(sequence, false);
      ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, peptide);

      peptidoforms.clear();
      if (parameters.variable_modifications.empty() || parameters.max_variable_modifications == 0)
      {
        peptidoforms.push_back(std::move(peptide));
      }
      else
      {
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, peptide,
                                                             parameters.max_variable_modifications, peptidoforms, true);
      }

      const std::vector<std::string> refs(protein_refs.begin(), protein_refs.end());
      protein_refs.clear();
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

          PeptDeepLibraryPrecursor precursor;
          precursor.peptide = peptidoform;
          precursor.id = peptidoform_id + "/" + std::to_string(charge);
          precursor.charge = charge;
          precursor.nce = static_cast<float>(parameters.nce);
          precursor.instrument_index = parameters.instrument_index;
          precursor.protein_refs = refs;
          batch.push_back(std::move(precursor));
          ++precursor_count;
          if (batch.size() == parameters.prediction_batch_size)
          {
            flush_batch();
          }
        }
      }
    }
    flush_batch();
    progress.endProgress();

    // The global stages below need neither the FASTA digest nor the ONNX sessions. Release
    // them before the library is duplicated into decoys and normalized.
    peptide_proteins.clear();
    protein_ids = std::unordered_set<std::string>();
    batch = std::vector<PeptDeepLibraryPrecursor>();
    predictor.reset();

    if (precursor_count == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_fasta,
                                    "FASTA digestion/modification/charge enumeration produced no supported precursor candidates "
                                    "inside the configured precursor m/z range.");
    }
    OPENMS_LOG_INFO << "Predicted library candidates: " << precursor_count
                    << " unique precursors from " << fasta_proteins << " FASTA proteins." << std::endl;

    if (assay_parameters.enable_ipf)
    {
      addUISTransitionsLight_(predicted_library, assay_parameters, uis_swathes, log_type_);
    }
    if (predicted_library.transitions.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Predicted assay preparation produced zero transitions.");
    }

    OpenSwath::LightTargetedExperiment predicted_decoys;
    generateDecoysLight_(predicted_library, predicted_decoys, decoy_parameters, log_type_);
    validateDecoyCoverage_(predicted_library.compounds.size(), predicted_library.proteins.size(),
                           predicted_decoys.compounds.size(), predicted_decoys.proteins.size(),
                           decoy_parameters.min_decoy_fraction, true, true);

    OpenSwath::LightTargetedExperiment prepared_library =
      mergeTargetDecoyLight_(std::move(predicted_library), std::move(predicted_decoys), decoy_parameters.separate);
    const LibraryStats stats = collectStats_(prepared_library);
    if (!stats.hasDecoys())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Predicted library decoy generation produced zero decoy transitions.");
    }

    TransitionPQPFile writer;
    writer.setLogType(log_type_);
    const auto source_ids = OpenSwathLibraryIDNormalizer::normalizeSourceIDs(prepared_library);
    writer.convertLightTargetedExperimentToPQP(output_pqp.c_str(), prepared_library, &source_ids);
    return stats;
#endif
  }
} // namespace OpenMS
