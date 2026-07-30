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

#include <cstdlib>
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
                           const TargetedExperiment& targeted_exp)
    {
      if (output_type == FileTypes::TSV)
      {
        TargetedExperiment writable_exp(targeted_exp);
        TransitionTSVFile tsv_writer;
        tsv_writer.setLogType(log_type);
        tsv_writer.convertTargetedExperimentToTSV(output_file.c_str(), writable_exp);
      }
      else if (output_type == FileTypes::PQP)
      {
        TargetedExperiment writable_exp(targeted_exp);
        TransitionPQPFile pqp_writer;
        pqp_writer.setLogType(log_type);
        pqp_writer.convertTargetedExperimentToPQP(output_file.c_str(), writable_exp);
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

    void ensureUnimodLoaded_(const OpenSwathLibraryPreparation::AssayGeneratorParameters& parameters)
    {
      if (!parameters.enable_ipf)
      {
        return;
      }

      if (parameters.unimod_file.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Please provide a valid Unimod XML file for IPF.");
      }

      if (!ModificationsDB::isInstantiated())
      {
        const ModificationsDB* ptr = ModificationsDB::initializeModificationsDB(parameters.unimod_file, std::string(""), std::string(""));
        OPENMS_LOG_INFO << "Unimod XML: " << ptr->getNumberOfModifications()
                        << " modification types and residue specificities imported from file: "
                        << parameters.unimod_file << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "ModificationsDB is already initialized; reusing the existing modification database for IPF.\n";
      }
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

    void validateDecoyCoverage_(const Size target_compounds,
                                const Size target_proteins,
                                const Size decoy_compounds,
                                const Size decoy_proteins,
                                const double min_decoy_fraction)
    {
      if (target_compounds == 0 || target_proteins == 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "The input experiment has no compounds or proteins.");
      }

      if (static_cast<double>(decoy_compounds) / static_cast<double>(target_compounds) < min_decoy_fraction ||
          static_cast<double>(decoy_proteins) / static_cast<double>(target_proteins) < min_decoy_fraction)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "The number of decoys for peptides or proteins is below the threshold of " +
                                         StringUtils::toStr(min_decoy_fraction * 100) + "% of the number of targets.");
      }
    }
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
      saveHeavyLibrary_(output_pqp, FileTypes::PQP, log_type_, targeted_exp);
      return collectStats_(targeted_exp);
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
    if (useLightPath_(input_type, output_type))
    {
      OpenSwath::LightTargetedExperiment light_exp;
      loadLightLibrary_(input_file, input_type, reader_parameters, log_type_, light_exp);

      MRMAssay assays;
      assays.setLogType(log_type_);
      assays.reannotateTransitionsLight(light_exp, parameters.precursor_mz_threshold, parameters.product_mz_threshold,
                                        parameters.allowed_fragment_types, parameters.allowed_fragment_charges,
                                        parameters.enable_detection_specific_losses, parameters.enable_detection_unspecific_losses);
      assays.restrictTransitionsLight(light_exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit,
                                      parameters.swathes);
      assays.detectingTransitionsLight(light_exp, parameters.min_transitions, parameters.max_transitions);

      if (parameters.enable_ipf)
      {
        ensureUnimodLoaded_(parameters);
        const auto [uis_seed, disable_decoy_transitions] = resolveIPFDecoySettings_(parameters);
        std::vector<std::pair<double, double>> uis_swathes = buildUISSwathes_(parameters);
        assays.uisTransitionsLight(light_exp,
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
        assays.restrictTransitionsLight(light_exp, parameters.product_lower_mz_limit, parameters.product_upper_mz_limit, {});
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
      ensureUnimodLoaded_(parameters);
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

    saveHeavyLibrary_(output_file, output_type, log_type_, targeted_exp);
    return collectStats_(targeted_exp);
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
      OpenSwath::LightTargetedExperiment light_merged;
      loadLightLibrary_(input_file, input_type, reader_parameters, log_type_, light_exp);

      MRMDecoy decoys;
      decoys.setLogType(log_type_);
      decoys.generateDecoysLight(light_exp, light_decoy, parameters.method,
                                 parameters.aim_decoy_fraction, parameters.switch_kr, parameters.decoy_tag,
                                 parameters.shuffle_max_attempts, parameters.shuffle_sequence_identity_threshold,
                                 parameters.shift_precursor_mz_shift, parameters.shift_product_mz_shift,
                                 parameters.product_mz_threshold, parameters.allowed_fragment_types,
                                 parameters.allowed_fragment_charges, parameters.enable_detection_specific_losses,
                                 parameters.enable_detection_unspecific_losses);

      validateDecoyCoverage_(light_exp.getCompounds().size(), light_exp.getProteins().size(),
                             light_decoy.getCompounds().size(), light_decoy.getProteins().size(),
                             parameters.min_decoy_fraction);

      if (parameters.separate)
      {
        light_merged = std::move(light_decoy);
      }
      else
      {
        light_merged = std::move(light_exp);
        light_merged.transitions.insert(light_merged.transitions.end(),
                                        light_decoy.transitions.begin(), light_decoy.transitions.end());
        light_merged.compounds.insert(light_merged.compounds.end(),
                                      light_decoy.compounds.begin(), light_decoy.compounds.end());
        light_merged.proteins.insert(light_merged.proteins.end(),
                                     light_decoy.proteins.begin(), light_decoy.proteins.end());
      }

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

    validateDecoyCoverage_(targeted_exp.getPeptides().size(), targeted_exp.getProteins().size(),
                           targeted_decoy.getPeptides().size(), targeted_decoy.getProteins().size(),
                           parameters.min_decoy_fraction);

    if (parameters.separate)
    {
      targeted_merged = std::move(targeted_decoy);
    }
    else
    {
      targeted_merged = std::move(targeted_exp);
      targeted_merged += std::move(targeted_decoy);
    }

    saveHeavyLibrary_(output_file, output_type, log_type_, targeted_merged);
    return collectStats_(targeted_merged);
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
    std::unique_ptr<File::TempDir> temp_dir;
    std::string working_dir = scratch_directory;
    if (working_dir.empty())
    {
      temp_dir = std::make_unique<File::TempDir>(true);
      working_dir = temp_dir->getPath();
    }
    else
    {
      File::makeDir(working_dir);
    }

    const std::string assay_output = File::absolutePath(working_dir + "/prepared_assays.pqp");
    const LibraryStats assay_stats = prepareAssays(input_file, input_type, assay_output, FileTypes::PQP, assay_parameters, reader_parameters);
    const auto removeAssayOutput = [&]()
    {
      if (assay_output != output_pqp && File::exists(assay_output))
      {
        File::remove(assay_output);
      }
    };

    if (assay_stats.transition_count == 0)
    {
      OPENMS_LOG_WARN << "Assay preparation produced an empty intermediate library; "
                      << "falling back to direct decoy generation on the empirical input library.\n";
      removeAssayOutput();
      return generateDecoys(input_file, input_type, output_pqp, FileTypes::PQP, decoy_parameters, reader_parameters);
    }

    try
    {
      const LibraryStats stats = generateDecoys(assay_output, FileTypes::PQP, output_pqp, FileTypes::PQP, decoy_parameters, reader_parameters);
      if (stats.hasDecoys())
      {
        removeAssayOutput();
        return stats;
      }

      OPENMS_LOG_WARN << "Decoy generation on the assay-prepared intermediate did not yield decoy transitions; "
                      << "falling back to direct decoy generation on the empirical input library.\n";
    }
    catch (const Exception::IllegalArgument& e)
    {
      OPENMS_LOG_WARN << "Decoy generation on the assay-prepared intermediate failed (" << e.what() << "); "
                      << "falling back to direct decoy generation on the empirical input library.\n";
    }

    removeAssayOutput();
    return generateDecoys(input_file, input_type, output_pqp, FileTypes::PQP, decoy_parameters, reader_parameters);
  }
} // namespace OpenMS
