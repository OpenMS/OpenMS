// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/config.h>

#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Shared spectral-library preparation helper for OpenSWATH-family workflows.

    This helper centralizes the file-format-aware spectral library preparation
    steps that historically lived directly inside TOPP front-ends:

      - prepared-library normalization to PQP
      - assay generation / preparation
      - decoy generation
      - empirical-library preparation for full DIA workflows

    The implementation mirrors the existing TOPP behavior, including the
    lightweight TSV/PQP/OSWPQ path where possible and the heavy TraML path when
    required.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathLibraryPreparation
  {
  public:
    struct OPENMS_DLLAPI LibraryStats
    {
      Size protein_count = 0;
      Size compound_count = 0;
      Size transition_count = 0;
      Size decoy_transition_count = 0;
      Size identifying_transition_count = 0;

      bool hasDecoys() const
      {
        return decoy_transition_count > 0;
      }

      bool hasIdentifyingTransitions() const
      {
        return identifying_transition_count > 0;
      }
    };

    struct OPENMS_DLLAPI AssayGeneratorParameters
    {
      Int min_transitions = 6;
      Int max_transitions = 6;
      std::vector<std::string> allowed_fragment_types{"b", "y"};
      std::vector<size_t> allowed_fragment_charges{1, 2, 3, 4};
      bool enable_detection_specific_losses = false;
      bool enable_detection_unspecific_losses = false;

      double precursor_mz_threshold = 0.025;
      double precursor_lower_mz_limit = 400.0;
      double precursor_upper_mz_limit = 1200.0;
      double product_mz_threshold = 0.025;
      double product_lower_mz_limit = 350.0;
      double product_upper_mz_limit = 2000.0;
      std::vector<std::pair<double, double>> swathes;

      bool enable_ipf = false;
      size_t max_num_alternative_localizations = 10000;
      bool enable_identification_ms2_precursors = true;
      bool enable_identification_specific_losses = true;
      bool enable_identification_unspecific_losses = false;
      bool enable_swath_specifity = false;
      bool disable_decoy_transitions = false;
      int ipf_decoy_seed = -1;
      bool test_mode = false;
      std::string unimod_file;
      bool reuse_existing_modifications_db = false;
    };

    struct OPENMS_DLLAPI DecoyGeneratorParameters
    {
      std::string method = "shuffle";
      std::string decoy_tag = "DECOY_";
      double min_decoy_fraction = 0.8;
      double aim_decoy_fraction = 1.0;
      Int shuffle_max_attempts = 30;
      double shuffle_sequence_identity_threshold = 0.5;
      double shift_precursor_mz_shift = 0.0;
      double shift_product_mz_shift = 20.0;
      double product_mz_threshold = 0.025;
      std::vector<std::string> allowed_fragment_types{"b", "y"};
      std::vector<size_t> allowed_fragment_charges{1, 2, 3, 4};
      bool enable_detection_specific_losses = false;
      bool enable_detection_unspecific_losses = false;
      bool switch_kr = true;
      bool separate = false;
    };

    OpenSwathLibraryPreparation();

    void setLogType(ProgressLogger::LogType log_type);

    ProgressLogger::LogType getLogType() const;

    /// Ensure the requested Unimod source is loaded before modified sequences are parsed.
    void ensureUnimodLoaded(const AssayGeneratorParameters& parameters) const;

    /// Normalize a prepared library into internal PQP format.
    LibraryStats normalizeLibraryToPQP(const std::string& input_file,
                                       FileTypes::Type input_type,
                                       const std::string& output_pqp,
                                       const Param& reader_parameters = Param()) const;

    /// Run assay generation / preparation and write the requested output format.
    LibraryStats prepareAssays(const std::string& input_file,
                               FileTypes::Type input_type,
                               const std::string& output_file,
                               FileTypes::Type output_type,
                               const AssayGeneratorParameters& parameters,
                               const Param& reader_parameters = Param()) const;

    /// Run decoy generation and write the requested output format.
    LibraryStats generateDecoys(const std::string& input_file,
                                FileTypes::Type input_type,
                                const std::string& output_file,
                                FileTypes::Type output_type,
                                const DecoyGeneratorParameters& parameters,
                                const Param& reader_parameters = Param()) const;

    /**
     * @brief Prepare an empirical library by running assay generation and decoy generation into a final PQP.
     *
     * @param input_file Input transition library file.
     * @param input_type File type of the input transition library.
     * @param output_pqp Output path for the prepared target/decoy PQP library.
     * @param assay_parameters Assay-generation parameters applied before decoy generation.
     * @param decoy_parameters Decoy-generation parameters used for the prepared assays.
     * @param reader_parameters Optional parameters passed to the transition-library reader.
     * @param scratch_directory Optional parent directory for temporary assay-preparation files.
     *        A unique per-invocation child directory is always created below this parent so concurrent
     *        callers cannot share or overwrite the same intermediate. The child directory is removed
     *        automatically on both success and failure.
     */
    LibraryStats prepareEmpiricalLibraryToPQP(const std::string& input_file,
                                              FileTypes::Type input_type,
                                              const std::string& output_pqp,
                                              const AssayGeneratorParameters& assay_parameters,
                                              const DecoyGeneratorParameters& decoy_parameters,
                                              const Param& reader_parameters = Param(),
                                              const std::string& scratch_directory = "") const;

  private:
    ProgressLogger::LogType log_type_ = ProgressLogger::CMD;
  };
} // namespace OpenMS
