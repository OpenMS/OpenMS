// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <memory>
#include <string>
#include <vector>

// Forward declarations
namespace arrow
{
  class Schema;
  class DataType;
  class Table;
}

namespace OpenMS
{

  /// @brief Utilities for validating Arrow table schemas against expected definitions
  namespace ArrowSchemaValidation
  {
    /// @brief Validation strictness: Strict requires exact match, Subset allows missing and extra columns
    enum class Mode
    {
      Strict,
      Subset
    };

    /// @brief Result of schema validation containing validity flag and error messages
    struct OPENMS_DLLAPI ValidationResult
    {
      bool valid = true;
      std::vector<std::string> errors;
      std::string toString() const;
    };

    /// @brief Validate an Arrow table's schema against an expected schema
    OPENMS_DLLAPI ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode = Mode::Strict);
  }

  /// @brief Schema for protein identification results table
  struct OPENMS_DLLAPI ProteinSchema
  {
    static constexpr const char* ACCESSION = "accession";
    static constexpr const char* SCORE = "score";
    static constexpr const char* RANK = "rank";
    static constexpr const char* COVERAGE = "coverage";
    static constexpr const char* SEQUENCE = "sequence";
    static constexpr const char* DESCRIPTION = "description";
    static constexpr const char* IS_DECOY = "is_decoy";
    static constexpr const char* RUN_IDENTIFIER = "run_identifier";
    static constexpr const char* MODIFICATIONS = "modifications";
    static constexpr const char* METAVALUES = "metavalues";

    static std::shared_ptr<arrow::DataType> modificationsType();
    static std::shared_ptr<arrow::DataType> metavaluesType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for protein group (indistinguishable group) results table
  struct OPENMS_DLLAPI ProteinGroupSchema
  {
    static constexpr const char* GROUP_TYPE = "group_type";
    static constexpr const char* PROBABILITY = "probability";
    static constexpr const char* ACCESSIONS = "accessions";
    static constexpr const char* RUN_IDENTIFIER = "run_identifier";
    static constexpr const char* GROUP_INDEX = "group_index";
    static constexpr const char* FLOAT_DATA = "float_data";
    static constexpr const char* STRING_DATA = "string_data";
    static constexpr const char* INTEGER_DATA = "integer_data";

    static std::shared_ptr<arrow::DataType> floatDataType();
    static std::shared_ptr<arrow::DataType> stringDataType();
    static std::shared_ptr<arrow::DataType> integerDataType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for search engine parameters and settings table
  struct OPENMS_DLLAPI SearchParamsSchema
  {
    static constexpr const char* RUN_IDENTIFIER = "run_identifier";
    static constexpr const char* SEARCH_ENGINE = "search_engine";
    static constexpr const char* SEARCH_ENGINE_VERSION = "search_engine_version";
    static constexpr const char* INFERENCE_ENGINE = "inference_engine";
    static constexpr const char* INFERENCE_ENGINE_VERSION = "inference_engine_version";
    static constexpr const char* DATE = "date";
    static constexpr const char* SCORE_TYPE = "score_type";
    static constexpr const char* HIGHER_SCORE_BETTER = "higher_score_better";
    static constexpr const char* SIGNIFICANCE_THRESHOLD = "significance_threshold";
    static constexpr const char* DB = "db";
    static constexpr const char* DB_VERSION = "db_version";
    static constexpr const char* TAXONOMY = "taxonomy";
    static constexpr const char* CHARGES = "charges";
    static constexpr const char* MASS_TYPE = "mass_type";
    static constexpr const char* PRECURSOR_MASS_TOLERANCE = "precursor_mass_tolerance";
    static constexpr const char* PRECURSOR_MASS_TOLERANCE_PPM = "precursor_mass_tolerance_ppm";
    static constexpr const char* FRAGMENT_MASS_TOLERANCE = "fragment_mass_tolerance";
    static constexpr const char* FRAGMENT_MASS_TOLERANCE_PPM = "fragment_mass_tolerance_ppm";
    static constexpr const char* DIGESTION_ENZYME = "digestion_enzyme";
    static constexpr const char* ENZYME_TERM_SPECIFICITY = "enzyme_term_specificity";
    static constexpr const char* MISSED_CLEAVAGES = "missed_cleavages";
    static constexpr const char* FIXED_MODIFICATIONS = "fixed_modifications";
    static constexpr const char* VARIABLE_MODIFICATIONS = "variable_modifications";
    static constexpr const char* PRIMARY_MS_RUN_PATHS = "primary_ms_run_paths";
    static constexpr const char* METAVALUES = "metavalues";
    static constexpr const char* SP_METAVALUES = "sp_metavalues";

    static std::shared_ptr<arrow::DataType> metavaluesType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for LC-MS feature table (FeatureMap features)
  struct OPENMS_DLLAPI FeatureSchema
  {
    static constexpr const char* UNIQUE_ID = "unique_id";
    static constexpr const char* PARENT_FEATURE_ID = "parent_feature_id";
    static constexpr const char* DEPTH = "depth";
    static constexpr const char* RT = "rt";
    static constexpr const char* MZ = "mz";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* QUALITY = "quality";
    static constexpr const char* QUALITY_RT = "quality_rt";
    static constexpr const char* QUALITY_MZ = "quality_mz";
    static constexpr const char* WIDTH = "width";
    static constexpr const char* RT_BB_MIN = "rt_bb_min";
    static constexpr const char* RT_BB_MAX = "rt_bb_max";
    static constexpr const char* MZ_BB_MIN = "mz_bb_min";
    static constexpr const char* MZ_BB_MAX = "mz_bb_max";
    static constexpr const char* CONVEX_HULLS = "convex_hulls";
    static constexpr const char* METAVALUES = "metavalues";

    static std::shared_ptr<arrow::DataType> convexHullType();
    static std::shared_ptr<arrow::DataType> metavaluesType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for consensus feature table (ConsensusMap features)
  struct OPENMS_DLLAPI ConsensusFeatureSchema
  {
    static constexpr const char* UNIQUE_ID = "unique_id";
    static constexpr const char* RT = "rt";
    static constexpr const char* MZ = "mz";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* QUALITY = "quality";
    static constexpr const char* WIDTH = "width";
    static constexpr const char* HANDLES = "handles";
    static constexpr const char* METAVALUES = "metavalues";

    static std::shared_ptr<arrow::DataType> handlesType();
    static std::shared_ptr<arrow::DataType> metavaluesType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for peptide-spectrum match (PSM) results table
  struct OPENMS_DLLAPI PSMSchema
  {
    static constexpr const char* SEQUENCE = "sequence";
    static constexpr const char* PEPTIDOFORM = "peptidoform";
    static constexpr const char* MODIFICATIONS = "modifications";
    static constexpr const char* PRECURSOR_CHARGE = "precursor_charge";
    static constexpr const char* POSTERIOR_ERROR_PROBABILITY = "posterior_error_probability";
    static constexpr const char* IS_DECOY = "is_decoy";
    static constexpr const char* CALCULATED_MZ = "calculated_mz";
    static constexpr const char* OBSERVED_MZ = "observed_mz";
    static constexpr const char* ADDITIONAL_SCORES = "additional_scores";
    static constexpr const char* PROTEIN_ACCESSIONS = "protein_accessions";
    static constexpr const char* PREDICTED_RT = "predicted_rt";
    static constexpr const char* REFERENCE_FILE_NAME = "reference_file_name";
    static constexpr const char* CV_PARAMS = "cv_params";
    static constexpr const char* SCAN = "scan";
    static constexpr const char* RT = "rt";
    static constexpr const char* ION_MOBILITY = "ion_mobility";
    static constexpr const char* SPECTRUM_REFERENCE = "spectrum_reference";
    static constexpr const char* SCORE = "score";
    static constexpr const char* SCORE_TYPE = "score_type";
    static constexpr const char* HIGHER_SCORE_BETTER = "higher_score_better";
    /// 0-based positional column: hit position within parent identification.
    /// Distinct from PeptideHit::getRank(), which is round-tripped via the
    /// "rank" UserParam carried in PSM_METAVALUES.
    static constexpr const char* HIT_INDEX = "hit_index";
    static constexpr const char* PEPTIDE_IDENTIFICATION_INDEX = "peptide_identification_index";
    static constexpr const char* PSM_METAVALUES = "psm_metavalues";
    static constexpr const char* SPECTRUM_METAVALUES = "spectrum_metavalues";
    static constexpr const char* RUN_IDENTIFIER = "run_identifier";
    static constexpr const char* MZ_ARRAY = "mz_array";
    static constexpr const char* INTENSITY_ARRAY = "intensity_array";
    static constexpr const char* CHARGE_ARRAY = "charge_array";
    static constexpr const char* ION_TYPE_ARRAY = "ion_type_array";

    static std::shared_ptr<arrow::DataType> modificationsType();
    static std::shared_ptr<arrow::DataType> additionalScoresType();
    static std::shared_ptr<arrow::DataType> metavaluesType();
    /**
      @brief Arrow type for protein_accessions: list<struct{accession, aa_before, aa_after, start, end}>
    */
    static std::shared_ptr<arrow::DataType> proteinAccessionsType();
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for QPX PSM export (quantms Parquet eXchange format, PSM table)
  ///
  /// Defines column names, nested Arrow types, and the complete schema for
  /// peptide-spectrum match data in the QPX format. Used by QPXFile for export.
  struct OPENMS_DLLAPI QPXPSMSchema
  {
    static constexpr const char* SEQUENCE = "sequence";
    static constexpr const char* PEPTIDOFORM = "peptidoform";
    static constexpr const char* MODIFICATIONS = "modifications";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* POSTERIOR_ERROR_PROBABILITY = "posterior_error_probability";
    static constexpr const char* IS_DECOY = "is_decoy";
    static constexpr const char* CALCULATED_MZ = "calculated_mz";
    static constexpr const char* OBSERVED_MZ = "observed_mz";
    static constexpr const char* MASS_ERROR_PPM = "mass_error_ppm";
    static constexpr const char* ADDITIONAL_SCORES = "additional_scores";
    static constexpr const char* PREDICTED_RT = "predicted_rt";
    static constexpr const char* RUN_FILE_NAME = "run_file_name";
    static constexpr const char* CV_PARAMS = "cv_params";
    static constexpr const char* SCAN = "scan";
    static constexpr const char* RT = "rt";
    static constexpr const char* ION_MOBILITY = "ion_mobility";
    static constexpr const char* MISSED_CLEAVAGES = "missed_cleavages";
    static constexpr const char* PROTEIN_ACCESSIONS = "protein_accessions";
    static constexpr const char* CROSS_LINKS = "cross_links";
    static constexpr const char* MZ_ARRAY = "mz_array";
    static constexpr const char* INTENSITY_ARRAY = "intensity_array";
    static constexpr const char* CHARGE_ARRAY = "charge_array";
    static constexpr const char* ION_TYPE_ARRAY = "ion_type_array";
    static constexpr const char* ION_MOBILITY_ARRAY = "ion_mobility_array";

    /// @brief Arrow type for modifications: list<struct{name, accession, positions: list<struct{position, amino_acid, scores}>}>
    static std::shared_ptr<arrow::DataType> modificationsType();
    /// @brief Arrow type for additional scores: list<struct{score_name, score_value, higher_better}>
    static std::shared_ptr<arrow::DataType> additionalScoresType();
    /// @brief Arrow type for CV params: list<struct{cv_name, cv_value}>
    static std::shared_ptr<arrow::DataType> cvParamsType();
    /// @brief Arrow type for cross-links: list<struct{xl_type, partner_sequence, ...}>
    static std::shared_ptr<arrow::DataType> crossLinksType();
    /// @brief Complete Arrow schema for QPX PSM table (24 fields)
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for QPX feature view (quantms Parquet eXchange format)
  ///
  /// Defines column names, nested Arrow types, and the complete schema for
  /// consensus feature data in the QPX format. Used by ConsensusMapArrowExport.
  struct OPENMS_DLLAPI QPXFeatureSchema
  {
    static constexpr const char* SEQUENCE = "sequence";
    static constexpr const char* PEPTIDOFORM = "peptidoform";
    static constexpr const char* MODIFICATIONS = "modifications";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* POSTERIOR_ERROR_PROBABILITY = "posterior_error_probability";
    static constexpr const char* IS_DECOY = "is_decoy";
    static constexpr const char* CALCULATED_MZ = "calculated_mz";
    static constexpr const char* OBSERVED_MZ = "observed_mz";
    static constexpr const char* MASS_ERROR_PPM = "mass_error_ppm";
    static constexpr const char* ADDITIONAL_SCORES = "additional_scores";
    static constexpr const char* PREDICTED_RT = "predicted_rt";
    static constexpr const char* RUN_FILE_NAME = "run_file_name";
    static constexpr const char* CV_PARAMS = "cv_params";
    static constexpr const char* SCAN = "scan";
    static constexpr const char* RT = "rt";
    static constexpr const char* ION_MOBILITY = "ion_mobility";
    static constexpr const char* MISSED_CLEAVAGES = "missed_cleavages";
    static constexpr const char* INTENSITIES = "intensities";
    static constexpr const char* ADDITIONAL_INTENSITIES = "additional_intensities";
    static constexpr const char* PG_ACCESSIONS = "pg_accessions";
    static constexpr const char* ANCHOR_PROTEIN = "anchor_protein";
    static constexpr const char* UNIQUE = "unique";
    static constexpr const char* PG_GLOBAL_QVALUE = "pg_global_qvalue";
    static constexpr const char* PG_POSITIONS = "pg_positions";
    static constexpr const char* ION_MOBILITY_START = "ion_mobility_start";
    static constexpr const char* ION_MOBILITY_STOP = "ion_mobility_stop";
    static constexpr const char* GG_ACCESSIONS = "gg_accessions";
    static constexpr const char* GG_NAMES = "gg_names";
    static constexpr const char* ID_RUN_FILE_NAME = "id_run_file_name";
    static constexpr const char* RT_START = "rt_start";
    static constexpr const char* RT_STOP = "rt_stop";

    /// @brief Arrow type for modifications (delegates to QPXPSMSchema::modificationsType)
    static std::shared_ptr<arrow::DataType> modificationsType();
    /// @brief Arrow type for additional scores (delegates to QPXPSMSchema::additionalScoresType)
    static std::shared_ptr<arrow::DataType> additionalScoresType();
    /// @brief Arrow type for CV params (delegates to QPXPSMSchema::cvParamsType)
    static std::shared_ptr<arrow::DataType> cvParamsType();
    /// @brief Arrow type for intensities: list<struct{label, intensity}>
    static std::shared_ptr<arrow::DataType> intensitiesType();
    /// @brief Arrow type for additional intensities: list<struct{label, intensities: list<struct{...}>}>
    static std::shared_ptr<arrow::DataType> additionalIntensitiesType();
    /// @brief Arrow type for protein group accessions: list<struct{accession, start, end, pre, post}>
    static std::shared_ptr<arrow::DataType> pgAccessionsType();
    /// @brief Arrow type for protein group positions: list<struct{protein_accession, start, end}>
    static std::shared_ptr<arrow::DataType> pgPositionsType();
    /// @brief Complete Arrow schema for QPX feature table (31 fields)
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for QPX protein group export (quantms Parquet eXchange format, pg table)
  ///
  /// Defines column names, nested Arrow types, and the complete schema for
  /// protein group data in the QPX format. Supports both quantified (ConsensusMap)
  /// and identification-only (search engine) output — quantification columns are nullable.
  struct OPENMS_DLLAPI QPXPgSchema
  {
    static constexpr const char* PG_ACCESSIONS = "pg_accessions";
    static constexpr const char* PG_NAMES = "pg_names";
    static constexpr const char* GG_ACCESSIONS = "gg_accessions";
    static constexpr const char* GG_NAMES = "gg_names";
    static constexpr const char* GG_QVALUE = "gg_qvalue";
    static constexpr const char* ANCHOR_PROTEIN = "anchor_protein";
    static constexpr const char* RUN_FILE_NAME = "run_file_name";
    static constexpr const char* GLOBAL_QVALUE = "global_qvalue";
    static constexpr const char* PG_QVALUE = "pg_qvalue";
    static constexpr const char* INTENSITIES = "intensities";
    static constexpr const char* ADDITIONAL_INTENSITIES = "additional_intensities";
    static constexpr const char* IS_DECOY = "is_decoy";
    static constexpr const char* CONTAMINANT = "contaminant";
    static constexpr const char* PEPTIDES = "peptides";
    static constexpr const char* PEPTIDE_COUNTS = "peptide_counts";
    static constexpr const char* FEATURE_COUNTS = "feature_counts";
    static constexpr const char* SEQUENCE_COVERAGE = "sequence_coverage";
    static constexpr const char* MOLECULAR_WEIGHT = "molecular_weight";
    static constexpr const char* ADDITIONAL_SCORES = "additional_scores";
    static constexpr const char* CV_PARAMS = "cv_params";

    /// @brief Arrow type for intensities: list<struct{label, intensity}> (nullable for search-engine output)
    static std::shared_ptr<arrow::DataType> intensitiesType();
    /// @brief Arrow type for additional intensities: list<struct{label, intensities: list<struct{...}>}>
    static std::shared_ptr<arrow::DataType> additionalIntensitiesType();
    /// @brief Arrow type for peptides: list<struct{protein_name, peptide_count}>
    static std::shared_ptr<arrow::DataType> peptidesType();
    /// @brief Arrow type for peptide_counts: struct{unique_sequences, total_sequences}
    static std::shared_ptr<arrow::DataType> peptideCountsType();
    /// @brief Arrow type for feature_counts: struct{unique_features, total_features}
    static std::shared_ptr<arrow::DataType> featureCountsType();
    /// @brief Arrow type for additional scores (delegates to QPXPSMSchema::additionalScoresType)
    static std::shared_ptr<arrow::DataType> additionalScoresType();
    /// @brief Arrow type for CV params (delegates to QPXPSMSchema::cvParamsType)
    static std::shared_ptr<arrow::DataType> cvParamsType();
    /// @brief Complete Arrow schema for QPX pg table (20 fields)
    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for spectra in long (one row per peak) format
  struct OPENMS_DLLAPI SpectraLongSchema
  {
    static constexpr const char* MZ = "mz";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* RT = "rt";
    static constexpr const char* ION_MOBILITY = "ion_mobility";
    static constexpr const char* SPECTRUM_INDEX = "spectrum_index";
    static constexpr const char* MS_LEVEL = "ms_level";
    static constexpr const char* NATIVE_ID = "native_id";
    static constexpr const char* PRECURSOR_MZ = "precursor_mz";
    static constexpr const char* PRECURSOR_CHARGE = "precursor_charge";
    static constexpr const char* PRECURSOR_INTENSITY = "precursor_intensity";
    static constexpr const char* ISOLATION_LOWER = "isolation_lower";
    static constexpr const char* ISOLATION_UPPER = "isolation_upper";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for spectra in semi-wide (one row per spectrum, list columns for peaks) format
  struct OPENMS_DLLAPI SpectraSemiWideSchema
  {
    static constexpr const char* SPECTRUM_INDEX = "spectrum_index";
    static constexpr const char* RT = "rt";
    static constexpr const char* MS_LEVEL = "ms_level";
    static constexpr const char* NATIVE_ID = "native_id";
    static constexpr const char* MZ = "mz";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* ION_MOBILITY = "ion_mobility";
    static constexpr const char* PRECURSOR_MZ = "precursor_mz";
    static constexpr const char* PRECURSOR_CHARGE = "precursor_charge";
    static constexpr const char* PRECURSOR_INTENSITY = "precursor_intensity";
    static constexpr const char* ISOLATION_LOWER = "isolation_lower";
    static constexpr const char* ISOLATION_UPPER = "isolation_upper";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for chromatograms in long (one row per data point) format
  struct OPENMS_DLLAPI ChromatogramSchema
  {
    static constexpr const char* RT = "rt";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* CHROMATOGRAM_INDEX = "chromatogram_index";
    static constexpr const char* NATIVE_ID = "native_id";
    static constexpr const char* PRECURSOR_MZ = "precursor_mz";
    static constexpr const char* PRODUCT_MZ = "product_mz";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for chromatograms in semi-wide (one row per chromatogram, list columns) format
  struct OPENMS_DLLAPI ChromatogramSemiWideSchema
  {
    static constexpr const char* CHROMATOGRAM_INDEX = "chromatogram_index";
    static constexpr const char* NATIVE_ID = "native_id";
    static constexpr const char* RT = "rt";
    static constexpr const char* INTENSITY = "intensity";
    static constexpr const char* PRECURSOR_MZ = "precursor_mz";
    static constexpr const char* PRODUCT_MZ = "product_mz";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH precursor (peptide query) table
  struct OPENMS_DLLAPI OSWPrecursorSchema
  {
    static constexpr const char* PRECURSOR_ID = "precursor_id";
    static constexpr const char* PRECURSOR_MZ = "precursor_mz";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* LIBRARY_RT = "library_rt";
    static constexpr const char* LIBRARY_DRIFT_TIME = "library_drift_time";
    static constexpr const char* DECOY = "decoy";
    static constexpr const char* TRAML_ID = "traml_id";
    static constexpr const char* MODIFIED_SEQUENCE = "modified_sequence";
    static constexpr const char* UNMODIFIED_SEQUENCE = "unmodified_sequence";
    static constexpr const char* PROTEIN_ACCESSIONS = "protein_accessions";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH transition (fragment ion) table
  struct OPENMS_DLLAPI OSWTransitionSchema
  {
    static constexpr const char* TRANSITION_ID = "transition_id";
    static constexpr const char* PRECURSOR_ID = "precursor_id";
    static constexpr const char* TRAML_ID = "traml_id";
    static constexpr const char* PRODUCT_MZ = "product_mz";
    static constexpr const char* CHARGE = "charge";
    static constexpr const char* TYPE = "type";
    static constexpr const char* ANNOTATION = "annotation";
    static constexpr const char* ORDINAL = "ordinal";
    static constexpr const char* DETECTING = "detecting";
    static constexpr const char* IDENTIFYING = "identifying";
    static constexpr const char* QUANTIFYING = "quantifying";
    static constexpr const char* LIBRARY_INTENSITY = "library_intensity";
    static constexpr const char* DECOY = "decoy";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH feature-level precursor intensity table
  struct OPENMS_DLLAPI OSWFeaturePrecursorSchema
  {
    static constexpr const char* FEATURE_ID = "feature_id";
    static constexpr const char* RUN_ID = "run_id";
    static constexpr const char* PRECURSOR_ISOTOPE = "precursor_isotope";
    static constexpr const char* PRECURSOR_AREA_INTENSITY = "precursor_area_intensity";
    static constexpr const char* PRECURSOR_APEX_INTENSITY = "precursor_apex_intensity";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH run metadata table
  struct OPENMS_DLLAPI OSWRunSchema
  {
    static constexpr const char* RUN_ID = "run_id";
    static constexpr const char* FILENAME = "filename";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH feature scoring results table
  struct OPENMS_DLLAPI OSWFeatureSchema
  {
    static constexpr const char* FEATURE_ID = "feature_id";
    static constexpr const char* RUN_ID = "run_id";
    static constexpr const char* PRECURSOR_ID = "precursor_id";
    static constexpr const char* EXP_RT = "exp_rt";
    static constexpr const char* EXP_IM = "exp_im";
    static constexpr const char* NORM_RT = "norm_rt";
    static constexpr const char* DELTA_RT = "delta_rt";
    static constexpr const char* LEFT_WIDTH = "left_width";
    static constexpr const char* RIGHT_WIDTH = "right_width";
    static constexpr const char* EXP_IM_LEFTWIDTH = "exp_im_leftwidth";
    static constexpr const char* EXP_IM_RIGHTWIDTH = "exp_im_rightwidth";
    static constexpr const char* MS1_AREA_INTENSITY = "ms1_area_intensity";
    static constexpr const char* MS1_APEX_INTENSITY = "ms1_apex_intensity";
    static constexpr const char* MS1_EXP_IM = "ms1_exp_im";
    static constexpr const char* MS1_DELTA_IM = "ms1_delta_im";
    static constexpr const char* VAR_MS1_MASSDEV_SCORE = "var_ms1_massdev_score";
    static constexpr const char* VAR_MS1_IM_MS1_DELTA_SCORE = "var_ms1_im_ms1_delta_score";
    static constexpr const char* VAR_MS1_MI_SCORE = "var_ms1_mi_score";
    static constexpr const char* VAR_MS1_MI_CONTRAST_SCORE = "var_ms1_mi_contrast_score";
    static constexpr const char* VAR_MS1_MI_COMBINED_SCORE = "var_ms1_mi_combined_score";
    static constexpr const char* VAR_MS1_ISOTOPE_CORRELATION_SCORE = "var_ms1_isotope_correlation_score";
    static constexpr const char* VAR_MS1_ISOTOPE_OVERLAP_SCORE = "var_ms1_isotope_overlap_score";
    static constexpr const char* VAR_MS1_XCORR_COELUTION = "var_ms1_xcorr_coelution";
    static constexpr const char* VAR_MS1_XCORR_COELUTION_CONTRAST = "var_ms1_xcorr_coelution_contrast";
    static constexpr const char* VAR_MS1_XCORR_COELUTION_COMBINED = "var_ms1_xcorr_coelution_combined";
    static constexpr const char* VAR_MS1_XCORR_SHAPE = "var_ms1_xcorr_shape";
    static constexpr const char* VAR_MS1_XCORR_SHAPE_CONTRAST = "var_ms1_xcorr_shape_contrast";
    static constexpr const char* VAR_MS1_XCORR_SHAPE_COMBINED = "var_ms1_xcorr_shape_combined";
    static constexpr const char* MS2_AREA_INTENSITY = "ms2_area_intensity";
    static constexpr const char* MS2_TOTAL_AREA_INTENSITY = "ms2_total_area_intensity";
    static constexpr const char* MS2_APEX_INTENSITY = "ms2_apex_intensity";
    static constexpr const char* MS2_EXP_IM = "ms2_exp_im";
    static constexpr const char* MS2_EXP_IM_LEFTWIDTH = "ms2_exp_im_leftwidth";
    static constexpr const char* MS2_EXP_IM_RIGHTWIDTH = "ms2_exp_im_rightwidth";
    static constexpr const char* MS2_DELTA_IM = "ms2_delta_im";
    static constexpr const char* MS2_TOTAL_MI = "ms2_total_mi";
    static constexpr const char* VAR_MS2_BSERIES_SCORE = "var_ms2_bseries_score";
    static constexpr const char* VAR_MS2_DOTPROD_SCORE = "var_ms2_dotprod_score";
    static constexpr const char* VAR_MS2_INTENSITY_SCORE = "var_ms2_intensity_score";
    static constexpr const char* VAR_MS2_ISOTOPE_CORRELATION_SCORE = "var_ms2_isotope_correlation_score";
    static constexpr const char* VAR_MS2_ISOTOPE_OVERLAP_SCORE = "var_ms2_isotope_overlap_score";
    static constexpr const char* VAR_MS2_LIBRARY_CORR = "var_ms2_library_corr";
    static constexpr const char* VAR_MS2_LIBRARY_DOTPROD = "var_ms2_library_dotprod";
    static constexpr const char* VAR_MS2_LIBRARY_MANHATTAN = "var_ms2_library_manhattan";
    static constexpr const char* VAR_MS2_LIBRARY_RMSD = "var_ms2_library_rmsd";
    static constexpr const char* VAR_MS2_LIBRARY_ROOTMEANSQUARE = "var_ms2_library_rootmeansquare";
    static constexpr const char* VAR_MS2_LIBRARY_SANGLE = "var_ms2_library_sangle";
    static constexpr const char* VAR_MS2_LOG_SN_SCORE = "var_ms2_log_sn_score";
    static constexpr const char* VAR_MS2_MANHATTAN_SCORE = "var_ms2_manhattan_score";
    static constexpr const char* VAR_MS2_MASSDEV_SCORE = "var_ms2_massdev_score";
    static constexpr const char* VAR_MS2_MASSDEV_SCORE_WEIGHTED = "var_ms2_massdev_score_weighted";
    static constexpr const char* VAR_MS2_MI_SCORE = "var_ms2_mi_score";
    static constexpr const char* VAR_MS2_MI_WEIGHTED_SCORE = "var_ms2_mi_weighted_score";
    static constexpr const char* VAR_MS2_MI_RATIO_SCORE = "var_ms2_mi_ratio_score";
    static constexpr const char* VAR_MS2_NORM_RT_SCORE = "var_ms2_norm_rt_score";
    static constexpr const char* VAR_MS2_XCORR_COELUTION = "var_ms2_xcorr_coelution";
    static constexpr const char* VAR_MS2_XCORR_COELUTION_WEIGHTED = "var_ms2_xcorr_coelution_weighted";
    static constexpr const char* VAR_MS2_XCORR_SHAPE = "var_ms2_xcorr_shape";
    static constexpr const char* VAR_MS2_XCORR_SHAPE_WEIGHTED = "var_ms2_xcorr_shape_weighted";
    static constexpr const char* VAR_MS2_YSERIES_SCORE = "var_ms2_yseries_score";
    static constexpr const char* VAR_MS2_ELUTION_MODEL_FIT_SCORE = "var_ms2_elution_model_fit_score";
    static constexpr const char* VAR_MS2_IM_XCORR_SHAPE = "var_ms2_im_xcorr_shape";
    static constexpr const char* VAR_MS2_IM_XCORR_COELUTION = "var_ms2_im_xcorr_coelution";
    static constexpr const char* VAR_MS2_IM_DELTA_SCORE = "var_ms2_im_delta_score";
    static constexpr const char* VAR_MS2_IM_LOG_INTENSITY = "var_ms2_im_log_intensity";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for OpenSWATH per-transition feature scoring results table
  struct OPENMS_DLLAPI OSWFeatureTransitionSchema
  {
    static constexpr const char* FEATURE_ID = "feature_id";
    static constexpr const char* RUN_ID = "run_id";
    static constexpr const char* TRANSITION_ID = "transition_id";
    static constexpr const char* AREA_INTENSITY = "area_intensity";
    static constexpr const char* TOTAL_AREA_INTENSITY = "total_area_intensity";
    static constexpr const char* APEX_INTENSITY = "apex_intensity";
    static constexpr const char* APEX_RT = "apex_rt";
    static constexpr const char* RT_FWHM = "rt_fwhm";
    static constexpr const char* MASSERROR_PPM = "masserror_ppm";
    static constexpr const char* TOTAL_MI = "total_mi";
    static constexpr const char* VAR_INTENSITY_SCORE = "var_intensity_score";
    static constexpr const char* VAR_INTENSITY_RATIO_SCORE = "var_intensity_ratio_score";
    static constexpr const char* VAR_LOG_INTENSITY = "var_log_intensity";
    static constexpr const char* VAR_XCORR_COELUTION = "var_xcorr_coelution";
    static constexpr const char* VAR_XCORR_SHAPE = "var_xcorr_shape";
    static constexpr const char* VAR_LOG_SN_SCORE = "var_log_sn_score";
    static constexpr const char* VAR_MASSDEV_SCORE = "var_massdev_score";
    static constexpr const char* VAR_MI_SCORE = "var_mi_score";
    static constexpr const char* VAR_MI_RATIO_SCORE = "var_mi_ratio_score";
    static constexpr const char* VAR_ISOTOPE_CORRELATION_SCORE = "var_isotope_correlation_score";
    static constexpr const char* VAR_ISOTOPE_OVERLAP_SCORE = "var_isotope_overlap_score";
    static constexpr const char* EXP_IM = "exp_im";
    static constexpr const char* EXP_IM_LEFTWIDTH = "exp_im_leftwidth";
    static constexpr const char* EXP_IM_RIGHTWIDTH = "exp_im_rightwidth";
    static constexpr const char* DELTA_IM = "delta_im";
    static constexpr const char* VAR_IM_DELTA_SCORE = "var_im_delta_score";
    static constexpr const char* VAR_IM_LOG_INTENSITY = "var_im_log_intensity";
    static constexpr const char* VAR_IM_XCORR_COELUTION_CONTRAST = "var_im_xcorr_coelution_contrast";
    static constexpr const char* VAR_IM_XCORR_SHAPE_CONTRAST = "var_im_xcorr_shape_contrast";
    static constexpr const char* VAR_IM_XCORR_COELUTION_COMBINED = "var_im_xcorr_coelution_combined";
    static constexpr const char* VAR_IM_XCORR_SHAPE_COMBINED = "var_im_xcorr_shape_combined";
    static constexpr const char* START_POSITION_AT_5 = "start_position_at_5";
    static constexpr const char* END_POSITION_AT_5 = "end_position_at_5";
    static constexpr const char* START_POSITION_AT_10 = "start_position_at_10";
    static constexpr const char* END_POSITION_AT_10 = "end_position_at_10";
    static constexpr const char* START_POSITION_AT_50 = "start_position_at_50";
    static constexpr const char* END_POSITION_AT_50 = "end_position_at_50";
    static constexpr const char* TOTAL_WIDTH = "total_width";
    static constexpr const char* TAILING_FACTOR = "tailing_factor";
    static constexpr const char* ASYMMETRY_FACTOR = "asymmetry_factor";
    static constexpr const char* SLOPE_OF_BASELINE = "slope_of_baseline";
    static constexpr const char* BASELINE_DELTA_2_HEIGHT = "baseline_delta_2_height";
    static constexpr const char* POINTS_ACROSS_BASELINE = "points_across_baseline";
    static constexpr const char* POINTS_ACROSS_HALF_HEIGHT = "points_across_half_height";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for extracted ion chromatogram (XIC) data table
  struct OPENMS_DLLAPI XICSchema
  {
    static constexpr const char* RUN_ID = "RUN_ID";
    static constexpr const char* SOURCE_FILE = "SOURCE_FILE";
    static constexpr const char* MS_LEVEL = "MS_LEVEL";
    static constexpr const char* PRECURSOR_ID = "PRECURSOR_ID";
    static constexpr const char* TRANSITION_ID = "TRANSITION_ID";
    static constexpr const char* MODIFIED_SEQUENCE = "MODIFIED_SEQUENCE";
    static constexpr const char* PRECURSOR_CHARGE = "PRECURSOR_CHARGE";
    static constexpr const char* PRODUCT_CHARGE = "PRODUCT_CHARGE";
    static constexpr const char* DETECTING_TRANSITION = "DETECTING_TRANSITION";
    static constexpr const char* PRECURSOR_DECOY = "PRECURSOR_DECOY";
    static constexpr const char* PRODUCT_DECOY = "PRODUCT_DECOY";
    static constexpr const char* TRANSITION_ORDINAL = "TRANSITION_ORDINAL";
    static constexpr const char* TRANSITION_TYPE = "TRANSITION_TYPE";
    static constexpr const char* ANNOTATION = "ANNOTATION";
    static constexpr const char* RT_DATA = "RT_DATA";
    static constexpr const char* INTENSITY_DATA = "INTENSITY_DATA";
    static constexpr const char* RT_COMPRESSION = "RT_COMPRESSION";
    static constexpr const char* INTENSITY_COMPRESSION = "INTENSITY_COMPRESSION";

    static std::shared_ptr<arrow::Schema> schema();
  };

  /// @brief Schema for extracted ion mobilogram (XIM) data table
  struct OPENMS_DLLAPI XIMSchema
  {
    static constexpr const char* RUN_ID = "RUN_ID";
    static constexpr const char* SOURCE_FILE = "SOURCE_FILE";
    static constexpr const char* MS_LEVEL = "MS_LEVEL";
    static constexpr const char* MOBILOGRAM_TYPE = "MOBILOGRAM_TYPE";
    static constexpr const char* PRECURSOR_ID = "PRECURSOR_ID";
    static constexpr const char* TRANSITION_ID = "TRANSITION_ID";
    static constexpr const char* FEATURE_ID = "FEATURE_ID";
    static constexpr const char* FEATURE_RT = "FEATURE_RT";
    static constexpr const char* MODIFIED_SEQUENCE = "MODIFIED_SEQUENCE";
    static constexpr const char* PRECURSOR_CHARGE = "PRECURSOR_CHARGE";
    static constexpr const char* PRODUCT_CHARGE = "PRODUCT_CHARGE";
    static constexpr const char* DETECTING_TRANSITION = "DETECTING_TRANSITION";
    static constexpr const char* PRECURSOR_DECOY = "PRECURSOR_DECOY";
    static constexpr const char* PRODUCT_DECOY = "PRODUCT_DECOY";
    static constexpr const char* TRANSITION_ORDINAL = "TRANSITION_ORDINAL";
    static constexpr const char* TRANSITION_TYPE = "TRANSITION_TYPE";
    static constexpr const char* ANNOTATION = "ANNOTATION";
    static constexpr const char* MOBILITY_DATA = "MOBILITY_DATA";
    static constexpr const char* INTENSITY_DATA = "INTENSITY_DATA";
    static constexpr const char* MOBILITY_COMPRESSION = "MOBILITY_COMPRESSION";
    static constexpr const char* INTENSITY_COMPRESSION = "INTENSITY_COMPRESSION";

    static std::shared_ptr<arrow::Schema> schema();
  };

} // namespace OpenMS
