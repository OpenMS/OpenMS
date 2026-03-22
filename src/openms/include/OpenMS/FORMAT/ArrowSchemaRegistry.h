#pragma once

#ifdef WITH_PARQUET

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

  namespace ArrowSchemaValidation
  {
    enum class Mode
    {
      Strict,
      Subset
    };

    struct OPENMS_DLLAPI ValidationResult
    {
      bool valid = true;
      std::vector<std::string> errors;
      std::string toString() const;
    };

    OPENMS_DLLAPI ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode = Mode::Strict);
  }

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
    static constexpr const char* RANK = "rank";
    static constexpr const char* PEPTIDE_IDENTIFICATION_INDEX = "peptide_identification_index";
    static constexpr const char* PSM_METAVALUES = "psm_metavalues";
    static constexpr const char* SPECTRUM_METAVALUES = "spectrum_metavalues";
    static constexpr const char* RUN_IDENTIFIER = "run_identifier";

    static std::shared_ptr<arrow::DataType> modificationsType();
    static std::shared_ptr<arrow::DataType> additionalScoresType();
    static std::shared_ptr<arrow::DataType> metavaluesType();
    static std::shared_ptr<arrow::Schema> schema();
  };

} // namespace OpenMS

#endif // WITH_PARQUET
