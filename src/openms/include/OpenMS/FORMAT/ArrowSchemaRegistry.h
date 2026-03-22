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

} // namespace OpenMS

#endif // WITH_PARQUET
