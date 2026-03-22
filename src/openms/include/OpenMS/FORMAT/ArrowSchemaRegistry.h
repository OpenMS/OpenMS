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

} // namespace OpenMS

#endif // WITH_PARQUET
