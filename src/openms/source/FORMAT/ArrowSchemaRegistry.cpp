#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>
#include <sstream>

namespace OpenMS
{

  namespace ArrowSchemaValidation
  {
    std::string ValidationResult::toString() const
    {
      std::ostringstream oss;
      for (size_t i = 0; i < errors.size(); ++i)
      {
        if (i > 0) oss << "; ";
        oss << errors[i];
      }
      return oss.str();
    }

    ValidationResult validate(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected_schema,
      Mode mode)
    {
      ValidationResult result;
      auto actual_schema = table->schema();
      int actual_n = actual_schema->num_fields();
      int expected_n = expected_schema->num_fields();

      if (mode == Mode::Strict)
      {
        if (actual_n != expected_n)
        {
          result.valid = false;
          result.errors.push_back("Field count mismatch: got " + std::to_string(actual_n) +
            ", expected " + std::to_string(expected_n));
        }

        for (int i = 0; i < std::min(actual_n, expected_n); ++i)
        {
          auto actual_field = actual_schema->field(i);
          auto expected_field = expected_schema->field(i);

          if (actual_field->name() != expected_field->name())
          {
            result.valid = false;
            result.errors.push_back("Field name mismatch at index " + std::to_string(i) +
              ": got '" + actual_field->name() + "', expected '" + expected_field->name() + "'");
            continue;
          }
          if (!actual_field->type()->Equals(expected_field->type()))
          {
            result.valid = false;
            result.errors.push_back("Type mismatch for field '" + actual_field->name() +
              "': got " + actual_field->type()->ToString() + ", expected " + expected_field->type()->ToString());
          }
          if (actual_field->nullable() != expected_field->nullable())
          {
            result.valid = false;
            result.errors.push_back("Nullability mismatch for field '" + actual_field->name() +
              "': got " + std::string(actual_field->nullable() ? "nullable" : "non-null") +
              ", expected " + std::string(expected_field->nullable() ? "nullable" : "non-null"));
          }
        }

        // Check for extra fields beyond expected count
        for (int i = expected_n; i < actual_n; ++i)
        {
          result.valid = false;
          result.errors.push_back("Unexpected field '" + actual_schema->field(i)->name() +
            "' not in expected schema");
        }

        // Check for missing fields beyond actual count
        for (int i = actual_n; i < expected_n; ++i)
        {
          result.valid = false;
          result.errors.push_back("Missing field '" + expected_schema->field(i)->name() + "'");
        }
      }
      else // Mode::Subset
      {
        for (int i = 0; i < actual_n; ++i)
        {
          auto actual_field = actual_schema->field(i);
          auto expected_field = expected_schema->GetFieldByName(actual_field->name());

          if (!expected_field)
          {
            result.valid = false;
            result.errors.push_back("Unexpected field '" + actual_field->name() +
              "' not in expected schema");
            continue;
          }
          if (!actual_field->type()->Equals(expected_field->type()))
          {
            result.valid = false;
            result.errors.push_back("Type mismatch for field '" + actual_field->name() +
              "': got " + actual_field->type()->ToString() + ", expected " + expected_field->type()->ToString());
          }
          if (actual_field->nullable() != expected_field->nullable())
          {
            result.valid = false;
            result.errors.push_back("Nullability mismatch for field '" + actual_field->name() +
              "': got " + std::string(actual_field->nullable() ? "nullable" : "non-null") +
              ", expected " + std::string(expected_field->nullable() ? "nullable" : "non-null"));
          }
        }
      }

      return result;
    }
  }

} // namespace OpenMS

#endif // WITH_PARQUET
