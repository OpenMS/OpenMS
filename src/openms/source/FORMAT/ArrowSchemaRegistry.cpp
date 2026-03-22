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
  } // namespace ArrowSchemaValidation

  // -- ProteinSchema --

  std::shared_ptr<arrow::DataType> ProteinSchema::modificationsType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("position", arrow::int32()),
      arrow::field("modification", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::DataType> ProteinSchema::metavaluesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("value", arrow::utf8()),
      arrow::field("value_type", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::Schema> ProteinSchema::schema()
  {
    return arrow::schema({
      arrow::field(ACCESSION, arrow::utf8(), /*nullable=*/false),
      arrow::field(SCORE, arrow::float64(), /*nullable=*/false),
      arrow::field(RANK, arrow::int32(), /*nullable=*/true),
      arrow::field(COVERAGE, arrow::float64(), /*nullable=*/true),
      arrow::field(SEQUENCE, arrow::utf8(), /*nullable=*/true),
      arrow::field(DESCRIPTION, arrow::utf8(), /*nullable=*/true),
      arrow::field(IS_DECOY, arrow::boolean(), /*nullable=*/true),
      arrow::field(RUN_IDENTIFIER, arrow::utf8(), /*nullable=*/false),
      arrow::field(MODIFICATIONS, modificationsType(), /*nullable=*/true),
      arrow::field(METAVALUES, metavaluesType(), /*nullable=*/false),
    });
  }

  // -- ProteinGroupSchema --

  std::shared_ptr<arrow::DataType> ProteinGroupSchema::floatDataType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("values", arrow::list(arrow::float64()))
    }));
  }

  std::shared_ptr<arrow::DataType> ProteinGroupSchema::stringDataType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("values", arrow::list(arrow::utf8()))
    }));
  }

  std::shared_ptr<arrow::DataType> ProteinGroupSchema::integerDataType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("values", arrow::list(arrow::int64()))
    }));
  }

  std::shared_ptr<arrow::Schema> ProteinGroupSchema::schema()
  {
    return arrow::schema({
      arrow::field(GROUP_TYPE, arrow::utf8(), /*nullable=*/false),
      arrow::field(PROBABILITY, arrow::float64(), /*nullable=*/false),
      arrow::field(ACCESSIONS, arrow::list(arrow::utf8()), /*nullable=*/false),
      arrow::field(RUN_IDENTIFIER, arrow::utf8(), /*nullable=*/false),
      arrow::field(GROUP_INDEX, arrow::int32(), /*nullable=*/false),
      arrow::field(FLOAT_DATA, floatDataType(), /*nullable=*/true),
      arrow::field(STRING_DATA, stringDataType(), /*nullable=*/true),
      arrow::field(INTEGER_DATA, integerDataType(), /*nullable=*/true),
    });
  }

  // -- SearchParamsSchema --

  std::shared_ptr<arrow::DataType> SearchParamsSchema::metavaluesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("value", arrow::utf8()),
      arrow::field("value_type", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::Schema> SearchParamsSchema::schema()
  {
    return arrow::schema({
      arrow::field(RUN_IDENTIFIER, arrow::utf8(), /*nullable=*/false),
      arrow::field(SEARCH_ENGINE, arrow::utf8(), /*nullable=*/false),
      arrow::field(SEARCH_ENGINE_VERSION, arrow::utf8(), /*nullable=*/true),
      arrow::field(INFERENCE_ENGINE, arrow::utf8(), /*nullable=*/true),
      arrow::field(INFERENCE_ENGINE_VERSION, arrow::utf8(), /*nullable=*/true),
      arrow::field(DATE, arrow::timestamp(arrow::TimeUnit::SECOND), /*nullable=*/true),
      arrow::field(SCORE_TYPE, arrow::utf8(), /*nullable=*/false),
      arrow::field(HIGHER_SCORE_BETTER, arrow::boolean(), /*nullable=*/false),
      arrow::field(SIGNIFICANCE_THRESHOLD, arrow::float64(), /*nullable=*/true),
      arrow::field(DB, arrow::utf8(), /*nullable=*/true),
      arrow::field(DB_VERSION, arrow::utf8(), /*nullable=*/true),
      arrow::field(TAXONOMY, arrow::utf8(), /*nullable=*/true),
      arrow::field(CHARGES, arrow::utf8(), /*nullable=*/true),
      arrow::field(MASS_TYPE, arrow::utf8(), /*nullable=*/false),
      arrow::field(PRECURSOR_MASS_TOLERANCE, arrow::float64(), /*nullable=*/false),
      arrow::field(PRECURSOR_MASS_TOLERANCE_PPM, arrow::boolean(), /*nullable=*/false),
      arrow::field(FRAGMENT_MASS_TOLERANCE, arrow::float64(), /*nullable=*/false),
      arrow::field(FRAGMENT_MASS_TOLERANCE_PPM, arrow::boolean(), /*nullable=*/false),
      arrow::field(DIGESTION_ENZYME, arrow::utf8(), /*nullable=*/true),
      arrow::field(ENZYME_TERM_SPECIFICITY, arrow::utf8(), /*nullable=*/true),
      arrow::field(MISSED_CLEAVAGES, arrow::int32(), /*nullable=*/false),
      arrow::field(FIXED_MODIFICATIONS, arrow::list(arrow::utf8()), /*nullable=*/false),
      arrow::field(VARIABLE_MODIFICATIONS, arrow::list(arrow::utf8()), /*nullable=*/false),
      arrow::field(PRIMARY_MS_RUN_PATHS, arrow::list(arrow::utf8()), /*nullable=*/false),
      arrow::field(METAVALUES, metavaluesType(), /*nullable=*/false),
      arrow::field(SP_METAVALUES, metavaluesType(), /*nullable=*/false),
    });
  }

  // -- FeatureSchema --

  std::shared_ptr<arrow::DataType> FeatureSchema::convexHullType()
  {
    auto point_struct_type = arrow::struct_({
      arrow::field("x", arrow::float64()),
      arrow::field("y", arrow::float64())
    });
    return arrow::list(arrow::struct_({
      arrow::field("hull_index", arrow::int32()),
      arrow::field("points", arrow::list(point_struct_type))
    }));
  }

  std::shared_ptr<arrow::DataType> FeatureSchema::metavaluesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("value", arrow::utf8()),
      arrow::field("value_type", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::Schema> FeatureSchema::schema()
  {
    return arrow::schema({
      arrow::field(UNIQUE_ID, arrow::int64(), /*nullable=*/false),
      arrow::field(PARENT_FEATURE_ID, arrow::int64(), /*nullable=*/true),
      arrow::field(DEPTH, arrow::int32(), /*nullable=*/false),
      arrow::field(RT, arrow::float64(), /*nullable=*/false),
      arrow::field(MZ, arrow::float64(), /*nullable=*/false),
      arrow::field(INTENSITY, arrow::float32(), /*nullable=*/false),
      arrow::field(CHARGE, arrow::int32(), /*nullable=*/false),
      arrow::field(QUALITY, arrow::float32(), /*nullable=*/false),
      arrow::field(QUALITY_RT, arrow::float32(), /*nullable=*/false),
      arrow::field(QUALITY_MZ, arrow::float32(), /*nullable=*/false),
      arrow::field(WIDTH, arrow::float32(), /*nullable=*/true),
      arrow::field(RT_BB_MIN, arrow::float64(), /*nullable=*/true),
      arrow::field(RT_BB_MAX, arrow::float64(), /*nullable=*/true),
      arrow::field(MZ_BB_MIN, arrow::float64(), /*nullable=*/true),
      arrow::field(MZ_BB_MAX, arrow::float64(), /*nullable=*/true),
      arrow::field(CONVEX_HULLS, convexHullType(), /*nullable=*/false),
      arrow::field(METAVALUES, metavaluesType(), /*nullable=*/false),
    });
  }

  // -- ConsensusFeatureSchema --

  std::shared_ptr<arrow::DataType> ConsensusFeatureSchema::handlesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("map_index", arrow::int64()),
      arrow::field("unique_id", arrow::int64()),
      arrow::field("rt", arrow::float64()),
      arrow::field("mz", arrow::float64()),
      arrow::field("intensity", arrow::float32()),
      arrow::field("charge", arrow::int32()),
      arrow::field("width", arrow::float32())
    }));
  }

  std::shared_ptr<arrow::DataType> ConsensusFeatureSchema::metavaluesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("value", arrow::utf8()),
      arrow::field("value_type", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::Schema> ConsensusFeatureSchema::schema()
  {
    return arrow::schema({
      arrow::field(UNIQUE_ID, arrow::int64()),
      arrow::field(RT, arrow::float64()),
      arrow::field(MZ, arrow::float64()),
      arrow::field(INTENSITY, arrow::float32()),
      arrow::field(CHARGE, arrow::int32()),
      arrow::field(QUALITY, arrow::float32()),
      arrow::field(WIDTH, arrow::float32()),
      arrow::field(HANDLES, handlesType()),
      arrow::field(METAVALUES, metavaluesType()),
    });
  }

  // -- PSMSchema --

  std::shared_ptr<arrow::DataType> PSMSchema::modificationsType()
  {
    auto position_struct_type = arrow::struct_({
      arrow::field("position", arrow::utf8()),
      arrow::field("scores", arrow::float64())
    });
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("accession", arrow::utf8()),
      arrow::field("positions", arrow::list(position_struct_type))
    }));
  }

  std::shared_ptr<arrow::DataType> PSMSchema::additionalScoresType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("score_name", arrow::utf8()),
      arrow::field("score_value", arrow::float64()),
      arrow::field("higher_better", arrow::boolean())
    }));
  }

  std::shared_ptr<arrow::DataType> PSMSchema::metavaluesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8()),
      arrow::field("value", arrow::utf8()),
      arrow::field("value_type", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::Schema> PSMSchema::schema()
  {
    return arrow::schema({
      arrow::field(SEQUENCE, arrow::utf8()),
      arrow::field(PEPTIDOFORM, arrow::utf8()),
      arrow::field(MODIFICATIONS, modificationsType()),
      arrow::field(PRECURSOR_CHARGE, arrow::int32()),
      arrow::field(POSTERIOR_ERROR_PROBABILITY, arrow::float64()),
      arrow::field(IS_DECOY, arrow::boolean()),
      arrow::field(CALCULATED_MZ, arrow::float64()),
      arrow::field(OBSERVED_MZ, arrow::float64()),
      arrow::field(ADDITIONAL_SCORES, additionalScoresType()),
      arrow::field(PROTEIN_ACCESSIONS, arrow::list(arrow::utf8())),
      arrow::field(PREDICTED_RT, arrow::float64()),
      arrow::field(REFERENCE_FILE_NAME, arrow::utf8()),
      arrow::field(CV_PARAMS, arrow::utf8()),
      arrow::field(SCAN, arrow::int32()),
      arrow::field(RT, arrow::float64()),
      arrow::field(ION_MOBILITY, arrow::float64()),
      arrow::field(SPECTRUM_REFERENCE, arrow::utf8()),
      arrow::field(SCORE, arrow::float64()),
      arrow::field(SCORE_TYPE, arrow::utf8()),
      arrow::field(HIGHER_SCORE_BETTER, arrow::boolean()),
      arrow::field(RANK, arrow::int32()),
      arrow::field(PEPTIDE_IDENTIFICATION_INDEX, arrow::int32()),
      arrow::field(PSM_METAVALUES, metavaluesType()),
      arrow::field(SPECTRUM_METAVALUES, metavaluesType()),
      arrow::field(RUN_IDENTIFIER, arrow::utf8()),
    });
  }

} // namespace OpenMS

#endif // WITH_PARQUET
