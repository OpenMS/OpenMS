// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>
#include <sstream>

namespace OpenMS
{

  namespace ArrowSchemaValidation
  {
    /// Check if two types are compatible, allowing timestamp precision differences.
    /// Parquet may change timestamp resolution (e.g., timestamp[s] -> timestamp[ms]) on round-trip.
    static bool areTypesCompatible(const std::shared_ptr<arrow::DataType>& actual,
                                   const std::shared_ptr<arrow::DataType>& expected)
    {
      if (actual->Equals(expected)) return true;
      // Allow timestamp precision differences (both are timestamp types)
      if (actual->id() == arrow::Type::TIMESTAMP && expected->id() == arrow::Type::TIMESTAMP) return true;
      return false;
    }

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
          if (!areTypesCompatible(actual_field->type(), expected_field->type()))
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
            continue;  // Ignore unknown fields for forward compatibility
          }
          if (!areTypesCompatible(actual_field->type(), expected_field->type()))
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
      arrow::field(RANK, arrow::int32(), /*nullable=*/false),
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

  std::shared_ptr<arrow::DataType> PSMSchema::proteinAccessionsType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("accession", arrow::utf8(), /*nullable=*/false),
      arrow::field("aa_before", arrow::utf8(), /*nullable=*/true),
      arrow::field("aa_after",  arrow::utf8(), /*nullable=*/true),
      arrow::field("start",     arrow::int32(), /*nullable=*/true),
      arrow::field("end",       arrow::int32(), /*nullable=*/true),
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
      arrow::field(PROTEIN_ACCESSIONS, proteinAccessionsType()),
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
      arrow::field(HIT_INDEX, arrow::int32()),
      arrow::field(PEPTIDE_IDENTIFICATION_INDEX, arrow::int32()),
      arrow::field(PSM_METAVALUES, metavaluesType()),
      arrow::field(SPECTRUM_METAVALUES, metavaluesType()),
      arrow::field(RUN_IDENTIFIER, arrow::utf8()),
      arrow::field(MZ_ARRAY, arrow::list(arrow::float32())),
      arrow::field(INTENSITY_ARRAY, arrow::list(arrow::float32())),
      arrow::field(CHARGE_ARRAY, arrow::list(arrow::int32())),
      arrow::field(ION_TYPE_ARRAY, arrow::list(arrow::utf8())),
    });
  }

  // -- QPXPSMSchema (quantms Parquet eXchange format, PSM table) --

  std::shared_ptr<arrow::DataType> QPXPSMSchema::modificationsType()
  {
    auto scores_struct = arrow::struct_({
      arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("score_value", arrow::float64(), /*nullable=*/false),
      arrow::field("higher_better", arrow::boolean())
    });
    auto position_struct = arrow::struct_({
      arrow::field("position", arrow::int32(), /*nullable=*/false),
      arrow::field("amino_acid", arrow::utf8()),
      arrow::field("scores", arrow::list(scores_struct))
    });
    return arrow::list(arrow::struct_({
      arrow::field("name", arrow::utf8(), /*nullable=*/false),
      arrow::field("accession", arrow::utf8()),
      arrow::field("positions", arrow::list(position_struct), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPSMSchema::additionalScoresType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("score_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("score_value", arrow::float64(), /*nullable=*/false),
      arrow::field("higher_better", arrow::boolean())
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPSMSchema::cvParamsType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("cv_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("cv_value", arrow::utf8(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPSMSchema::crossLinksType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("xl_type", arrow::utf8(), /*nullable=*/false),
      arrow::field("partner_sequence", arrow::utf8()),
      arrow::field("partner_peptidoform", arrow::utf8()),
      arrow::field("donor_position", arrow::int32(), /*nullable=*/false),
      arrow::field("acceptor_position", arrow::int32()),
      arrow::field("linker_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("linker_accession", arrow::utf8()),
      arrow::field("linker_mass", arrow::float64(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::Schema> QPXPSMSchema::schema()
  {
    return arrow::schema({
      arrow::field(SEQUENCE, arrow::utf8(), /*nullable=*/false),
      arrow::field(PEPTIDOFORM, arrow::utf8(), /*nullable=*/false),
      arrow::field(MODIFICATIONS, modificationsType()),
      arrow::field(CHARGE, arrow::int16(), /*nullable=*/false),
      arrow::field(POSTERIOR_ERROR_PROBABILITY, arrow::float64()),
      arrow::field(IS_DECOY, arrow::boolean(), /*nullable=*/false),
      arrow::field(CALCULATED_MZ, arrow::float32(), /*nullable=*/false),
      arrow::field(OBSERVED_MZ, arrow::float32(), /*nullable=*/false),
      arrow::field(MASS_ERROR_PPM, arrow::float32()),
      arrow::field(ADDITIONAL_SCORES, additionalScoresType()),
      arrow::field(PREDICTED_RT, arrow::float32()),
      arrow::field(RUN_FILE_NAME, arrow::utf8(), /*nullable=*/false),
      arrow::field(CV_PARAMS, cvParamsType()),
      arrow::field(SCAN, arrow::list(arrow::int32()), /*nullable=*/false),
      arrow::field(RT, arrow::float32()),
      arrow::field(ION_MOBILITY, arrow::float32()),
      arrow::field(MISSED_CLEAVAGES, arrow::int16()),
      arrow::field(PROTEIN_ACCESSIONS, arrow::list(arrow::utf8())),
      arrow::field(CROSS_LINKS, crossLinksType()),
      arrow::field(MZ_ARRAY, arrow::list(arrow::float32())),
      arrow::field(INTENSITY_ARRAY, arrow::list(arrow::float32())),
      arrow::field(CHARGE_ARRAY, arrow::list(arrow::int32())),
      arrow::field(ION_TYPE_ARRAY, arrow::list(arrow::utf8())),
      arrow::field(ION_MOBILITY_ARRAY, arrow::list(arrow::float32())),
    });
  }

  // -- QPXPgSchema (quantms Parquet eXchange format, protein group table) --

  std::shared_ptr<arrow::DataType> QPXPgSchema::intensitiesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("label", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensity", arrow::float32(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::additionalIntensitiesType()
  {
    auto int_pair_type = arrow::struct_({
      arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
    });
    return arrow::list(arrow::struct_({
      arrow::field("label", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensities", arrow::list(int_pair_type), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::peptidesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("protein_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("peptide_count", arrow::int32(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::peptideCountsType()
  {
    return arrow::struct_({
      arrow::field("unique_sequences", arrow::int32(), /*nullable=*/false),
      arrow::field("total_sequences", arrow::int32(), /*nullable=*/false)
    });
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::featureCountsType()
  {
    return arrow::struct_({
      arrow::field("unique_features", arrow::int32(), /*nullable=*/false),
      arrow::field("total_features", arrow::int32(), /*nullable=*/false)
    });
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::additionalScoresType()
  {
    return QPXPSMSchema::additionalScoresType();
  }

  std::shared_ptr<arrow::DataType> QPXPgSchema::cvParamsType()
  {
    return QPXPSMSchema::cvParamsType();
  }

  std::shared_ptr<arrow::Schema> QPXPgSchema::schema()
  {
    return arrow::schema({
      arrow::field(PG_ACCESSIONS, arrow::list(arrow::utf8()), /*nullable=*/false),
      arrow::field(PG_NAMES, arrow::list(arrow::utf8())),
      arrow::field(GG_ACCESSIONS, arrow::list(arrow::utf8())),
      arrow::field(GG_NAMES, arrow::list(arrow::utf8())),
      arrow::field(GG_QVALUE, arrow::float64()),
      arrow::field(ANCHOR_PROTEIN, arrow::utf8(), /*nullable=*/false),
      arrow::field(RUN_FILE_NAME, arrow::utf8(), /*nullable=*/false),
      arrow::field(GLOBAL_QVALUE, arrow::float64()),
      arrow::field(PG_QVALUE, arrow::float64()),
      arrow::field(INTENSITIES, intensitiesType()),
      arrow::field(ADDITIONAL_INTENSITIES, additionalIntensitiesType()),
      arrow::field(IS_DECOY, arrow::boolean(), /*nullable=*/false),
      arrow::field(CONTAMINANT, arrow::boolean()),
      arrow::field(PEPTIDES, peptidesType(), /*nullable=*/false),
      arrow::field(PEPTIDE_COUNTS, peptideCountsType()),
      arrow::field(FEATURE_COUNTS, featureCountsType()),
      arrow::field(SEQUENCE_COVERAGE, arrow::float32()),
      arrow::field(MOLECULAR_WEIGHT, arrow::float32()),
      arrow::field(ADDITIONAL_SCORES, additionalScoresType()),
      arrow::field(CV_PARAMS, cvParamsType()),
    });
  }

  // -- QPXFeatureSchema (quantms Parquet eXchange format) --
  // Shared types delegate to QPXPSMSchema to avoid duplication

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::modificationsType()
  {
    return QPXPSMSchema::modificationsType();
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::additionalScoresType()
  {
    return QPXPSMSchema::additionalScoresType();
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::cvParamsType()
  {
    return QPXPSMSchema::cvParamsType();
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::intensitiesType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("label", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensity", arrow::float32(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::additionalIntensitiesType()
  {
    auto intensity_entry = arrow::struct_({
      arrow::field("intensity_name", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensity_value", arrow::float32(), /*nullable=*/false)
    });
    return arrow::list(arrow::struct_({
      arrow::field("label", arrow::utf8(), /*nullable=*/false),
      arrow::field("intensities", arrow::list(intensity_entry), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::pgAccessionsType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("accession", arrow::utf8(), /*nullable=*/false),
      arrow::field("start", arrow::int32()),
      arrow::field("end", arrow::int32()),
      arrow::field("pre", arrow::utf8()),
      arrow::field("post", arrow::utf8())
    }));
  }

  std::shared_ptr<arrow::DataType> QPXFeatureSchema::pgPositionsType()
  {
    return arrow::list(arrow::struct_({
      arrow::field("protein_accession", arrow::utf8(), /*nullable=*/false),
      arrow::field("start", arrow::int32(), /*nullable=*/false),
      arrow::field("end", arrow::int32(), /*nullable=*/false)
    }));
  }

  std::shared_ptr<arrow::Schema> QPXFeatureSchema::schema()
  {
    return arrow::schema({
      arrow::field(SEQUENCE, arrow::utf8(), /*nullable=*/false),
      arrow::field(PEPTIDOFORM, arrow::utf8(), /*nullable=*/false),
      arrow::field(MODIFICATIONS, modificationsType()),
      arrow::field(CHARGE, arrow::int16(), /*nullable=*/false),
      arrow::field(POSTERIOR_ERROR_PROBABILITY, arrow::float64()),
      arrow::field(IS_DECOY, arrow::boolean(), /*nullable=*/false),
      arrow::field(CALCULATED_MZ, arrow::float32(), /*nullable=*/false),
      arrow::field(OBSERVED_MZ, arrow::float32(), /*nullable=*/false),
      arrow::field(MASS_ERROR_PPM, arrow::float32()),
      arrow::field(ADDITIONAL_SCORES, additionalScoresType()),
      arrow::field(PREDICTED_RT, arrow::float32()),
      arrow::field(RUN_FILE_NAME, arrow::utf8(), /*nullable=*/false),
      arrow::field(CV_PARAMS, cvParamsType()),
      arrow::field(SCAN, arrow::list(arrow::int32()), /*nullable=*/false),
      arrow::field(RT, arrow::float32()),
      arrow::field(ION_MOBILITY, arrow::float32()),
      arrow::field(MISSED_CLEAVAGES, arrow::int16()),
      arrow::field(INTENSITIES, intensitiesType()),
      arrow::field(ADDITIONAL_INTENSITIES, additionalIntensitiesType()),
      arrow::field(PG_ACCESSIONS, pgAccessionsType()),
      arrow::field(ANCHOR_PROTEIN, arrow::utf8(), /*nullable=*/false),
      arrow::field(UNIQUE, arrow::boolean()),
      arrow::field(PG_GLOBAL_QVALUE, arrow::float64()),
      arrow::field(PG_POSITIONS, pgPositionsType()),
      arrow::field(ION_MOBILITY_START, arrow::float32()),
      arrow::field(ION_MOBILITY_STOP, arrow::float32()),
      arrow::field(GG_ACCESSIONS, arrow::list(arrow::utf8())),
      arrow::field(GG_NAMES, arrow::list(arrow::utf8())),
      arrow::field(ID_RUN_FILE_NAME, arrow::utf8()),
      arrow::field(RT_START, arrow::float32()),
      arrow::field(RT_STOP, arrow::float32()),
    });
  }

  // -- SpectraLongSchema --

  std::shared_ptr<arrow::Schema> SpectraLongSchema::schema()
  {
    return arrow::schema({
      arrow::field(MZ, arrow::float64()),
      arrow::field(INTENSITY, arrow::float32()),
      arrow::field(RT, arrow::float32()),
      arrow::field(ION_MOBILITY, arrow::float32()),
      arrow::field(SPECTRUM_INDEX, arrow::uint32()),
      arrow::field(MS_LEVEL, arrow::uint8()),
      arrow::field(NATIVE_ID, arrow::utf8()),
      arrow::field(PRECURSOR_MZ, arrow::float64()),
      arrow::field(PRECURSOR_CHARGE, arrow::int16()),
      arrow::field(PRECURSOR_INTENSITY, arrow::float32()),
      arrow::field(ISOLATION_LOWER, arrow::float64()),
      arrow::field(ISOLATION_UPPER, arrow::float64()),
    });
  }

  // -- SpectraSemiWideSchema --

  std::shared_ptr<arrow::Schema> SpectraSemiWideSchema::schema()
  {
    return arrow::schema({
      arrow::field(SPECTRUM_INDEX, arrow::uint32()),
      arrow::field(RT, arrow::float32()),
      arrow::field(MS_LEVEL, arrow::uint8()),
      arrow::field(NATIVE_ID, arrow::utf8()),
      arrow::field(MZ, arrow::list(arrow::float64())),
      arrow::field(INTENSITY, arrow::list(arrow::float32())),
      arrow::field(ION_MOBILITY, arrow::list(arrow::float32())),
      arrow::field(PRECURSOR_MZ, arrow::float64()),
      arrow::field(PRECURSOR_CHARGE, arrow::int16()),
      arrow::field(PRECURSOR_INTENSITY, arrow::float32()),
      arrow::field(ISOLATION_LOWER, arrow::float64()),
      arrow::field(ISOLATION_UPPER, arrow::float64()),
    });
  }

  // -- ChromatogramSchema --

  std::shared_ptr<arrow::Schema> ChromatogramSchema::schema()
  {
    return arrow::schema({
      arrow::field(RT, arrow::float64()),
      arrow::field(INTENSITY, arrow::float32()),
      arrow::field(CHROMATOGRAM_INDEX, arrow::uint32()),
      arrow::field(NATIVE_ID, arrow::utf8()),
      arrow::field(PRECURSOR_MZ, arrow::float64()),
      arrow::field(PRODUCT_MZ, arrow::float64()),
    });
  }

  // -- ChromatogramSemiWideSchema --

  std::shared_ptr<arrow::Schema> ChromatogramSemiWideSchema::schema()
  {
    return arrow::schema({
      arrow::field(CHROMATOGRAM_INDEX, arrow::uint32()),
      arrow::field(NATIVE_ID, arrow::utf8()),
      arrow::field(RT, arrow::list(arrow::float64())),
      arrow::field(INTENSITY, arrow::list(arrow::float32())),
      arrow::field(PRECURSOR_MZ, arrow::float64()),
      arrow::field(PRODUCT_MZ, arrow::float64()),
    });
  }

  // -- OSWPrecursorSchema --

  std::shared_ptr<arrow::Schema> OSWPrecursorSchema::schema()
  {
    return arrow::schema({
      arrow::field(PRECURSOR_ID, arrow::int64()),
      arrow::field(PRECURSOR_MZ, arrow::float64()),
      arrow::field(CHARGE, arrow::int32()),
      arrow::field(LIBRARY_RT, arrow::float64()),
      arrow::field(LIBRARY_DRIFT_TIME, arrow::float64()),
      arrow::field(DECOY, arrow::boolean()),
      arrow::field(TRAML_ID, arrow::utf8()),
      arrow::field(MODIFIED_SEQUENCE, arrow::utf8()),
      arrow::field(UNMODIFIED_SEQUENCE, arrow::utf8()),
      arrow::field(PROTEIN_ACCESSIONS, arrow::utf8()),
    });
  }

  // -- OSWTransitionSchema --

  std::shared_ptr<arrow::Schema> OSWTransitionSchema::schema()
  {
    return arrow::schema({
      arrow::field(TRANSITION_ID, arrow::int64()),
      arrow::field(PRECURSOR_ID, arrow::int64()),
      arrow::field(TRAML_ID, arrow::utf8()),
      arrow::field(PRODUCT_MZ, arrow::float64()),
      arrow::field(CHARGE, arrow::int32()),
      arrow::field(TYPE, arrow::utf8()),
      arrow::field(ANNOTATION, arrow::utf8()),
      arrow::field(ORDINAL, arrow::int32()),
      arrow::field(DETECTING, arrow::boolean()),
      arrow::field(IDENTIFYING, arrow::boolean()),
      arrow::field(QUANTIFYING, arrow::boolean()),
      arrow::field(LIBRARY_INTENSITY, arrow::float64()),
      arrow::field(DECOY, arrow::boolean()),
    });
  }

  // -- OSWFeaturePrecursorSchema --

  std::shared_ptr<arrow::Schema> OSWFeaturePrecursorSchema::schema()
  {
    return arrow::schema({
      arrow::field(FEATURE_ID, arrow::int64()),
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(PRECURSOR_ISOTOPE, arrow::int32()),
      arrow::field(PRECURSOR_AREA_INTENSITY, arrow::float64()),
      arrow::field(PRECURSOR_APEX_INTENSITY, arrow::float64()),
    });
  }

  // -- OSWRunSchema --

  std::shared_ptr<arrow::Schema> OSWRunSchema::schema()
  {
    return arrow::schema({
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(FILENAME, arrow::utf8()),
    });
  }

  // -- OSWFeatureSchema --

  std::shared_ptr<arrow::Schema> OSWFeatureSchema::schema()
  {
    return arrow::schema({
      arrow::field(FEATURE_ID, arrow::int64()),
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(PRECURSOR_ID, arrow::int64()),
      arrow::field(EXP_RT, arrow::float64()),
      arrow::field(EXP_IM, arrow::float64()),
      arrow::field(NORM_RT, arrow::float64()),
      arrow::field(DELTA_RT, arrow::float64()),
      arrow::field(LEFT_WIDTH, arrow::float64()),
      arrow::field(RIGHT_WIDTH, arrow::float64()),
      arrow::field(EXP_IM_LEFTWIDTH, arrow::float64()),
      arrow::field(EXP_IM_RIGHTWIDTH, arrow::float64()),
      arrow::field(MS1_AREA_INTENSITY, arrow::float64()),
      arrow::field(MS1_APEX_INTENSITY, arrow::float64()),
      arrow::field(MS1_EXP_IM, arrow::float64()),
      arrow::field(MS1_DELTA_IM, arrow::float64()),
      arrow::field(VAR_MS1_MASSDEV_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_IM_MS1_DELTA_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_MI_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_MI_CONTRAST_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_MI_COMBINED_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_ISOTOPE_CORRELATION_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_ISOTOPE_OVERLAP_SCORE, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_COELUTION, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_COELUTION_CONTRAST, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_COELUTION_COMBINED, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_SHAPE, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_SHAPE_CONTRAST, arrow::float64()),
      arrow::field(VAR_MS1_XCORR_SHAPE_COMBINED, arrow::float64()),
      arrow::field(MS2_AREA_INTENSITY, arrow::float64()),
      arrow::field(MS2_TOTAL_AREA_INTENSITY, arrow::float64()),
      arrow::field(MS2_APEX_INTENSITY, arrow::float64()),
      arrow::field(MS2_EXP_IM, arrow::float64()),
      arrow::field(MS2_EXP_IM_LEFTWIDTH, arrow::float64()),
      arrow::field(MS2_EXP_IM_RIGHTWIDTH, arrow::float64()),
      arrow::field(MS2_DELTA_IM, arrow::float64()),
      arrow::field(MS2_TOTAL_MI, arrow::float64()),
      arrow::field(VAR_MS2_BSERIES_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_DOTPROD_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_INTENSITY_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_ISOTOPE_CORRELATION_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_ISOTOPE_OVERLAP_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_CORR, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_DOTPROD, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_MANHATTAN, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_RMSD, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_ROOTMEANSQUARE, arrow::float64()),
      arrow::field(VAR_MS2_LIBRARY_SANGLE, arrow::float64()),
      arrow::field(VAR_MS2_LOG_SN_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_MANHATTAN_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_MASSDEV_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_MASSDEV_SCORE_WEIGHTED, arrow::float64()),
      arrow::field(VAR_MS2_MI_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_MI_WEIGHTED_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_MI_RATIO_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_NORM_RT_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_XCORR_COELUTION, arrow::float64()),
      arrow::field(VAR_MS2_XCORR_COELUTION_WEIGHTED, arrow::float64()),
      arrow::field(VAR_MS2_XCORR_SHAPE, arrow::float64()),
      arrow::field(VAR_MS2_XCORR_SHAPE_WEIGHTED, arrow::float64()),
      arrow::field(VAR_MS2_YSERIES_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_ELUTION_MODEL_FIT_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_IM_XCORR_SHAPE, arrow::float64()),
      arrow::field(VAR_MS2_IM_XCORR_COELUTION, arrow::float64()),
      arrow::field(VAR_MS2_IM_DELTA_SCORE, arrow::float64()),
      arrow::field(VAR_MS2_IM_LOG_INTENSITY, arrow::float64()),
    });
  }

  // -- OSWFeatureTransitionSchema --

  std::shared_ptr<arrow::Schema> OSWFeatureTransitionSchema::schema()
  {
    return arrow::schema({
      arrow::field(FEATURE_ID, arrow::int64()),
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(TRANSITION_ID, arrow::int64()),
      arrow::field(AREA_INTENSITY, arrow::float64()),
      arrow::field(TOTAL_AREA_INTENSITY, arrow::float64()),
      arrow::field(APEX_INTENSITY, arrow::float64()),
      arrow::field(APEX_RT, arrow::float64()),
      arrow::field(RT_FWHM, arrow::float64()),
      arrow::field(MASSERROR_PPM, arrow::float64()),
      arrow::field(TOTAL_MI, arrow::float64()),
      arrow::field(VAR_INTENSITY_SCORE, arrow::float64()),
      arrow::field(VAR_INTENSITY_RATIO_SCORE, arrow::float64()),
      arrow::field(VAR_LOG_INTENSITY, arrow::float64()),
      arrow::field(VAR_XCORR_COELUTION, arrow::float64()),
      arrow::field(VAR_XCORR_SHAPE, arrow::float64()),
      arrow::field(VAR_LOG_SN_SCORE, arrow::float64()),
      arrow::field(VAR_MASSDEV_SCORE, arrow::float64()),
      arrow::field(VAR_MI_SCORE, arrow::float64()),
      arrow::field(VAR_MI_RATIO_SCORE, arrow::float64()),
      arrow::field(VAR_ISOTOPE_CORRELATION_SCORE, arrow::float64()),
      arrow::field(VAR_ISOTOPE_OVERLAP_SCORE, arrow::float64()),
      arrow::field(EXP_IM, arrow::float64()),
      arrow::field(EXP_IM_LEFTWIDTH, arrow::float64()),
      arrow::field(EXP_IM_RIGHTWIDTH, arrow::float64()),
      arrow::field(DELTA_IM, arrow::float64()),
      arrow::field(VAR_IM_DELTA_SCORE, arrow::float64()),
      arrow::field(VAR_IM_LOG_INTENSITY, arrow::float64()),
      arrow::field(VAR_IM_XCORR_COELUTION_CONTRAST, arrow::float64()),
      arrow::field(VAR_IM_XCORR_SHAPE_CONTRAST, arrow::float64()),
      arrow::field(VAR_IM_XCORR_COELUTION_COMBINED, arrow::float64()),
      arrow::field(VAR_IM_XCORR_SHAPE_COMBINED, arrow::float64()),
      arrow::field(START_POSITION_AT_5, arrow::float64()),
      arrow::field(END_POSITION_AT_5, arrow::float64()),
      arrow::field(START_POSITION_AT_10, arrow::float64()),
      arrow::field(END_POSITION_AT_10, arrow::float64()),
      arrow::field(START_POSITION_AT_50, arrow::float64()),
      arrow::field(END_POSITION_AT_50, arrow::float64()),
      arrow::field(TOTAL_WIDTH, arrow::float64()),
      arrow::field(TAILING_FACTOR, arrow::float64()),
      arrow::field(ASYMMETRY_FACTOR, arrow::float64()),
      arrow::field(SLOPE_OF_BASELINE, arrow::float64()),
      arrow::field(BASELINE_DELTA_2_HEIGHT, arrow::float64()),
      arrow::field(POINTS_ACROSS_BASELINE, arrow::float64()),
      arrow::field(POINTS_ACROSS_HALF_HEIGHT, arrow::float64()),
    });
  }

  // -- XICSchema --

  std::shared_ptr<arrow::Schema> XICSchema::schema()
  {
    return arrow::schema({
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(SOURCE_FILE, arrow::utf8()),
      arrow::field(MS_LEVEL, arrow::int64()),
      arrow::field(PRECURSOR_ID, arrow::int64()),
      arrow::field(TRANSITION_ID, arrow::int64()),
      arrow::field(MODIFIED_SEQUENCE, arrow::utf8()),
      arrow::field(PRECURSOR_CHARGE, arrow::int64()),
      arrow::field(PRODUCT_CHARGE, arrow::int64()),
      arrow::field(DETECTING_TRANSITION, arrow::int64()),
      arrow::field(PRECURSOR_DECOY, arrow::int64()),
      arrow::field(PRODUCT_DECOY, arrow::int64()),
      arrow::field(TRANSITION_ORDINAL, arrow::int64()),
      arrow::field(TRANSITION_TYPE, arrow::utf8()),
      arrow::field(ANNOTATION, arrow::utf8()),
      arrow::field(RT_DATA, arrow::binary()),
      arrow::field(INTENSITY_DATA, arrow::binary()),
      arrow::field(RT_COMPRESSION, arrow::int64()),
      arrow::field(INTENSITY_COMPRESSION, arrow::int64()),
    });
  }

  // -- XIMSchema --

  std::shared_ptr<arrow::Schema> XIMSchema::schema()
  {
    return arrow::schema({
      arrow::field(RUN_ID, arrow::int64()),
      arrow::field(SOURCE_FILE, arrow::utf8()),
      arrow::field(MS_LEVEL, arrow::int64()),
      arrow::field(MOBILOGRAM_TYPE, arrow::utf8()),
      arrow::field(PRECURSOR_ID, arrow::int64()),
      arrow::field(TRANSITION_ID, arrow::int64()),
      arrow::field(FEATURE_ID, arrow::int64()),
      arrow::field(FEATURE_RT, arrow::float64()),
      arrow::field(MODIFIED_SEQUENCE, arrow::utf8()),
      arrow::field(PRECURSOR_CHARGE, arrow::int64()),
      arrow::field(PRODUCT_CHARGE, arrow::int64()),
      arrow::field(DETECTING_TRANSITION, arrow::int64()),
      arrow::field(PRECURSOR_DECOY, arrow::int64()),
      arrow::field(PRODUCT_DECOY, arrow::int64()),
      arrow::field(TRANSITION_ORDINAL, arrow::int64()),
      arrow::field(TRANSITION_TYPE, arrow::utf8()),
      arrow::field(ANNOTATION, arrow::utf8()),
      arrow::field(MOBILITY_DATA, arrow::binary()),
      arrow::field(INTENSITY_DATA, arrow::binary()),
      arrow::field(MOBILITY_COMPRESSION, arrow::int64()),
      arrow::field(INTENSITY_COMPRESSION, arrow::int64()),
    });
  }

} // namespace OpenMS
