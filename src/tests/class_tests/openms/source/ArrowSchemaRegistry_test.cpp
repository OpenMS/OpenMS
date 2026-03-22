#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>
#include <arrow/json/from_string.h>

using namespace OpenMS;

START_TEST(ArrowSchemaRegistry, "$Id$")

START_SECTION(ValidationResult::toString())
{
  ArrowSchemaValidation::ValidationResult result;
  result.valid = false;
  result.errors = {"error1", "error2"};
  TEST_STRING_EQUAL(result.toString(), "error1; error2")
}
END_SECTION

START_SECTION(validate - Strict mode - exact match passes)
{
  auto schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x","y"])").ValueOrDie();
  auto arr_b = arrow::json::ArrayFromJSONString(arrow::float64(), "[1.0, 2.0]").ValueOrDie();
  auto table = arrow::Table::Make(schema, {arr_a, arr_b});

  auto result = ArrowSchemaValidation::validate(table, schema);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(validate - Strict mode - missing field)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
  TEST_EQUAL(result.errors.empty(), false)
}
END_SECTION

START_SECTION(validate - Strict mode - extra field)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto arr_b = arrow::json::ArrayFromJSONString(arrow::float64(), "[1.0]").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_b});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate - Strict mode - type mismatch)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::float64(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::int32(), false),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::int32(), "[1]").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate - Strict mode - nullability mismatch)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), true),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate - Subset mode - valid subset passes)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("b", arrow::float64(), true),
    arrow::field("c", arrow::int32(), true),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("c", arrow::int32(), true),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto arr_c = arrow::json::ArrayFromJSONString(arrow::int32(), "[1]").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_c});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
}
END_SECTION

START_SECTION(validate - Subset mode - unknown field fails)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("z", arrow::float64(), true),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto arr_z = arrow::json::ArrayFromJSONString(arrow::float64(), "[1.0]").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_z});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate - Subset mode - type mismatch fails)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::float64(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::int32(), false),
  });
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::int32(), "[1]").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, false)
}
END_SECTION

START_SECTION(validate - metadata is ignored)
{
  auto metadata = arrow::key_value_metadata({"key"}, {"value"});
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  })->WithMetadata(metadata);
  auto arr_a = arrow::json::ArrayFromJSONString(arrow::utf8(), R"(["x"])").ValueOrDie();
  auto table = arrow::Table::Make(actual_schema, {arr_a});

  auto result = ArrowSchemaValidation::validate(table, expected);
  TEST_EQUAL(result.valid, true)
}
END_SECTION

// ========== ProteinSchema ==========

START_SECTION(ProteinSchema::schema() returns non-null with 10 fields)
{
  auto s = ProteinSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 10)
}
END_SECTION

START_SECTION(ProteinSchema column name constants)
{
  TEST_STRING_EQUAL(ProteinSchema::ACCESSION, "accession")
  TEST_STRING_EQUAL(ProteinSchema::SCORE, "score")
  TEST_STRING_EQUAL(ProteinSchema::RANK, "rank")
  TEST_STRING_EQUAL(ProteinSchema::COVERAGE, "coverage")
  TEST_STRING_EQUAL(ProteinSchema::SEQUENCE, "sequence")
  TEST_STRING_EQUAL(ProteinSchema::DESCRIPTION, "description")
  TEST_STRING_EQUAL(ProteinSchema::IS_DECOY, "is_decoy")
  TEST_STRING_EQUAL(ProteinSchema::RUN_IDENTIFIER, "run_identifier")
  TEST_STRING_EQUAL(ProteinSchema::MODIFICATIONS, "modifications")
  TEST_STRING_EQUAL(ProteinSchema::METAVALUES, "metavalues")
}
END_SECTION

START_SECTION(ProteinSchema field types and nullability)
{
  auto s = ProteinSchema::schema();
  // accession: utf8, non-null
  TEST_EQUAL(s->field(0)->name(), "accession")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // score: float64, non-null
  TEST_EQUAL(s->field(1)->name(), "score")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(1)->nullable(), false)
  // rank: int32, nullable
  TEST_EQUAL(s->field(2)->name(), "rank")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // coverage: float64, nullable
  TEST_EQUAL(s->field(3)->name(), "coverage")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // sequence: utf8, nullable
  TEST_EQUAL(s->field(4)->name(), "sequence")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // description: utf8, nullable
  TEST_EQUAL(s->field(5)->name(), "description")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // is_decoy: boolean, nullable
  TEST_EQUAL(s->field(6)->name(), "is_decoy")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // run_identifier: utf8, non-null
  TEST_EQUAL(s->field(7)->name(), "run_identifier")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(7)->nullable(), false)
  // modifications: list<struct>, nullable
  TEST_EQUAL(s->field(8)->name(), "modifications")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(8)->nullable(), true)
  TEST_EQUAL(s->field(8)->type()->Equals(ProteinSchema::modificationsType()), true)
  // metavalues: list<struct>, non-null
  TEST_EQUAL(s->field(9)->name(), "metavalues")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(9)->nullable(), false)
  TEST_EQUAL(s->field(9)->type()->Equals(ProteinSchema::metavaluesType()), true)
}
END_SECTION

START_SECTION(ProteinSchema validation with table)
{
  auto s = ProteinSchema::schema();
  // Build zero-length arrays for each field
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== ProteinGroupSchema ==========

START_SECTION(ProteinGroupSchema::schema() returns non-null with 8 fields)
{
  auto s = ProteinGroupSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 8)
}
END_SECTION

START_SECTION(ProteinGroupSchema column name constants)
{
  TEST_STRING_EQUAL(ProteinGroupSchema::GROUP_TYPE, "group_type")
  TEST_STRING_EQUAL(ProteinGroupSchema::PROBABILITY, "probability")
  TEST_STRING_EQUAL(ProteinGroupSchema::ACCESSIONS, "accessions")
  TEST_STRING_EQUAL(ProteinGroupSchema::RUN_IDENTIFIER, "run_identifier")
  TEST_STRING_EQUAL(ProteinGroupSchema::GROUP_INDEX, "group_index")
  TEST_STRING_EQUAL(ProteinGroupSchema::FLOAT_DATA, "float_data")
  TEST_STRING_EQUAL(ProteinGroupSchema::STRING_DATA, "string_data")
  TEST_STRING_EQUAL(ProteinGroupSchema::INTEGER_DATA, "integer_data")
}
END_SECTION

START_SECTION(ProteinGroupSchema field types and nullability)
{
  auto s = ProteinGroupSchema::schema();
  // group_type: utf8, non-null
  TEST_EQUAL(s->field(0)->name(), "group_type")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // probability: float64, non-null
  TEST_EQUAL(s->field(1)->name(), "probability")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(1)->nullable(), false)
  // accessions: list<utf8>, non-null
  TEST_EQUAL(s->field(2)->name(), "accessions")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), false)
  // run_identifier: utf8, non-null
  TEST_EQUAL(s->field(3)->name(), "run_identifier")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(3)->nullable(), false)
  // group_index: int32, non-null
  TEST_EQUAL(s->field(4)->name(), "group_index")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(4)->nullable(), false)
  // float_data: list<struct{name:utf8, values:list<float64>}>, nullable
  TEST_EQUAL(s->field(5)->name(), "float_data")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(5)->nullable(), true)
  TEST_EQUAL(s->field(5)->type()->Equals(ProteinGroupSchema::floatDataType()), true)
  // string_data: list<struct{name:utf8, values:list<utf8>}>, nullable
  TEST_EQUAL(s->field(6)->name(), "string_data")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(6)->nullable(), true)
  TEST_EQUAL(s->field(6)->type()->Equals(ProteinGroupSchema::stringDataType()), true)
  // integer_data: list<struct{name:utf8, values:list<int64>}>, nullable
  TEST_EQUAL(s->field(7)->name(), "integer_data")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(7)->nullable(), true)
  TEST_EQUAL(s->field(7)->type()->Equals(ProteinGroupSchema::integerDataType()), true)
}
END_SECTION

START_SECTION(ProteinGroupSchema validation with table)
{
  auto s = ProteinGroupSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== SearchParamsSchema ==========

START_SECTION(SearchParamsSchema::schema() returns non-null with 26 fields)
{
  auto s = SearchParamsSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 26)
}
END_SECTION

START_SECTION(SearchParamsSchema column name constants)
{
  TEST_STRING_EQUAL(SearchParamsSchema::RUN_IDENTIFIER, "run_identifier")
  TEST_STRING_EQUAL(SearchParamsSchema::SEARCH_ENGINE, "search_engine")
  TEST_STRING_EQUAL(SearchParamsSchema::SEARCH_ENGINE_VERSION, "search_engine_version")
  TEST_STRING_EQUAL(SearchParamsSchema::INFERENCE_ENGINE, "inference_engine")
  TEST_STRING_EQUAL(SearchParamsSchema::INFERENCE_ENGINE_VERSION, "inference_engine_version")
  TEST_STRING_EQUAL(SearchParamsSchema::DATE, "date")
  TEST_STRING_EQUAL(SearchParamsSchema::SCORE_TYPE, "score_type")
  TEST_STRING_EQUAL(SearchParamsSchema::HIGHER_SCORE_BETTER, "higher_score_better")
  TEST_STRING_EQUAL(SearchParamsSchema::SIGNIFICANCE_THRESHOLD, "significance_threshold")
  TEST_STRING_EQUAL(SearchParamsSchema::DB, "db")
  TEST_STRING_EQUAL(SearchParamsSchema::DB_VERSION, "db_version")
  TEST_STRING_EQUAL(SearchParamsSchema::TAXONOMY, "taxonomy")
  TEST_STRING_EQUAL(SearchParamsSchema::CHARGES, "charges")
  TEST_STRING_EQUAL(SearchParamsSchema::MASS_TYPE, "mass_type")
  TEST_STRING_EQUAL(SearchParamsSchema::PRECURSOR_MASS_TOLERANCE, "precursor_mass_tolerance")
  TEST_STRING_EQUAL(SearchParamsSchema::PRECURSOR_MASS_TOLERANCE_PPM, "precursor_mass_tolerance_ppm")
  TEST_STRING_EQUAL(SearchParamsSchema::FRAGMENT_MASS_TOLERANCE, "fragment_mass_tolerance")
  TEST_STRING_EQUAL(SearchParamsSchema::FRAGMENT_MASS_TOLERANCE_PPM, "fragment_mass_tolerance_ppm")
  TEST_STRING_EQUAL(SearchParamsSchema::DIGESTION_ENZYME, "digestion_enzyme")
  TEST_STRING_EQUAL(SearchParamsSchema::ENZYME_TERM_SPECIFICITY, "enzyme_term_specificity")
  TEST_STRING_EQUAL(SearchParamsSchema::MISSED_CLEAVAGES, "missed_cleavages")
  TEST_STRING_EQUAL(SearchParamsSchema::FIXED_MODIFICATIONS, "fixed_modifications")
  TEST_STRING_EQUAL(SearchParamsSchema::VARIABLE_MODIFICATIONS, "variable_modifications")
  TEST_STRING_EQUAL(SearchParamsSchema::PRIMARY_MS_RUN_PATHS, "primary_ms_run_paths")
  TEST_STRING_EQUAL(SearchParamsSchema::METAVALUES, "metavalues")
  TEST_STRING_EQUAL(SearchParamsSchema::SP_METAVALUES, "sp_metavalues")
}
END_SECTION

START_SECTION(SearchParamsSchema field types and nullability)
{
  auto s = SearchParamsSchema::schema();
  // run_identifier: utf8, non-null
  TEST_EQUAL(s->field(0)->name(), "run_identifier")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // search_engine: utf8, non-null
  TEST_EQUAL(s->field(1)->name(), "search_engine")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), false)
  // search_engine_version: utf8, nullable
  TEST_EQUAL(s->field(2)->name(), "search_engine_version")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // inference_engine: utf8, nullable
  TEST_EQUAL(s->field(3)->name(), "inference_engine")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // inference_engine_version: utf8, nullable
  TEST_EQUAL(s->field(4)->name(), "inference_engine_version")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // date: timestamp[s], nullable
  TEST_EQUAL(s->field(5)->name(), "date")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::TIMESTAMP)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // score_type: utf8, non-null
  TEST_EQUAL(s->field(6)->name(), "score_type")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(6)->nullable(), false)
  // higher_score_better: boolean, non-null
  TEST_EQUAL(s->field(7)->name(), "higher_score_better")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(7)->nullable(), false)
  // significance_threshold: float64, nullable
  TEST_EQUAL(s->field(8)->name(), "significance_threshold")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // db: utf8, nullable
  TEST_EQUAL(s->field(9)->name(), "db")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(9)->nullable(), true)
  // db_version: utf8, nullable
  TEST_EQUAL(s->field(10)->name(), "db_version")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // taxonomy: utf8, nullable
  TEST_EQUAL(s->field(11)->name(), "taxonomy")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(11)->nullable(), true)
  // charges: utf8, nullable
  TEST_EQUAL(s->field(12)->name(), "charges")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(12)->nullable(), true)
  // mass_type: utf8, non-null
  TEST_EQUAL(s->field(13)->name(), "mass_type")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(13)->nullable(), false)
  // precursor_mass_tolerance: float64, non-null
  TEST_EQUAL(s->field(14)->name(), "precursor_mass_tolerance")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(14)->nullable(), false)
  // precursor_mass_tolerance_ppm: boolean, non-null
  TEST_EQUAL(s->field(15)->name(), "precursor_mass_tolerance_ppm")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(15)->nullable(), false)
  // fragment_mass_tolerance: float64, non-null
  TEST_EQUAL(s->field(16)->name(), "fragment_mass_tolerance")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(16)->nullable(), false)
  // fragment_mass_tolerance_ppm: boolean, non-null
  TEST_EQUAL(s->field(17)->name(), "fragment_mass_tolerance_ppm")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(17)->nullable(), false)
  // digestion_enzyme: utf8, nullable
  TEST_EQUAL(s->field(18)->name(), "digestion_enzyme")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(18)->nullable(), true)
  // enzyme_term_specificity: utf8, nullable
  TEST_EQUAL(s->field(19)->name(), "enzyme_term_specificity")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(19)->nullable(), true)
  // missed_cleavages: int32, non-null
  TEST_EQUAL(s->field(20)->name(), "missed_cleavages")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(20)->nullable(), false)
  // fixed_modifications: list<utf8>, non-null
  TEST_EQUAL(s->field(21)->name(), "fixed_modifications")
  TEST_EQUAL(s->field(21)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(21)->nullable(), false)
  // variable_modifications: list<utf8>, non-null
  TEST_EQUAL(s->field(22)->name(), "variable_modifications")
  TEST_EQUAL(s->field(22)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(22)->nullable(), false)
  // primary_ms_run_paths: list<utf8>, non-null
  TEST_EQUAL(s->field(23)->name(), "primary_ms_run_paths")
  TEST_EQUAL(s->field(23)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(23)->nullable(), false)
  // metavalues: list<struct>, non-null
  TEST_EQUAL(s->field(24)->name(), "metavalues")
  TEST_EQUAL(s->field(24)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(24)->nullable(), false)
  TEST_EQUAL(s->field(24)->type()->Equals(SearchParamsSchema::metavaluesType()), true)
  // sp_metavalues: list<struct>, non-null
  TEST_EQUAL(s->field(25)->name(), "sp_metavalues")
  TEST_EQUAL(s->field(25)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(25)->nullable(), false)
  TEST_EQUAL(s->field(25)->type()->Equals(SearchParamsSchema::metavaluesType()), true)
}
END_SECTION

START_SECTION(SearchParamsSchema validation with table)
{
  auto s = SearchParamsSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== FeatureSchema ==========

START_SECTION(FeatureSchema::schema() returns non-null with 17 fields)
{
  auto s = FeatureSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 17)
}
END_SECTION

START_SECTION(FeatureSchema column name constants)
{
  TEST_STRING_EQUAL(FeatureSchema::UNIQUE_ID, "unique_id")
  TEST_STRING_EQUAL(FeatureSchema::PARENT_FEATURE_ID, "parent_feature_id")
  TEST_STRING_EQUAL(FeatureSchema::DEPTH, "depth")
  TEST_STRING_EQUAL(FeatureSchema::RT, "rt")
  TEST_STRING_EQUAL(FeatureSchema::MZ, "mz")
  TEST_STRING_EQUAL(FeatureSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(FeatureSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(FeatureSchema::QUALITY, "quality")
  TEST_STRING_EQUAL(FeatureSchema::QUALITY_RT, "quality_rt")
  TEST_STRING_EQUAL(FeatureSchema::QUALITY_MZ, "quality_mz")
  TEST_STRING_EQUAL(FeatureSchema::WIDTH, "width")
  TEST_STRING_EQUAL(FeatureSchema::RT_BB_MIN, "rt_bb_min")
  TEST_STRING_EQUAL(FeatureSchema::RT_BB_MAX, "rt_bb_max")
  TEST_STRING_EQUAL(FeatureSchema::MZ_BB_MIN, "mz_bb_min")
  TEST_STRING_EQUAL(FeatureSchema::MZ_BB_MAX, "mz_bb_max")
  TEST_STRING_EQUAL(FeatureSchema::CONVEX_HULLS, "convex_hulls")
  TEST_STRING_EQUAL(FeatureSchema::METAVALUES, "metavalues")
}
END_SECTION

START_SECTION(FeatureSchema field types and nullability)
{
  auto s = FeatureSchema::schema();
  // unique_id: int64, non-null
  TEST_EQUAL(s->field(0)->name(), "unique_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // parent_feature_id: int64, nullable
  TEST_EQUAL(s->field(1)->name(), "parent_feature_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // depth: int32, non-null
  TEST_EQUAL(s->field(2)->name(), "depth")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(2)->nullable(), false)
  // rt: float64, non-null
  TEST_EQUAL(s->field(3)->name(), "rt")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), false)
  // mz: float64, non-null
  TEST_EQUAL(s->field(4)->name(), "mz")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), false)
  // intensity: float32, non-null
  TEST_EQUAL(s->field(5)->name(), "intensity")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(5)->nullable(), false)
  // charge: int32, non-null
  TEST_EQUAL(s->field(6)->name(), "charge")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(6)->nullable(), false)
  // quality: float32, non-null
  TEST_EQUAL(s->field(7)->name(), "quality")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(7)->nullable(), false)
  // quality_rt: float32, non-null
  TEST_EQUAL(s->field(8)->name(), "quality_rt")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(8)->nullable(), false)
  // quality_mz: float32, non-null
  TEST_EQUAL(s->field(9)->name(), "quality_mz")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(9)->nullable(), false)
  // width: float32, nullable
  TEST_EQUAL(s->field(10)->name(), "width")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // rt_bb_min: float64, nullable
  TEST_EQUAL(s->field(11)->name(), "rt_bb_min")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(11)->nullable(), true)
  // rt_bb_max: float64, nullable
  TEST_EQUAL(s->field(12)->name(), "rt_bb_max")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(12)->nullable(), true)
  // mz_bb_min: float64, nullable
  TEST_EQUAL(s->field(13)->name(), "mz_bb_min")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(13)->nullable(), true)
  // mz_bb_max: float64, nullable
  TEST_EQUAL(s->field(14)->name(), "mz_bb_max")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(14)->nullable(), true)
  // convex_hulls: list<struct>, non-null
  TEST_EQUAL(s->field(15)->name(), "convex_hulls")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(15)->nullable(), false)
  TEST_EQUAL(s->field(15)->type()->Equals(FeatureSchema::convexHullType()), true)
  // metavalues: list<struct>, non-null
  TEST_EQUAL(s->field(16)->name(), "metavalues")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(16)->nullable(), false)
  TEST_EQUAL(s->field(16)->type()->Equals(FeatureSchema::metavaluesType()), true)
}
END_SECTION

START_SECTION(FeatureSchema validation with table)
{
  auto s = FeatureSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== ConsensusFeatureSchema ==========

START_SECTION(ConsensusFeatureSchema::schema() returns non-null with 9 fields)
{
  auto s = ConsensusFeatureSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 9)
}
END_SECTION

START_SECTION(ConsensusFeatureSchema column name constants)
{
  TEST_STRING_EQUAL(ConsensusFeatureSchema::UNIQUE_ID, "unique_id")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::RT, "rt")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::MZ, "mz")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::QUALITY, "quality")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::WIDTH, "width")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::HANDLES, "handles")
  TEST_STRING_EQUAL(ConsensusFeatureSchema::METAVALUES, "metavalues")
}
END_SECTION

START_SECTION(ConsensusFeatureSchema field types and nullability)
{
  auto s = ConsensusFeatureSchema::schema();
  // unique_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "unique_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // rt: float64, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "rt")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // mz: float64, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "mz")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // intensity: float32, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "intensity")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // charge: int32, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "charge")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // quality: float32, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "quality")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // width: float32, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "width")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // handles: list<struct>, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "handles")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(7)->nullable(), true)
  TEST_EQUAL(s->field(7)->type()->Equals(ConsensusFeatureSchema::handlesType()), true)
  // metavalues: list<struct>, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "metavalues")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(8)->nullable(), true)
  TEST_EQUAL(s->field(8)->type()->Equals(ConsensusFeatureSchema::metavaluesType()), true)
}
END_SECTION

START_SECTION(ConsensusFeatureSchema validation with table)
{
  auto s = ConsensusFeatureSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== PSMSchema ==========

START_SECTION(PSMSchema::schema() returns non-null with 25 fields)
{
  auto s = PSMSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 25)
}
END_SECTION

START_SECTION(PSMSchema column name constants)
{
  TEST_STRING_EQUAL(PSMSchema::SEQUENCE, "sequence")
  TEST_STRING_EQUAL(PSMSchema::PEPTIDOFORM, "peptidoform")
  TEST_STRING_EQUAL(PSMSchema::MODIFICATIONS, "modifications")
  TEST_STRING_EQUAL(PSMSchema::PRECURSOR_CHARGE, "precursor_charge")
  TEST_STRING_EQUAL(PSMSchema::POSTERIOR_ERROR_PROBABILITY, "posterior_error_probability")
  TEST_STRING_EQUAL(PSMSchema::IS_DECOY, "is_decoy")
  TEST_STRING_EQUAL(PSMSchema::CALCULATED_MZ, "calculated_mz")
  TEST_STRING_EQUAL(PSMSchema::OBSERVED_MZ, "observed_mz")
  TEST_STRING_EQUAL(PSMSchema::ADDITIONAL_SCORES, "additional_scores")
  TEST_STRING_EQUAL(PSMSchema::PROTEIN_ACCESSIONS, "protein_accessions")
  TEST_STRING_EQUAL(PSMSchema::PREDICTED_RT, "predicted_rt")
  TEST_STRING_EQUAL(PSMSchema::REFERENCE_FILE_NAME, "reference_file_name")
  TEST_STRING_EQUAL(PSMSchema::CV_PARAMS, "cv_params")
  TEST_STRING_EQUAL(PSMSchema::SCAN, "scan")
  TEST_STRING_EQUAL(PSMSchema::RT, "rt")
  TEST_STRING_EQUAL(PSMSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(PSMSchema::SPECTRUM_REFERENCE, "spectrum_reference")
  TEST_STRING_EQUAL(PSMSchema::SCORE, "score")
  TEST_STRING_EQUAL(PSMSchema::SCORE_TYPE, "score_type")
  TEST_STRING_EQUAL(PSMSchema::HIGHER_SCORE_BETTER, "higher_score_better")
  TEST_STRING_EQUAL(PSMSchema::RANK, "rank")
  TEST_STRING_EQUAL(PSMSchema::PEPTIDE_IDENTIFICATION_INDEX, "peptide_identification_index")
  TEST_STRING_EQUAL(PSMSchema::PSM_METAVALUES, "psm_metavalues")
  TEST_STRING_EQUAL(PSMSchema::SPECTRUM_METAVALUES, "spectrum_metavalues")
  TEST_STRING_EQUAL(PSMSchema::RUN_IDENTIFIER, "run_identifier")
}
END_SECTION

START_SECTION(PSMSchema field types and nullability)
{
  auto s = PSMSchema::schema();
  // sequence: utf8, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "sequence")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // peptidoform: utf8, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "peptidoform")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // modifications: list<struct>, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "modifications")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), true)
  TEST_EQUAL(s->field(2)->type()->Equals(PSMSchema::modificationsType()), true)
  // precursor_charge: int32, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "precursor_charge")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // posterior_error_probability: float64, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // is_decoy: boolean, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "is_decoy")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // calculated_mz: float64, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "calculated_mz")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // observed_mz: float64, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "observed_mz")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // additional_scores: list<struct>, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "additional_scores")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(8)->nullable(), true)
  TEST_EQUAL(s->field(8)->type()->Equals(PSMSchema::additionalScoresType()), true)
  // protein_accessions: list<utf8>, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "protein_accessions")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(9)->nullable(), true)
  // predicted_rt: float64, nullable (default)
  TEST_EQUAL(s->field(10)->name(), "predicted_rt")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // reference_file_name: utf8, nullable (default)
  TEST_EQUAL(s->field(11)->name(), "reference_file_name")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(11)->nullable(), true)
  // cv_params: utf8, nullable (default)
  TEST_EQUAL(s->field(12)->name(), "cv_params")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(12)->nullable(), true)
  // scan: int32, nullable (default)
  TEST_EQUAL(s->field(13)->name(), "scan")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(13)->nullable(), true)
  // rt: float64, nullable (default)
  TEST_EQUAL(s->field(14)->name(), "rt")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(14)->nullable(), true)
  // ion_mobility: float64, nullable (default)
  TEST_EQUAL(s->field(15)->name(), "ion_mobility")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(15)->nullable(), true)
  // spectrum_reference: utf8, nullable (default)
  TEST_EQUAL(s->field(16)->name(), "spectrum_reference")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(16)->nullable(), true)
  // score: float64, nullable (default)
  TEST_EQUAL(s->field(17)->name(), "score")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(17)->nullable(), true)
  // score_type: utf8, nullable (default)
  TEST_EQUAL(s->field(18)->name(), "score_type")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(18)->nullable(), true)
  // higher_score_better: boolean, nullable (default)
  TEST_EQUAL(s->field(19)->name(), "higher_score_better")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(19)->nullable(), true)
  // rank: int32, nullable (default)
  TEST_EQUAL(s->field(20)->name(), "rank")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(20)->nullable(), true)
  // peptide_identification_index: int32, nullable (default)
  TEST_EQUAL(s->field(21)->name(), "peptide_identification_index")
  TEST_EQUAL(s->field(21)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(21)->nullable(), true)
  // psm_metavalues: list<struct>, nullable (default)
  TEST_EQUAL(s->field(22)->name(), "psm_metavalues")
  TEST_EQUAL(s->field(22)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(22)->nullable(), true)
  TEST_EQUAL(s->field(22)->type()->Equals(PSMSchema::metavaluesType()), true)
  // spectrum_metavalues: list<struct>, nullable (default)
  TEST_EQUAL(s->field(23)->name(), "spectrum_metavalues")
  TEST_EQUAL(s->field(23)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(23)->nullable(), true)
  TEST_EQUAL(s->field(23)->type()->Equals(PSMSchema::metavaluesType()), true)
  // run_identifier: utf8, nullable (default)
  TEST_EQUAL(s->field(24)->name(), "run_identifier")
  TEST_EQUAL(s->field(24)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(24)->nullable(), true)
}
END_SECTION

START_SECTION(PSMSchema validation with table)
{
  auto s = PSMSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== ConsensusFeatureExportSchema ==========

START_SECTION(ConsensusFeatureExportSchema::schema() returns non-null with 33 fields)
{
  auto s = ConsensusFeatureExportSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 33)
}
END_SECTION

START_SECTION(ConsensusFeatureExportSchema column name constants)
{
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SEQUENCE, "sequence")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::PEPTIDOFORM, "peptidoform")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::MODIFICATIONS, "modifications")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::PRECURSOR_CHARGE, "precursor_charge")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::CALCULATED_MZ, "calculated_mz")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::OBSERVED_MZ, "observed_mz")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::RT, "rt")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::POSTERIOR_ERROR_PROBABILITY, "posterior_error_probability")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::IS_DECOY, "is_decoy")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::ADDITIONAL_SCORES, "additional_scores")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::PREDICTED_RT, "predicted_rt")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::REFERENCE_FILE_NAME, "reference_file_name")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::CV_PARAMS, "cv_params")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SCAN, "scan")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::START_ION_MOBILITY, "start_ion_mobility")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::STOP_ION_MOBILITY, "stop_ion_mobility")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::INTENSITIES, "intensities")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::ADDITIONAL_INTENSITIES, "additional_intensities")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::PG_ACCESSIONS, "pg_accessions")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::ANCHOR_PROTEIN, "anchor_protein")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::UNIQUE, "unique")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::PG_GLOBAL_QVALUE, "pg_global_qvalue")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::GG_ACCESSIONS, "gg_accessions")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::GG_NAMES, "gg_names")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SCAN_REFERENCE_FILE_NAME, "scan_reference_file_name")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::RT_START, "rt_start")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::RT_STOP, "rt_stop")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::QUALITY, "quality")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SCORE, "score")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SCORE_TYPE, "score_type")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::SPECTRUM_REFERENCE, "spectrum_reference")
  TEST_STRING_EQUAL(ConsensusFeatureExportSchema::FEATURE_METAVALUES, "feature_metavalues")
}
END_SECTION

START_SECTION(ConsensusFeatureExportSchema field types and nullability)
{
  auto s = ConsensusFeatureExportSchema::schema();
  // sequence: utf8, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "sequence")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // peptidoform: utf8, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "peptidoform")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // modifications: list<struct>, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "modifications")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), true)
  TEST_EQUAL(s->field(2)->type()->Equals(ConsensusFeatureExportSchema::modificationsType()), true)
  // precursor_charge: int32, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "precursor_charge")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // calculated_mz: float32, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "calculated_mz")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // observed_mz: float32, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "observed_mz")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // rt: float32, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "rt")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // posterior_error_probability: float64, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "posterior_error_probability")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // is_decoy: int32, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "is_decoy")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // additional_scores: list<struct>, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "additional_scores")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(9)->nullable(), true)
  TEST_EQUAL(s->field(9)->type()->Equals(ConsensusFeatureExportSchema::additionalScoresType()), true)
  // predicted_rt: float32, nullable (default)
  TEST_EQUAL(s->field(10)->name(), "predicted_rt")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // reference_file_name: utf8, nullable (default)
  TEST_EQUAL(s->field(11)->name(), "reference_file_name")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(11)->nullable(), true)
  // cv_params: utf8, nullable (default)
  TEST_EQUAL(s->field(12)->name(), "cv_params")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(12)->nullable(), true)
  // scan: utf8, nullable (default)
  TEST_EQUAL(s->field(13)->name(), "scan")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(13)->nullable(), true)
  // ion_mobility: float32, nullable (default)
  TEST_EQUAL(s->field(14)->name(), "ion_mobility")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(14)->nullable(), true)
  // start_ion_mobility: float32, nullable (default)
  TEST_EQUAL(s->field(15)->name(), "start_ion_mobility")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(15)->nullable(), true)
  // stop_ion_mobility: float32, nullable (default)
  TEST_EQUAL(s->field(16)->name(), "stop_ion_mobility")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(16)->nullable(), true)
  // intensities: list<struct>, nullable (default)
  TEST_EQUAL(s->field(17)->name(), "intensities")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(17)->nullable(), true)
  TEST_EQUAL(s->field(17)->type()->Equals(ConsensusFeatureExportSchema::intensitiesType()), true)
  // additional_intensities: utf8, nullable (default)
  TEST_EQUAL(s->field(18)->name(), "additional_intensities")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(18)->nullable(), true)
  // pg_accessions: list<utf8>, nullable (default)
  TEST_EQUAL(s->field(19)->name(), "pg_accessions")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(19)->nullable(), true)
  // anchor_protein: utf8, nullable (default)
  TEST_EQUAL(s->field(20)->name(), "anchor_protein")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(20)->nullable(), true)
  // unique: int32, nullable (default)
  TEST_EQUAL(s->field(21)->name(), "unique")
  TEST_EQUAL(s->field(21)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(21)->nullable(), true)
  // pg_global_qvalue: float64, nullable (default)
  TEST_EQUAL(s->field(22)->name(), "pg_global_qvalue")
  TEST_EQUAL(s->field(22)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(22)->nullable(), true)
  // gg_accessions: list<utf8>, nullable (default)
  TEST_EQUAL(s->field(23)->name(), "gg_accessions")
  TEST_EQUAL(s->field(23)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(23)->nullable(), true)
  // gg_names: list<utf8>, nullable (default)
  TEST_EQUAL(s->field(24)->name(), "gg_names")
  TEST_EQUAL(s->field(24)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(24)->nullable(), true)
  // scan_reference_file_name: utf8, nullable (default)
  TEST_EQUAL(s->field(25)->name(), "scan_reference_file_name")
  TEST_EQUAL(s->field(25)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(25)->nullable(), true)
  // rt_start: float32, nullable (default)
  TEST_EQUAL(s->field(26)->name(), "rt_start")
  TEST_EQUAL(s->field(26)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(26)->nullable(), true)
  // rt_stop: float32, nullable (default)
  TEST_EQUAL(s->field(27)->name(), "rt_stop")
  TEST_EQUAL(s->field(27)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(27)->nullable(), true)
  // quality: float32, nullable (default)
  TEST_EQUAL(s->field(28)->name(), "quality")
  TEST_EQUAL(s->field(28)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(28)->nullable(), true)
  // score: float64, nullable (default)
  TEST_EQUAL(s->field(29)->name(), "score")
  TEST_EQUAL(s->field(29)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(29)->nullable(), true)
  // score_type: utf8, nullable (default)
  TEST_EQUAL(s->field(30)->name(), "score_type")
  TEST_EQUAL(s->field(30)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(30)->nullable(), true)
  // spectrum_reference: utf8, nullable (default)
  TEST_EQUAL(s->field(31)->name(), "spectrum_reference")
  TEST_EQUAL(s->field(31)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(31)->nullable(), true)
  // feature_metavalues: list<struct>, nullable (default)
  TEST_EQUAL(s->field(32)->name(), "feature_metavalues")
  TEST_EQUAL(s->field(32)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(32)->nullable(), true)
  TEST_EQUAL(s->field(32)->type()->Equals(ConsensusFeatureExportSchema::metavaluesType()), true)
}
END_SECTION

START_SECTION(ConsensusFeatureExportSchema nested types match PSMSchema)
{
  // modificationsType and additionalScoresType should be the same as PSMSchema's
  TEST_EQUAL(ConsensusFeatureExportSchema::modificationsType()->Equals(PSMSchema::modificationsType()), true)
  TEST_EQUAL(ConsensusFeatureExportSchema::additionalScoresType()->Equals(PSMSchema::additionalScoresType()), true)
}
END_SECTION

START_SECTION(ConsensusFeatureExportSchema validation with table)
{
  auto s = ConsensusFeatureExportSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== SpectraLongSchema ==========

START_SECTION(SpectraLongSchema::schema() returns non-null with 12 fields)
{
  auto s = SpectraLongSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 12)
}
END_SECTION

START_SECTION(SpectraLongSchema column name constants)
{
  TEST_STRING_EQUAL(SpectraLongSchema::MZ, "mz")
  TEST_STRING_EQUAL(SpectraLongSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(SpectraLongSchema::RT, "rt")
  TEST_STRING_EQUAL(SpectraLongSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(SpectraLongSchema::SPECTRUM_INDEX, "spectrum_index")
  TEST_STRING_EQUAL(SpectraLongSchema::MS_LEVEL, "ms_level")
  TEST_STRING_EQUAL(SpectraLongSchema::NATIVE_ID, "native_id")
  TEST_STRING_EQUAL(SpectraLongSchema::PRECURSOR_MZ, "precursor_mz")
  TEST_STRING_EQUAL(SpectraLongSchema::PRECURSOR_CHARGE, "precursor_charge")
  TEST_STRING_EQUAL(SpectraLongSchema::PRECURSOR_INTENSITY, "precursor_intensity")
  TEST_STRING_EQUAL(SpectraLongSchema::ISOLATION_LOWER, "isolation_lower")
  TEST_STRING_EQUAL(SpectraLongSchema::ISOLATION_UPPER, "isolation_upper")
}
END_SECTION

START_SECTION(SpectraLongSchema field types and nullability)
{
  auto s = SpectraLongSchema::schema();
  // mz: float64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "mz")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // intensity: float32, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "intensity")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // rt: float32, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "rt")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // ion_mobility: float32, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "ion_mobility")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // spectrum_index: uint32, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "spectrum_index")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::UINT32)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // ms_level: uint8, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "ms_level")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::UINT8)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // native_id: utf8, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "native_id")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // precursor_mz: float64, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "precursor_mz")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // precursor_charge: int16, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "precursor_charge")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // precursor_intensity: float32, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "precursor_intensity")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(9)->nullable(), true)
  // isolation_lower: float64, nullable (default)
  TEST_EQUAL(s->field(10)->name(), "isolation_lower")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // isolation_upper: float64, nullable (default)
  TEST_EQUAL(s->field(11)->name(), "isolation_upper")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(11)->nullable(), true)
}
END_SECTION

START_SECTION(SpectraLongSchema subset validation with table)
{
  auto s = SpectraLongSchema::schema();
  // Build a table with only a subset of fields (mz, rt, intensity)
  auto subset_schema = arrow::schema({
    arrow::field("mz", arrow::float64()),
    arrow::field("intensity", arrow::float32()),
    arrow::field("rt", arrow::float32()),
  });
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < subset_schema->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(subset_schema->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(subset_schema, columns);
  auto result = ArrowSchemaValidation::validate(table, s,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(SpectraLongSchema validation with full table)
{
  auto s = SpectraLongSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== SpectraSemiWideSchema ==========

START_SECTION(SpectraSemiWideSchema::schema() returns non-null with 12 fields)
{
  auto s = SpectraSemiWideSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 12)
}
END_SECTION

START_SECTION(SpectraSemiWideSchema column name constants)
{
  TEST_STRING_EQUAL(SpectraSemiWideSchema::SPECTRUM_INDEX, "spectrum_index")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::RT, "rt")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::MS_LEVEL, "ms_level")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::NATIVE_ID, "native_id")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::MZ, "mz")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::PRECURSOR_MZ, "precursor_mz")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::PRECURSOR_CHARGE, "precursor_charge")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::PRECURSOR_INTENSITY, "precursor_intensity")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::ISOLATION_LOWER, "isolation_lower")
  TEST_STRING_EQUAL(SpectraSemiWideSchema::ISOLATION_UPPER, "isolation_upper")
}
END_SECTION

START_SECTION(SpectraSemiWideSchema field types and nullability)
{
  auto s = SpectraSemiWideSchema::schema();
  // spectrum_index: uint32, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "spectrum_index")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::UINT32)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // rt: float32, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "rt")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // ms_level: uint8, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "ms_level")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::UINT8)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // native_id: utf8, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "native_id")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // mz: list<float64>, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "mz")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(4)->nullable(), true)
  TEST_EQUAL(s->field(4)->type()->Equals(arrow::list(arrow::float64())), true)
  // intensity: list<float32>, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "intensity")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(5)->nullable(), true)
  TEST_EQUAL(s->field(5)->type()->Equals(arrow::list(arrow::float32())), true)
  // ion_mobility: list<float32>, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "ion_mobility")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(6)->nullable(), true)
  TEST_EQUAL(s->field(6)->type()->Equals(arrow::list(arrow::float32())), true)
  // precursor_mz: float64, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "precursor_mz")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // precursor_charge: int16, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "precursor_charge")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // precursor_intensity: float32, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "precursor_intensity")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(9)->nullable(), true)
  // isolation_lower: float64, nullable (default)
  TEST_EQUAL(s->field(10)->name(), "isolation_lower")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // isolation_upper: float64, nullable (default)
  TEST_EQUAL(s->field(11)->name(), "isolation_upper")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(11)->nullable(), true)
}
END_SECTION

START_SECTION(SpectraSemiWideSchema differs from SpectraLongSchema in mz/intensity/ion_mobility types)
{
  auto long_s = SpectraLongSchema::schema();
  auto semi_s = SpectraSemiWideSchema::schema();
  // Long format: mz is scalar float64
  TEST_EQUAL(long_s->GetFieldByName("mz")->type()->id(), arrow::Type::DOUBLE)
  // Semi-wide format: mz is list<float64>
  TEST_EQUAL(semi_s->GetFieldByName("mz")->type()->id(), arrow::Type::LIST)
  // Long format: intensity is scalar float32
  TEST_EQUAL(long_s->GetFieldByName("intensity")->type()->id(), arrow::Type::FLOAT)
  // Semi-wide format: intensity is list<float32>
  TEST_EQUAL(semi_s->GetFieldByName("intensity")->type()->id(), arrow::Type::LIST)
  // Long format: ion_mobility is scalar float32
  TEST_EQUAL(long_s->GetFieldByName("ion_mobility")->type()->id(), arrow::Type::FLOAT)
  // Semi-wide format: ion_mobility is list<float32>
  TEST_EQUAL(semi_s->GetFieldByName("ion_mobility")->type()->id(), arrow::Type::LIST)
}
END_SECTION

START_SECTION(SpectraSemiWideSchema subset validation with table)
{
  auto s = SpectraSemiWideSchema::schema();
  // Build a table with only a subset of fields
  auto subset_schema = arrow::schema({
    arrow::field("spectrum_index", arrow::uint32()),
    arrow::field("mz", arrow::list(arrow::float64())),
    arrow::field("intensity", arrow::list(arrow::float32())),
  });
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < subset_schema->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(subset_schema->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(subset_schema, columns);
  auto result = ArrowSchemaValidation::validate(table, s,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(SpectraSemiWideSchema validation with full table)
{
  auto s = SpectraSemiWideSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== ChromatogramSchema ==========

START_SECTION(ChromatogramSchema::schema() returns non-null with 6 fields)
{
  auto s = ChromatogramSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 6)
}
END_SECTION

START_SECTION(ChromatogramSchema column name constants)
{
  TEST_STRING_EQUAL(ChromatogramSchema::RT, "rt")
  TEST_STRING_EQUAL(ChromatogramSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(ChromatogramSchema::CHROMATOGRAM_INDEX, "chromatogram_index")
  TEST_STRING_EQUAL(ChromatogramSchema::NATIVE_ID, "native_id")
  TEST_STRING_EQUAL(ChromatogramSchema::PRECURSOR_MZ, "precursor_mz")
  TEST_STRING_EQUAL(ChromatogramSchema::PRODUCT_MZ, "product_mz")
}
END_SECTION

START_SECTION(ChromatogramSchema field types and nullability)
{
  auto s = ChromatogramSchema::schema();
  // rt: float64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "rt")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // intensity: float32, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "intensity")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // chromatogram_index: uint32, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "chromatogram_index")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::UINT32)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // native_id: utf8, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "native_id")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // precursor_mz: float64, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "precursor_mz")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // product_mz: float64, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "product_mz")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(5)->nullable(), true)
}
END_SECTION

START_SECTION(ChromatogramSchema subset validation with table)
{
  auto s = ChromatogramSchema::schema();
  // Build a table with only a subset of fields (rt, intensity)
  auto subset_schema = arrow::schema({
    arrow::field("rt", arrow::float64()),
    arrow::field("intensity", arrow::float32()),
  });
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < subset_schema->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(subset_schema->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(subset_schema, columns);
  auto result = ArrowSchemaValidation::validate(table, s,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(ChromatogramSchema validation with full table)
{
  auto s = ChromatogramSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== OSWPrecursorSchema ==========

START_SECTION(OSWPrecursorSchema::schema() returns non-null with 10 fields)
{
  auto s = OSWPrecursorSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 10)
}
END_SECTION

START_SECTION(OSWPrecursorSchema column name constants)
{
  TEST_STRING_EQUAL(OSWPrecursorSchema::PRECURSOR_ID, "precursor_id")
  TEST_STRING_EQUAL(OSWPrecursorSchema::PRECURSOR_MZ, "precursor_mz")
  TEST_STRING_EQUAL(OSWPrecursorSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(OSWPrecursorSchema::LIBRARY_RT, "library_rt")
  TEST_STRING_EQUAL(OSWPrecursorSchema::LIBRARY_DRIFT_TIME, "library_drift_time")
  TEST_STRING_EQUAL(OSWPrecursorSchema::DECOY, "decoy")
  TEST_STRING_EQUAL(OSWPrecursorSchema::TRAML_ID, "traml_id")
  TEST_STRING_EQUAL(OSWPrecursorSchema::MODIFIED_SEQUENCE, "modified_sequence")
  TEST_STRING_EQUAL(OSWPrecursorSchema::UNMODIFIED_SEQUENCE, "unmodified_sequence")
  TEST_STRING_EQUAL(OSWPrecursorSchema::PROTEIN_ACCESSIONS, "protein_accessions")
}
END_SECTION

START_SECTION(OSWPrecursorSchema field types and nullability)
{
  auto s = OSWPrecursorSchema::schema();
  // precursor_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "precursor_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // precursor_mz: float64, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "precursor_mz")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // charge: int32, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "charge")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // library_rt: float64, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "library_rt")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // library_drift_time: float64, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "library_drift_time")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // decoy: boolean, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "decoy")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // traml_id: utf8, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "traml_id")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // modified_sequence: utf8, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "modified_sequence")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // unmodified_sequence: utf8, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "unmodified_sequence")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // protein_accessions: utf8, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "protein_accessions")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(9)->nullable(), true)
}
END_SECTION

START_SECTION(OSWPrecursorSchema validation with table)
{
  auto s = OSWPrecursorSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== OSWTransitionSchema ==========

START_SECTION(OSWTransitionSchema::schema() returns non-null with 13 fields)
{
  auto s = OSWTransitionSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 13)
}
END_SECTION

START_SECTION(OSWTransitionSchema column name constants)
{
  TEST_STRING_EQUAL(OSWTransitionSchema::TRANSITION_ID, "transition_id")
  TEST_STRING_EQUAL(OSWTransitionSchema::PRECURSOR_ID, "precursor_id")
  TEST_STRING_EQUAL(OSWTransitionSchema::TRAML_ID, "traml_id")
  TEST_STRING_EQUAL(OSWTransitionSchema::PRODUCT_MZ, "product_mz")
  TEST_STRING_EQUAL(OSWTransitionSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(OSWTransitionSchema::TYPE, "type")
  TEST_STRING_EQUAL(OSWTransitionSchema::ANNOTATION, "annotation")
  TEST_STRING_EQUAL(OSWTransitionSchema::ORDINAL, "ordinal")
  TEST_STRING_EQUAL(OSWTransitionSchema::DETECTING, "detecting")
  TEST_STRING_EQUAL(OSWTransitionSchema::IDENTIFYING, "identifying")
  TEST_STRING_EQUAL(OSWTransitionSchema::QUANTIFYING, "quantifying")
  TEST_STRING_EQUAL(OSWTransitionSchema::LIBRARY_INTENSITY, "library_intensity")
  TEST_STRING_EQUAL(OSWTransitionSchema::DECOY, "decoy")
}
END_SECTION

START_SECTION(OSWTransitionSchema field types and nullability)
{
  auto s = OSWTransitionSchema::schema();
  // transition_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "transition_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // precursor_id: int64, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "precursor_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // traml_id: utf8, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "traml_id")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // product_mz: float64, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "product_mz")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // charge: int32, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "charge")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // type: utf8, nullable (default)
  TEST_EQUAL(s->field(5)->name(), "type")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(5)->nullable(), true)
  // annotation: utf8, nullable (default)
  TEST_EQUAL(s->field(6)->name(), "annotation")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(6)->nullable(), true)
  // ordinal: int32, nullable (default)
  TEST_EQUAL(s->field(7)->name(), "ordinal")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(7)->nullable(), true)
  // detecting: boolean, nullable (default)
  TEST_EQUAL(s->field(8)->name(), "detecting")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // identifying: boolean, nullable (default)
  TEST_EQUAL(s->field(9)->name(), "identifying")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(9)->nullable(), true)
  // quantifying: boolean, nullable (default)
  TEST_EQUAL(s->field(10)->name(), "quantifying")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // library_intensity: float64, nullable (default)
  TEST_EQUAL(s->field(11)->name(), "library_intensity")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(11)->nullable(), true)
  // decoy: boolean, nullable (default)
  TEST_EQUAL(s->field(12)->name(), "decoy")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(12)->nullable(), true)
}
END_SECTION

START_SECTION(OSWTransitionSchema validation with table)
{
  auto s = OSWTransitionSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== OSWFeaturePrecursorSchema ==========

START_SECTION(OSWFeaturePrecursorSchema::schema() returns non-null with 5 fields)
{
  auto s = OSWFeaturePrecursorSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 5)
}
END_SECTION

START_SECTION(OSWFeaturePrecursorSchema column name constants)
{
  TEST_STRING_EQUAL(OSWFeaturePrecursorSchema::FEATURE_ID, "feature_id")
  TEST_STRING_EQUAL(OSWFeaturePrecursorSchema::RUN_ID, "run_id")
  TEST_STRING_EQUAL(OSWFeaturePrecursorSchema::PRECURSOR_ISOTOPE, "precursor_isotope")
  TEST_STRING_EQUAL(OSWFeaturePrecursorSchema::PRECURSOR_AREA_INTENSITY, "precursor_area_intensity")
  TEST_STRING_EQUAL(OSWFeaturePrecursorSchema::PRECURSOR_APEX_INTENSITY, "precursor_apex_intensity")
}
END_SECTION

START_SECTION(OSWFeaturePrecursorSchema field types and nullability)
{
  auto s = OSWFeaturePrecursorSchema::schema();
  // feature_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "feature_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // run_id: int64, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "run_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // precursor_isotope: int32, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "precursor_isotope")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT32)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // precursor_area_intensity: float64, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "precursor_area_intensity")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // precursor_apex_intensity: float64, nullable (default)
  TEST_EQUAL(s->field(4)->name(), "precursor_apex_intensity")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
}
END_SECTION

START_SECTION(OSWFeaturePrecursorSchema validation with table)
{
  auto s = OSWFeaturePrecursorSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

// ========== OSWRunSchema ==========

START_SECTION(OSWRunSchema::schema() returns non-null with 2 fields)
{
  auto s = OSWRunSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 2)
}
END_SECTION

START_SECTION(OSWRunSchema column name constants)
{
  TEST_STRING_EQUAL(OSWRunSchema::RUN_ID, "run_id")
  TEST_STRING_EQUAL(OSWRunSchema::FILENAME, "filename")
}
END_SECTION

START_SECTION(OSWRunSchema field types and nullability)
{
  auto s = OSWRunSchema::schema();
  // run_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "run_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // filename: utf8, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "filename")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), true)
}
END_SECTION

START_SECTION(OSWRunSchema validation with table)
{
  auto s = OSWRunSchema::schema();
  std::vector<std::shared_ptr<arrow::Array>> columns;
  for (int i = 0; i < s->num_fields(); ++i)
  {
    auto result = arrow::MakeEmptyArray(s->field(i)->type());
    columns.push_back(result.ValueOrDie());
  }
  auto table = arrow::Table::Make(s, columns);
  auto result = ArrowSchemaValidation::validate(table, s);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

END_TEST

#else // WITH_PARQUET

START_TEST(ArrowSchemaRegistry, "$Id$")
END_TEST

#endif
