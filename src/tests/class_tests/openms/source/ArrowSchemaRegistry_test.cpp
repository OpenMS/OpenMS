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

END_TEST

#else // WITH_PARQUET

START_TEST(ArrowSchemaRegistry, "$Id$")
END_TEST

#endif
