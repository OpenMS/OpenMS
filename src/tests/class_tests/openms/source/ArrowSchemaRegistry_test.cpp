#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>
using namespace OpenMS;

// Helper functions to build small test arrays without arrow::json::ArrayFromJSONString
// (which is not available in static Arrow builds on Windows).
namespace {
  std::shared_ptr<arrow::Array> makeUtf8Array(const std::vector<std::string>& values)
  {
    arrow::StringBuilder builder;
    for (const auto& v : values) { (void)builder.Append(v); }
    return builder.Finish().ValueOrDie();
  }

  std::shared_ptr<arrow::Array> makeFloat64Array(const std::vector<double>& values)
  {
    arrow::DoubleBuilder builder;
    for (auto v : values) { (void)builder.Append(v); }
    return builder.Finish().ValueOrDie();
  }

  std::shared_ptr<arrow::Array> makeInt32Array(const std::vector<int32_t>& values)
  {
    arrow::Int32Builder builder;
    for (auto v : values) { (void)builder.Append(v); }
    return builder.Finish().ValueOrDie();
  }

  std::shared_ptr<arrow::Array> makeTimestampArray(
    const std::shared_ptr<arrow::DataType>& type,
    const std::vector<int64_t>& values)
  {
    arrow::TimestampBuilder builder(type, arrow::default_memory_pool());
    for (auto v : values) { (void)builder.Append(v); }
    return builder.Finish().ValueOrDie();
  }
} // anonymous namespace

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
  auto arr_a = makeUtf8Array({"x", "y"});
  auto arr_b = makeFloat64Array({1.0, 2.0});
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
  auto arr_a = makeUtf8Array({"x"});
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
  auto arr_a = makeUtf8Array({"x"});
  auto arr_b = makeFloat64Array({1.0});
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
  auto arr_a = makeInt32Array({1});
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
  auto arr_a = makeUtf8Array({"x"});
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
  auto arr_a = makeUtf8Array({"x"});
  auto arr_c = makeInt32Array({1});
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_c});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
}
END_SECTION

START_SECTION(validate - Subset mode - unknown field ignored for forward compatibility)
{
  auto expected = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
  });
  auto actual_schema = arrow::schema({
    arrow::field("a", arrow::utf8(), false),
    arrow::field("z", arrow::float64(), true),
  });
  auto arr_a = makeUtf8Array({"x"});
  auto arr_z = makeFloat64Array({1.0});
  auto table = arrow::Table::Make(actual_schema, {arr_a, arr_z});

  auto result = ArrowSchemaValidation::validate(table, expected,
    ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
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
  auto arr_a = makeInt32Array({1});
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
  auto arr_a = makeUtf8Array({"x"});
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
  TEST_EQUAL(s->num_fields(), 29)
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

// ========== QPXPSMSchema ==========

START_SECTION(QPXPSMSchema::schema() returns non-null with 24 fields)
{
  auto s = QPXPSMSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 24)
}
END_SECTION

START_SECTION(QPXPSMSchema column name constants)
{
  TEST_STRING_EQUAL(QPXPSMSchema::SEQUENCE, "sequence")
  TEST_STRING_EQUAL(QPXPSMSchema::PEPTIDOFORM, "peptidoform")
  TEST_STRING_EQUAL(QPXPSMSchema::MODIFICATIONS, "modifications")
  TEST_STRING_EQUAL(QPXPSMSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(QPXPSMSchema::POSTERIOR_ERROR_PROBABILITY, "posterior_error_probability")
  TEST_STRING_EQUAL(QPXPSMSchema::IS_DECOY, "is_decoy")
  TEST_STRING_EQUAL(QPXPSMSchema::CALCULATED_MZ, "calculated_mz")
  TEST_STRING_EQUAL(QPXPSMSchema::OBSERVED_MZ, "observed_mz")
  TEST_STRING_EQUAL(QPXPSMSchema::MASS_ERROR_PPM, "mass_error_ppm")
  TEST_STRING_EQUAL(QPXPSMSchema::ADDITIONAL_SCORES, "additional_scores")
  TEST_STRING_EQUAL(QPXPSMSchema::PREDICTED_RT, "predicted_rt")
  TEST_STRING_EQUAL(QPXPSMSchema::RUN_FILE_NAME, "run_file_name")
  TEST_STRING_EQUAL(QPXPSMSchema::CV_PARAMS, "cv_params")
  TEST_STRING_EQUAL(QPXPSMSchema::SCAN, "scan")
  TEST_STRING_EQUAL(QPXPSMSchema::RT, "rt")
  TEST_STRING_EQUAL(QPXPSMSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(QPXPSMSchema::MISSED_CLEAVAGES, "missed_cleavages")
  TEST_STRING_EQUAL(QPXPSMSchema::PROTEIN_ACCESSIONS, "protein_accessions")
  TEST_STRING_EQUAL(QPXPSMSchema::CROSS_LINKS, "cross_links")
  TEST_STRING_EQUAL(QPXPSMSchema::MZ_ARRAY, "mz_array")
  TEST_STRING_EQUAL(QPXPSMSchema::INTENSITY_ARRAY, "intensity_array")
  TEST_STRING_EQUAL(QPXPSMSchema::CHARGE_ARRAY, "charge_array")
  TEST_STRING_EQUAL(QPXPSMSchema::ION_TYPE_ARRAY, "ion_type_array")
  TEST_STRING_EQUAL(QPXPSMSchema::ION_MOBILITY_ARRAY, "ion_mobility_array")
}
END_SECTION

START_SECTION(QPXPSMSchema field types and nullability)
{
  auto s = QPXPSMSchema::schema();
  // sequence: utf8, not null
  TEST_EQUAL(s->field(0)->name(), "sequence")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // peptidoform: utf8, not null
  TEST_EQUAL(s->field(1)->name(), "peptidoform")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), false)
  // modifications: list<struct>, nullable
  TEST_EQUAL(s->field(2)->name(), "modifications")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), true)
  TEST_EQUAL(s->field(2)->type()->Equals(QPXPSMSchema::modificationsType()), true)
  // charge: int16, not null
  TEST_EQUAL(s->field(3)->name(), "charge")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(3)->nullable(), false)
  // posterior_error_probability: float64, nullable
  TEST_EQUAL(s->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // is_decoy: bool, not null
  TEST_EQUAL(s->field(5)->name(), "is_decoy")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(5)->nullable(), false)
  // calculated_mz: float32, not null
  TEST_EQUAL(s->field(6)->name(), "calculated_mz")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(6)->nullable(), false)
  // observed_mz: float32, not null
  TEST_EQUAL(s->field(7)->name(), "observed_mz")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(7)->nullable(), false)
  // mass_error_ppm: float32, nullable
  TEST_EQUAL(s->field(8)->name(), "mass_error_ppm")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // additional_scores: list<struct>, nullable
  TEST_EQUAL(s->field(9)->name(), "additional_scores")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(9)->nullable(), true)
  TEST_EQUAL(s->field(9)->type()->Equals(QPXPSMSchema::additionalScoresType()), true)
  // predicted_rt: float32, nullable
  TEST_EQUAL(s->field(10)->name(), "predicted_rt")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // run_file_name: utf8, not null
  TEST_EQUAL(s->field(11)->name(), "run_file_name")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(11)->nullable(), false)
  // cv_params: list<struct>, nullable
  TEST_EQUAL(s->field(12)->name(), "cv_params")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(12)->nullable(), true)
  TEST_EQUAL(s->field(12)->type()->Equals(QPXPSMSchema::cvParamsType()), true)
  // scan: list<int32>, not null
  TEST_EQUAL(s->field(13)->name(), "scan")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(13)->nullable(), false)
  TEST_EQUAL(std::static_pointer_cast<arrow::ListType>(s->field(13)->type())->value_type()->id(), arrow::Type::INT32)
  // rt: float32, nullable
  TEST_EQUAL(s->field(14)->name(), "rt")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(14)->nullable(), true)
  // ion_mobility: float32, nullable
  TEST_EQUAL(s->field(15)->name(), "ion_mobility")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(15)->nullable(), true)
  // missed_cleavages: int16, nullable
  TEST_EQUAL(s->field(16)->name(), "missed_cleavages")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(16)->nullable(), true)
  // protein_accessions: list<utf8>, nullable
  TEST_EQUAL(s->field(17)->name(), "protein_accessions")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(17)->nullable(), true)
  // cross_links: list<struct>, nullable
  TEST_EQUAL(s->field(18)->name(), "cross_links")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(18)->nullable(), true)
  TEST_EQUAL(s->field(18)->type()->Equals(QPXPSMSchema::crossLinksType()), true)
  // mz_array: list<float>, nullable
  TEST_EQUAL(s->field(19)->name(), "mz_array")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(19)->nullable(), true)
  // intensity_array: list<float>, nullable
  TEST_EQUAL(s->field(20)->name(), "intensity_array")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(20)->nullable(), true)
  // charge_array: list<int32>, nullable
  TEST_EQUAL(s->field(21)->name(), "charge_array")
  TEST_EQUAL(s->field(21)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(21)->nullable(), true)
  // ion_type_array: list<string>, nullable
  TEST_EQUAL(s->field(22)->name(), "ion_type_array")
  TEST_EQUAL(s->field(22)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(22)->nullable(), true)
  // ion_mobility_array: list<float>, nullable
  TEST_EQUAL(s->field(23)->name(), "ion_mobility_array")
  TEST_EQUAL(s->field(23)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(23)->nullable(), true)
}
END_SECTION

START_SECTION(QPXPSMSchema validation with table)
{
  auto s = QPXPSMSchema::schema();
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

// ========== QPXFeatureSchema ==========

START_SECTION(QPXFeatureSchema::schema() returns non-null with 31 fields)
{
  auto s = QPXFeatureSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 31)
}
END_SECTION

START_SECTION(QPXFeatureSchema column name constants)
{
  TEST_STRING_EQUAL(QPXFeatureSchema::SEQUENCE, "sequence")
  TEST_STRING_EQUAL(QPXFeatureSchema::PEPTIDOFORM, "peptidoform")
  TEST_STRING_EQUAL(QPXFeatureSchema::MODIFICATIONS, "modifications")
  TEST_STRING_EQUAL(QPXFeatureSchema::CHARGE, "charge")
  TEST_STRING_EQUAL(QPXFeatureSchema::POSTERIOR_ERROR_PROBABILITY, "posterior_error_probability")
  TEST_STRING_EQUAL(QPXFeatureSchema::IS_DECOY, "is_decoy")
  TEST_STRING_EQUAL(QPXFeatureSchema::CALCULATED_MZ, "calculated_mz")
  TEST_STRING_EQUAL(QPXFeatureSchema::OBSERVED_MZ, "observed_mz")
  TEST_STRING_EQUAL(QPXFeatureSchema::MASS_ERROR_PPM, "mass_error_ppm")
  TEST_STRING_EQUAL(QPXFeatureSchema::ADDITIONAL_SCORES, "additional_scores")
  TEST_STRING_EQUAL(QPXFeatureSchema::PREDICTED_RT, "predicted_rt")
  TEST_STRING_EQUAL(QPXFeatureSchema::RUN_FILE_NAME, "run_file_name")
  TEST_STRING_EQUAL(QPXFeatureSchema::CV_PARAMS, "cv_params")
  TEST_STRING_EQUAL(QPXFeatureSchema::SCAN, "scan")
  TEST_STRING_EQUAL(QPXFeatureSchema::RT, "rt")
  TEST_STRING_EQUAL(QPXFeatureSchema::ION_MOBILITY, "ion_mobility")
  TEST_STRING_EQUAL(QPXFeatureSchema::MISSED_CLEAVAGES, "missed_cleavages")
  TEST_STRING_EQUAL(QPXFeatureSchema::INTENSITIES, "intensities")
  TEST_STRING_EQUAL(QPXFeatureSchema::ADDITIONAL_INTENSITIES, "additional_intensities")
  TEST_STRING_EQUAL(QPXFeatureSchema::PG_ACCESSIONS, "pg_accessions")
  TEST_STRING_EQUAL(QPXFeatureSchema::ANCHOR_PROTEIN, "anchor_protein")
  TEST_STRING_EQUAL(QPXFeatureSchema::UNIQUE, "unique")
  TEST_STRING_EQUAL(QPXFeatureSchema::PG_GLOBAL_QVALUE, "pg_global_qvalue")
  TEST_STRING_EQUAL(QPXFeatureSchema::PG_POSITIONS, "pg_positions")
  TEST_STRING_EQUAL(QPXFeatureSchema::ION_MOBILITY_START, "ion_mobility_start")
  TEST_STRING_EQUAL(QPXFeatureSchema::ION_MOBILITY_STOP, "ion_mobility_stop")
  TEST_STRING_EQUAL(QPXFeatureSchema::GG_ACCESSIONS, "gg_accessions")
  TEST_STRING_EQUAL(QPXFeatureSchema::GG_NAMES, "gg_names")
  TEST_STRING_EQUAL(QPXFeatureSchema::ID_RUN_FILE_NAME, "id_run_file_name")
  TEST_STRING_EQUAL(QPXFeatureSchema::RT_START, "rt_start")
  TEST_STRING_EQUAL(QPXFeatureSchema::RT_STOP, "rt_stop")
}
END_SECTION

START_SECTION(QPXFeatureSchema field types and nullability)
{
  auto s = QPXFeatureSchema::schema();
  // sequence: utf8, not null
  TEST_EQUAL(s->field(0)->name(), "sequence")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(0)->nullable(), false)
  // peptidoform: utf8, not null
  TEST_EQUAL(s->field(1)->name(), "peptidoform")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), false)
  // modifications: list<struct>, nullable
  TEST_EQUAL(s->field(2)->name(), "modifications")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), true)
  TEST_EQUAL(s->field(2)->type()->Equals(QPXFeatureSchema::modificationsType()), true)
  // charge: int16, not null
  TEST_EQUAL(s->field(3)->name(), "charge")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(3)->nullable(), false)
  // posterior_error_probability: float64, nullable
  TEST_EQUAL(s->field(4)->name(), "posterior_error_probability")
  TEST_EQUAL(s->field(4)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(4)->nullable(), true)
  // is_decoy: bool, not null
  TEST_EQUAL(s->field(5)->name(), "is_decoy")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(5)->nullable(), false)
  // calculated_mz: float32, not null
  TEST_EQUAL(s->field(6)->name(), "calculated_mz")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(6)->nullable(), false)
  // observed_mz: float32, not null
  TEST_EQUAL(s->field(7)->name(), "observed_mz")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(7)->nullable(), false)
  // mass_error_ppm: float32, nullable
  TEST_EQUAL(s->field(8)->name(), "mass_error_ppm")
  TEST_EQUAL(s->field(8)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(8)->nullable(), true)
  // additional_scores: list<struct>, nullable
  TEST_EQUAL(s->field(9)->name(), "additional_scores")
  TEST_EQUAL(s->field(9)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(9)->nullable(), true)
  TEST_EQUAL(s->field(9)->type()->Equals(QPXFeatureSchema::additionalScoresType()), true)
  // predicted_rt: float32, nullable
  TEST_EQUAL(s->field(10)->name(), "predicted_rt")
  TEST_EQUAL(s->field(10)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(10)->nullable(), true)
  // run_file_name: utf8, not null
  TEST_EQUAL(s->field(11)->name(), "run_file_name")
  TEST_EQUAL(s->field(11)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(11)->nullable(), false)
  // cv_params: list<struct>, nullable
  TEST_EQUAL(s->field(12)->name(), "cv_params")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(12)->nullable(), true)
  TEST_EQUAL(s->field(12)->type()->Equals(QPXFeatureSchema::cvParamsType()), true)
  // scan: list<int32>, not null
  TEST_EQUAL(s->field(13)->name(), "scan")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(13)->nullable(), false)
  TEST_EQUAL(std::static_pointer_cast<arrow::ListType>(s->field(13)->type())->value_type()->id(), arrow::Type::INT32)
  // rt: float32, nullable
  TEST_EQUAL(s->field(14)->name(), "rt")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(14)->nullable(), true)
  // ion_mobility: float32, nullable
  TEST_EQUAL(s->field(15)->name(), "ion_mobility")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(15)->nullable(), true)
  // missed_cleavages: int16, nullable
  TEST_EQUAL(s->field(16)->name(), "missed_cleavages")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::INT16)
  TEST_EQUAL(s->field(16)->nullable(), true)
  // intensities: list<struct>, nullable
  TEST_EQUAL(s->field(17)->name(), "intensities")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(17)->nullable(), true)
  TEST_EQUAL(s->field(17)->type()->Equals(QPXFeatureSchema::intensitiesType()), true)
  // additional_intensities: list<struct>, nullable
  TEST_EQUAL(s->field(18)->name(), "additional_intensities")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(18)->nullable(), true)
  TEST_EQUAL(s->field(18)->type()->Equals(QPXFeatureSchema::additionalIntensitiesType()), true)
  // pg_accessions: list<struct>, nullable
  TEST_EQUAL(s->field(19)->name(), "pg_accessions")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(19)->nullable(), true)
  TEST_EQUAL(s->field(19)->type()->Equals(QPXFeatureSchema::pgAccessionsType()), true)
  // anchor_protein: utf8, not null
  TEST_EQUAL(s->field(20)->name(), "anchor_protein")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(20)->nullable(), false)
  // unique: bool, nullable
  TEST_EQUAL(s->field(21)->name(), "unique")
  TEST_EQUAL(s->field(21)->type()->id(), arrow::Type::BOOL)
  TEST_EQUAL(s->field(21)->nullable(), true)
  // pg_global_qvalue: float64, nullable
  TEST_EQUAL(s->field(22)->name(), "pg_global_qvalue")
  TEST_EQUAL(s->field(22)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(22)->nullable(), true)
  // pg_positions: list<struct>, nullable
  TEST_EQUAL(s->field(23)->name(), "pg_positions")
  TEST_EQUAL(s->field(23)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(23)->nullable(), true)
  TEST_EQUAL(s->field(23)->type()->Equals(QPXFeatureSchema::pgPositionsType()), true)
  // ion_mobility_start: float32, nullable
  TEST_EQUAL(s->field(24)->name(), "ion_mobility_start")
  TEST_EQUAL(s->field(24)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(24)->nullable(), true)
  // ion_mobility_stop: float32, nullable
  TEST_EQUAL(s->field(25)->name(), "ion_mobility_stop")
  TEST_EQUAL(s->field(25)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(25)->nullable(), true)
  // gg_accessions: list<utf8>, nullable
  TEST_EQUAL(s->field(26)->name(), "gg_accessions")
  TEST_EQUAL(s->field(26)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(26)->nullable(), true)
  // gg_names: list<utf8>, nullable
  TEST_EQUAL(s->field(27)->name(), "gg_names")
  TEST_EQUAL(s->field(27)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(27)->nullable(), true)
  // id_run_file_name: utf8, nullable
  TEST_EQUAL(s->field(28)->name(), "id_run_file_name")
  TEST_EQUAL(s->field(28)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(28)->nullable(), true)
  // rt_start: float32, nullable
  TEST_EQUAL(s->field(29)->name(), "rt_start")
  TEST_EQUAL(s->field(29)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(29)->nullable(), true)
  // rt_stop: float32, nullable
  TEST_EQUAL(s->field(30)->name(), "rt_stop")
  TEST_EQUAL(s->field(30)->type()->id(), arrow::Type::FLOAT)
  TEST_EQUAL(s->field(30)->nullable(), true)
}
END_SECTION

START_SECTION(QPXFeatureSchema validation with table)
{
  auto s = QPXFeatureSchema::schema();
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

// ========== ChromatogramSemiWideSchema ==========

START_SECTION(ChromatogramSemiWideSchema::schema() returns non-null with 6 fields)
{
  auto s = ChromatogramSemiWideSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 6)
}
END_SECTION

START_SECTION(ChromatogramSemiWideSchema column name constants)
{
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::CHROMATOGRAM_INDEX, "chromatogram_index")
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::NATIVE_ID, "native_id")
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::RT, "rt")
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::INTENSITY, "intensity")
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::PRECURSOR_MZ, "precursor_mz")
  TEST_STRING_EQUAL(ChromatogramSemiWideSchema::PRODUCT_MZ, "product_mz")
}
END_SECTION

START_SECTION(ChromatogramSemiWideSchema field types and nullability)
{
  auto s = ChromatogramSemiWideSchema::schema();
  // chromatogram_index: uint32, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "chromatogram_index")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::UINT32)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // native_id: utf8, nullable (default)
  TEST_EQUAL(s->field(1)->name(), "native_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // rt: list<float64>, nullable (default)
  TEST_EQUAL(s->field(2)->name(), "rt")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::LIST)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // intensity: list<float32>, nullable (default)
  TEST_EQUAL(s->field(3)->name(), "intensity")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::LIST)
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

START_SECTION(ChromatogramSemiWideSchema differs from ChromatogramSchema in rt/intensity types)
{
  auto long_s = ChromatogramSchema::schema();
  auto semi_s = ChromatogramSemiWideSchema::schema();
  // Long format: rt is scalar float64
  TEST_EQUAL(long_s->GetFieldByName("rt")->type()->id(), arrow::Type::DOUBLE)
  // Semi-wide format: rt is list<float64>
  TEST_EQUAL(semi_s->GetFieldByName("rt")->type()->id(), arrow::Type::LIST)
  // Long format: intensity is scalar float32
  TEST_EQUAL(long_s->GetFieldByName("intensity")->type()->id(), arrow::Type::FLOAT)
  // Semi-wide format: intensity is list<float32>
  TEST_EQUAL(semi_s->GetFieldByName("intensity")->type()->id(), arrow::Type::LIST)
}
END_SECTION

START_SECTION(ChromatogramSemiWideSchema subset validation with table)
{
  auto s = ChromatogramSemiWideSchema::schema();
  // Build a table with only a subset of fields
  auto subset_schema = arrow::schema({
    arrow::field("chromatogram_index", arrow::uint32()),
    arrow::field("rt", arrow::list(arrow::float64())),
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

START_SECTION(ChromatogramSemiWideSchema validation with full table)
{
  auto s = ChromatogramSemiWideSchema::schema();
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

// ========== OSWFeatureSchema ==========

START_SECTION(OSWFeatureSchema::schema() returns non-null with 65 fields)
{
  auto s = OSWFeatureSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 65)
}
END_SECTION

START_SECTION(OSWFeatureSchema column name constants)
{
  TEST_STRING_EQUAL(OSWFeatureSchema::FEATURE_ID, "feature_id")
  TEST_STRING_EQUAL(OSWFeatureSchema::RUN_ID, "run_id")
  TEST_STRING_EQUAL(OSWFeatureSchema::PRECURSOR_ID, "precursor_id")
  TEST_STRING_EQUAL(OSWFeatureSchema::EXP_RT, "exp_rt")
  TEST_STRING_EQUAL(OSWFeatureSchema::EXP_IM, "exp_im")
  TEST_STRING_EQUAL(OSWFeatureSchema::NORM_RT, "norm_rt")
  TEST_STRING_EQUAL(OSWFeatureSchema::DELTA_RT, "delta_rt")
  TEST_STRING_EQUAL(OSWFeatureSchema::LEFT_WIDTH, "left_width")
  TEST_STRING_EQUAL(OSWFeatureSchema::RIGHT_WIDTH, "right_width")
  TEST_STRING_EQUAL(OSWFeatureSchema::EXP_IM_LEFTWIDTH, "exp_im_leftwidth")
  TEST_STRING_EQUAL(OSWFeatureSchema::EXP_IM_RIGHTWIDTH, "exp_im_rightwidth")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS1_AREA_INTENSITY, "ms1_area_intensity")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS1_APEX_INTENSITY, "ms1_apex_intensity")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS1_EXP_IM, "ms1_exp_im")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS1_DELTA_IM, "ms1_delta_im")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_MASSDEV_SCORE, "var_ms1_massdev_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_IM_MS1_DELTA_SCORE, "var_ms1_im_ms1_delta_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_MI_SCORE, "var_ms1_mi_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_MI_CONTRAST_SCORE, "var_ms1_mi_contrast_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_MI_COMBINED_SCORE, "var_ms1_mi_combined_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_ISOTOPE_CORRELATION_SCORE, "var_ms1_isotope_correlation_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_ISOTOPE_OVERLAP_SCORE, "var_ms1_isotope_overlap_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_COELUTION, "var_ms1_xcorr_coelution")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_COELUTION_CONTRAST, "var_ms1_xcorr_coelution_contrast")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_COELUTION_COMBINED, "var_ms1_xcorr_coelution_combined")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_SHAPE, "var_ms1_xcorr_shape")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_SHAPE_CONTRAST, "var_ms1_xcorr_shape_contrast")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS1_XCORR_SHAPE_COMBINED, "var_ms1_xcorr_shape_combined")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_AREA_INTENSITY, "ms2_area_intensity")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_TOTAL_AREA_INTENSITY, "ms2_total_area_intensity")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_APEX_INTENSITY, "ms2_apex_intensity")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_EXP_IM, "ms2_exp_im")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_EXP_IM_LEFTWIDTH, "ms2_exp_im_leftwidth")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_EXP_IM_RIGHTWIDTH, "ms2_exp_im_rightwidth")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_DELTA_IM, "ms2_delta_im")
  TEST_STRING_EQUAL(OSWFeatureSchema::MS2_TOTAL_MI, "ms2_total_mi")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_BSERIES_SCORE, "var_ms2_bseries_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_DOTPROD_SCORE, "var_ms2_dotprod_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_INTENSITY_SCORE, "var_ms2_intensity_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_ISOTOPE_CORRELATION_SCORE, "var_ms2_isotope_correlation_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_ISOTOPE_OVERLAP_SCORE, "var_ms2_isotope_overlap_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_CORR, "var_ms2_library_corr")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_DOTPROD, "var_ms2_library_dotprod")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_MANHATTAN, "var_ms2_library_manhattan")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_RMSD, "var_ms2_library_rmsd")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_ROOTMEANSQUARE, "var_ms2_library_rootmeansquare")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LIBRARY_SANGLE, "var_ms2_library_sangle")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_LOG_SN_SCORE, "var_ms2_log_sn_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MANHATTAN_SCORE, "var_ms2_manhattan_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MASSDEV_SCORE, "var_ms2_massdev_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MASSDEV_SCORE_WEIGHTED, "var_ms2_massdev_score_weighted")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MI_SCORE, "var_ms2_mi_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MI_WEIGHTED_SCORE, "var_ms2_mi_weighted_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_MI_RATIO_SCORE, "var_ms2_mi_ratio_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_NORM_RT_SCORE, "var_ms2_norm_rt_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_XCORR_COELUTION, "var_ms2_xcorr_coelution")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_XCORR_COELUTION_WEIGHTED, "var_ms2_xcorr_coelution_weighted")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_XCORR_SHAPE, "var_ms2_xcorr_shape")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_XCORR_SHAPE_WEIGHTED, "var_ms2_xcorr_shape_weighted")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_YSERIES_SCORE, "var_ms2_yseries_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_ELUTION_MODEL_FIT_SCORE, "var_ms2_elution_model_fit_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_IM_XCORR_SHAPE, "var_ms2_im_xcorr_shape")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_IM_XCORR_COELUTION, "var_ms2_im_xcorr_coelution")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_IM_DELTA_SCORE, "var_ms2_im_delta_score")
  TEST_STRING_EQUAL(OSWFeatureSchema::VAR_MS2_IM_LOG_INTENSITY, "var_ms2_im_log_intensity")
}
END_SECTION

START_SECTION(OSWFeatureSchema field types - representative sample)
{
  auto s = OSWFeatureSchema::schema();
  // feature_id: int64, nullable (default)
  TEST_EQUAL(s->field(0)->name(), "feature_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // run_id: int64
  TEST_EQUAL(s->field(1)->name(), "run_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(1)->nullable(), true)
  // precursor_id: int64
  TEST_EQUAL(s->field(2)->name(), "precursor_id")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(2)->nullable(), true)
  // exp_rt: float64
  TEST_EQUAL(s->field(3)->name(), "exp_rt")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  TEST_EQUAL(s->field(3)->nullable(), true)
  // ms2_area_intensity: float64 (field 28)
  TEST_EQUAL(s->field(28)->name(), "ms2_area_intensity")
  TEST_EQUAL(s->field(28)->type()->id(), arrow::Type::DOUBLE)
  // var_ms2_im_log_intensity: float64 (last field, index 64)
  TEST_EQUAL(s->field(64)->name(), "var_ms2_im_log_intensity")
  TEST_EQUAL(s->field(64)->type()->id(), arrow::Type::DOUBLE)
}
END_SECTION

START_SECTION(OSWFeatureSchema validation with table)
{
  auto s = OSWFeatureSchema::schema();
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

// ========== OSWFeatureTransitionSchema ==========

START_SECTION(OSWFeatureTransitionSchema::schema() returns non-null with 44 fields)
{
  auto s = OSWFeatureTransitionSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 44)
}
END_SECTION

START_SECTION(OSWFeatureTransitionSchema column name constants)
{
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::FEATURE_ID, "feature_id")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::RUN_ID, "run_id")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::TRANSITION_ID, "transition_id")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::AREA_INTENSITY, "area_intensity")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::TOTAL_AREA_INTENSITY, "total_area_intensity")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::APEX_INTENSITY, "apex_intensity")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::APEX_RT, "apex_rt")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::RT_FWHM, "rt_fwhm")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::MASSERROR_PPM, "masserror_ppm")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::TOTAL_MI, "total_mi")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_INTENSITY_SCORE, "var_intensity_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_INTENSITY_RATIO_SCORE, "var_intensity_ratio_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_LOG_INTENSITY, "var_log_intensity")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_XCORR_COELUTION, "var_xcorr_coelution")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_XCORR_SHAPE, "var_xcorr_shape")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_LOG_SN_SCORE, "var_log_sn_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_MASSDEV_SCORE, "var_massdev_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_MI_SCORE, "var_mi_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_MI_RATIO_SCORE, "var_mi_ratio_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_ISOTOPE_CORRELATION_SCORE, "var_isotope_correlation_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_ISOTOPE_OVERLAP_SCORE, "var_isotope_overlap_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::EXP_IM, "exp_im")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::EXP_IM_LEFTWIDTH, "exp_im_leftwidth")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::EXP_IM_RIGHTWIDTH, "exp_im_rightwidth")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::DELTA_IM, "delta_im")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_DELTA_SCORE, "var_im_delta_score")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_LOG_INTENSITY, "var_im_log_intensity")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_XCORR_COELUTION_CONTRAST, "var_im_xcorr_coelution_contrast")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_XCORR_SHAPE_CONTRAST, "var_im_xcorr_shape_contrast")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_XCORR_COELUTION_COMBINED, "var_im_xcorr_coelution_combined")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::VAR_IM_XCORR_SHAPE_COMBINED, "var_im_xcorr_shape_combined")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::START_POSITION_AT_5, "start_position_at_5")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::END_POSITION_AT_5, "end_position_at_5")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::START_POSITION_AT_10, "start_position_at_10")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::END_POSITION_AT_10, "end_position_at_10")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::START_POSITION_AT_50, "start_position_at_50")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::END_POSITION_AT_50, "end_position_at_50")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::TOTAL_WIDTH, "total_width")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::TAILING_FACTOR, "tailing_factor")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::ASYMMETRY_FACTOR, "asymmetry_factor")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::SLOPE_OF_BASELINE, "slope_of_baseline")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::BASELINE_DELTA_2_HEIGHT, "baseline_delta_2_height")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::POINTS_ACROSS_BASELINE, "points_across_baseline")
  TEST_STRING_EQUAL(OSWFeatureTransitionSchema::POINTS_ACROSS_HALF_HEIGHT, "points_across_half_height")
}
END_SECTION

START_SECTION(OSWFeatureTransitionSchema field types - representative sample)
{
  auto s = OSWFeatureTransitionSchema::schema();
  // feature_id: int64
  TEST_EQUAL(s->field(0)->name(), "feature_id")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // run_id: int64
  TEST_EQUAL(s->field(1)->name(), "run_id")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::INT64)
  // transition_id: int64
  TEST_EQUAL(s->field(2)->name(), "transition_id")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT64)
  // area_intensity: float64
  TEST_EQUAL(s->field(3)->name(), "area_intensity")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::DOUBLE)
  // points_across_half_height: float64 (last field, index 43)
  TEST_EQUAL(s->field(43)->name(), "points_across_half_height")
  TEST_EQUAL(s->field(43)->type()->id(), arrow::Type::DOUBLE)
}
END_SECTION

START_SECTION(OSWFeatureTransitionSchema validation with table)
{
  auto s = OSWFeatureTransitionSchema::schema();
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

// ========== XICSchema ==========

START_SECTION(XICSchema::schema() returns non-null with 18 fields)
{
  auto s = XICSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 18)
}
END_SECTION

START_SECTION(XICSchema column name constants)
{
  TEST_STRING_EQUAL(XICSchema::RUN_ID, "RUN_ID")
  TEST_STRING_EQUAL(XICSchema::SOURCE_FILE, "SOURCE_FILE")
  TEST_STRING_EQUAL(XICSchema::MS_LEVEL, "MS_LEVEL")
  TEST_STRING_EQUAL(XICSchema::PRECURSOR_ID, "PRECURSOR_ID")
  TEST_STRING_EQUAL(XICSchema::TRANSITION_ID, "TRANSITION_ID")
  TEST_STRING_EQUAL(XICSchema::MODIFIED_SEQUENCE, "MODIFIED_SEQUENCE")
  TEST_STRING_EQUAL(XICSchema::PRECURSOR_CHARGE, "PRECURSOR_CHARGE")
  TEST_STRING_EQUAL(XICSchema::PRODUCT_CHARGE, "PRODUCT_CHARGE")
  TEST_STRING_EQUAL(XICSchema::DETECTING_TRANSITION, "DETECTING_TRANSITION")
  TEST_STRING_EQUAL(XICSchema::PRECURSOR_DECOY, "PRECURSOR_DECOY")
  TEST_STRING_EQUAL(XICSchema::PRODUCT_DECOY, "PRODUCT_DECOY")
  TEST_STRING_EQUAL(XICSchema::TRANSITION_ORDINAL, "TRANSITION_ORDINAL")
  TEST_STRING_EQUAL(XICSchema::TRANSITION_TYPE, "TRANSITION_TYPE")
  TEST_STRING_EQUAL(XICSchema::ANNOTATION, "ANNOTATION")
  TEST_STRING_EQUAL(XICSchema::RT_DATA, "RT_DATA")
  TEST_STRING_EQUAL(XICSchema::INTENSITY_DATA, "INTENSITY_DATA")
  TEST_STRING_EQUAL(XICSchema::RT_COMPRESSION, "RT_COMPRESSION")
  TEST_STRING_EQUAL(XICSchema::INTENSITY_COMPRESSION, "INTENSITY_COMPRESSION")
}
END_SECTION

START_SECTION(XICSchema field types)
{
  auto s = XICSchema::schema();
  // RUN_ID: int64
  TEST_EQUAL(s->field(0)->name(), "RUN_ID")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // SOURCE_FILE: utf8
  TEST_EQUAL(s->field(1)->name(), "SOURCE_FILE")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  // MS_LEVEL: int64
  TEST_EQUAL(s->field(2)->name(), "MS_LEVEL")
  TEST_EQUAL(s->field(2)->type()->id(), arrow::Type::INT64)
  // MODIFIED_SEQUENCE: utf8
  TEST_EQUAL(s->field(5)->name(), "MODIFIED_SEQUENCE")
  TEST_EQUAL(s->field(5)->type()->id(), arrow::Type::STRING)
  // TRANSITION_TYPE: utf8
  TEST_EQUAL(s->field(12)->name(), "TRANSITION_TYPE")
  TEST_EQUAL(s->field(12)->type()->id(), arrow::Type::STRING)
  // ANNOTATION: utf8
  TEST_EQUAL(s->field(13)->name(), "ANNOTATION")
  TEST_EQUAL(s->field(13)->type()->id(), arrow::Type::STRING)
  // RT_DATA: binary
  TEST_EQUAL(s->field(14)->name(), "RT_DATA")
  TEST_EQUAL(s->field(14)->type()->id(), arrow::Type::BINARY)
  // INTENSITY_DATA: binary
  TEST_EQUAL(s->field(15)->name(), "INTENSITY_DATA")
  TEST_EQUAL(s->field(15)->type()->id(), arrow::Type::BINARY)
  // RT_COMPRESSION: int64
  TEST_EQUAL(s->field(16)->name(), "RT_COMPRESSION")
  TEST_EQUAL(s->field(16)->type()->id(), arrow::Type::INT64)
  // INTENSITY_COMPRESSION: int64
  TEST_EQUAL(s->field(17)->name(), "INTENSITY_COMPRESSION")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::INT64)
}
END_SECTION

START_SECTION(XICSchema validation with table)
{
  auto s = XICSchema::schema();
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

// ========== XIMSchema ==========

START_SECTION(XIMSchema::schema() returns non-null with 21 fields)
{
  auto s = XIMSchema::schema();
  TEST_NOT_EQUAL(s, nullptr)
  TEST_EQUAL(s->num_fields(), 21)
}
END_SECTION

START_SECTION(XIMSchema column name constants)
{
  TEST_STRING_EQUAL(XIMSchema::RUN_ID, "RUN_ID")
  TEST_STRING_EQUAL(XIMSchema::SOURCE_FILE, "SOURCE_FILE")
  TEST_STRING_EQUAL(XIMSchema::MS_LEVEL, "MS_LEVEL")
  TEST_STRING_EQUAL(XIMSchema::MOBILOGRAM_TYPE, "MOBILOGRAM_TYPE")
  TEST_STRING_EQUAL(XIMSchema::PRECURSOR_ID, "PRECURSOR_ID")
  TEST_STRING_EQUAL(XIMSchema::TRANSITION_ID, "TRANSITION_ID")
  TEST_STRING_EQUAL(XIMSchema::FEATURE_ID, "FEATURE_ID")
  TEST_STRING_EQUAL(XIMSchema::FEATURE_RT, "FEATURE_RT")
  TEST_STRING_EQUAL(XIMSchema::MODIFIED_SEQUENCE, "MODIFIED_SEQUENCE")
  TEST_STRING_EQUAL(XIMSchema::PRECURSOR_CHARGE, "PRECURSOR_CHARGE")
  TEST_STRING_EQUAL(XIMSchema::PRODUCT_CHARGE, "PRODUCT_CHARGE")
  TEST_STRING_EQUAL(XIMSchema::DETECTING_TRANSITION, "DETECTING_TRANSITION")
  TEST_STRING_EQUAL(XIMSchema::PRECURSOR_DECOY, "PRECURSOR_DECOY")
  TEST_STRING_EQUAL(XIMSchema::PRODUCT_DECOY, "PRODUCT_DECOY")
  TEST_STRING_EQUAL(XIMSchema::TRANSITION_ORDINAL, "TRANSITION_ORDINAL")
  TEST_STRING_EQUAL(XIMSchema::TRANSITION_TYPE, "TRANSITION_TYPE")
  TEST_STRING_EQUAL(XIMSchema::ANNOTATION, "ANNOTATION")
  TEST_STRING_EQUAL(XIMSchema::MOBILITY_DATA, "MOBILITY_DATA")
  TEST_STRING_EQUAL(XIMSchema::INTENSITY_DATA, "INTENSITY_DATA")
  TEST_STRING_EQUAL(XIMSchema::MOBILITY_COMPRESSION, "MOBILITY_COMPRESSION")
  TEST_STRING_EQUAL(XIMSchema::INTENSITY_COMPRESSION, "INTENSITY_COMPRESSION")
}
END_SECTION

START_SECTION(XIMSchema field types)
{
  auto s = XIMSchema::schema();
  // RUN_ID: int64
  TEST_EQUAL(s->field(0)->name(), "RUN_ID")
  TEST_EQUAL(s->field(0)->type()->id(), arrow::Type::INT64)
  TEST_EQUAL(s->field(0)->nullable(), true)
  // SOURCE_FILE: utf8
  TEST_EQUAL(s->field(1)->name(), "SOURCE_FILE")
  TEST_EQUAL(s->field(1)->type()->id(), arrow::Type::STRING)
  // MOBILOGRAM_TYPE: utf8
  TEST_EQUAL(s->field(3)->name(), "MOBILOGRAM_TYPE")
  TEST_EQUAL(s->field(3)->type()->id(), arrow::Type::STRING)
  // FEATURE_ID: int64
  TEST_EQUAL(s->field(6)->name(), "FEATURE_ID")
  TEST_EQUAL(s->field(6)->type()->id(), arrow::Type::INT64)
  // FEATURE_RT: float64
  TEST_EQUAL(s->field(7)->name(), "FEATURE_RT")
  TEST_EQUAL(s->field(7)->type()->id(), arrow::Type::DOUBLE)
  // MOBILITY_DATA: binary
  TEST_EQUAL(s->field(17)->name(), "MOBILITY_DATA")
  TEST_EQUAL(s->field(17)->type()->id(), arrow::Type::BINARY)
  // INTENSITY_DATA: binary
  TEST_EQUAL(s->field(18)->name(), "INTENSITY_DATA")
  TEST_EQUAL(s->field(18)->type()->id(), arrow::Type::BINARY)
  // MOBILITY_COMPRESSION: int64
  TEST_EQUAL(s->field(19)->name(), "MOBILITY_COMPRESSION")
  TEST_EQUAL(s->field(19)->type()->id(), arrow::Type::INT64)
  // INTENSITY_COMPRESSION: int64
  TEST_EQUAL(s->field(20)->name(), "INTENSITY_COMPRESSION")
  TEST_EQUAL(s->field(20)->type()->id(), arrow::Type::INT64)
}
END_SECTION

START_SECTION(XIMSchema validation with table)
{
  auto s = XIMSchema::schema();
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

// ========== Timestamp unit compatibility ==========

START_SECTION(validate - Strict mode - timestamp unit compatibility)
{
  // Expected schema uses timestamp[s]
  auto expected = arrow::schema({
    arrow::field("ts", arrow::timestamp(arrow::TimeUnit::SECOND), true),
  });
  // Actual table uses timestamp[ms] (e.g., Parquet round-trip can change precision)
  auto actual_schema = arrow::schema({
    arrow::field("ts", arrow::timestamp(arrow::TimeUnit::MILLI), true),
  });
  auto arr = makeTimestampArray(arrow::timestamp(arrow::TimeUnit::MILLI), {1000});
  auto table = arrow::Table::Make(actual_schema, {arr});

  // Strict mode: different timestamp units should still pass (areTypesCompatible)
  auto result = ArrowSchemaValidation::validate(table, expected, ArrowSchemaValidation::Mode::Strict);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

START_SECTION(validate - Subset mode - timestamp unit compatibility)
{
  // Expected schema has multiple fields including a timestamp[s]
  auto expected = arrow::schema({
    arrow::field("id", arrow::utf8(), false),
    arrow::field("ts", arrow::timestamp(arrow::TimeUnit::SECOND), true),
  });
  // Actual table uses timestamp[us] for the timestamp field
  auto actual_schema = arrow::schema({
    arrow::field("ts", arrow::timestamp(arrow::TimeUnit::MICRO), true),
  });
  auto arr = makeTimestampArray(arrow::timestamp(arrow::TimeUnit::MICRO), {1000000});
  auto table = arrow::Table::Make(actual_schema, {arr});

  // Subset mode: timestamp with different unit should pass
  auto result = ArrowSchemaValidation::validate(table, expected, ArrowSchemaValidation::Mode::Subset);
  TEST_EQUAL(result.valid, true)
  TEST_EQUAL(result.errors.size(), 0)
}
END_SECTION

END_TEST
