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

END_TEST

#else // WITH_PARQUET

START_TEST(ArrowSchemaRegistry, "$Id$")
END_TEST

#endif
