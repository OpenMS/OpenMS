#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#ifdef WITH_PARQUET

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <arrow/api.h>

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

END_TEST

#else // WITH_PARQUET

START_TEST(ArrowSchemaRegistry, "$Id$")
END_TEST

#endif
