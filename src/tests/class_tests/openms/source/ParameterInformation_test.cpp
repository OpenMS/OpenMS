// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/ParameterInformation.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(ParameterInformation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ParameterInformation* ptr = nullptr;
ParameterInformation* null_ptr = nullptr;
START_SECTION(ParameterInformation())
{
	ptr = new ParameterInformation();
	TEST_NOT_EQUAL(ptr, null_ptr)

	TEST_EQUAL(ptr->name, "")
  TEST_EQUAL(ptr->type, ParameterInformation::NONE)
  TEST_EQUAL(ptr->default_value, ParamValue())
	TEST_EQUAL(ptr->description, "")
  TEST_EQUAL(ptr->argument, "")
  TEST_EQUAL(ptr->required, true)
  TEST_EQUAL(ptr->advanced, false)
  TEST_EQUAL(ptr->tags.size(), 0)
  TEST_EQUAL(ptr->valid_strings.size(), 0)
  TEST_EQUAL(ptr->min_int, -std::numeric_limits<Int>::max())
  TEST_EQUAL(ptr->max_int, std::numeric_limits<Int>::max())
  TEST_EQUAL(ptr->min_float, -std::numeric_limits<double>::max())
  TEST_EQUAL(ptr->max_float, std::numeric_limits<double>::max())
}
END_SECTION

START_SECTION(~ParameterInformation())
{
	delete ptr;
}
END_SECTION

START_SECTION((ParameterInformation(const String &n, ParameterTypes t, const String &arg, const DataValue &def, const String &desc, bool req, bool adv, const StringList &tag_values=StringList())))
{
  ParameterInformation pi1("pi1_name", ParameterInformation::STRING, "<STRING>", "def_value", "this is a description", false, true, ListUtils::create<String>("tag1,tag2"));

	TEST_EQUAL(pi1.name, "pi1_name")
  TEST_EQUAL(pi1.type, ParameterInformation::STRING)
  TEST_EQUAL(pi1.default_value, "def_value")
	TEST_EQUAL(pi1.description, "this is a description")
  TEST_EQUAL(pi1.argument, "<STRING>")
  TEST_EQUAL(pi1.required, false)
  TEST_EQUAL(pi1.advanced, true)
  TEST_EQUAL(pi1.tags.size(), 2)
  ABORT_IF(pi1.tags.size() != 2)
  TEST_EQUAL(pi1.tags[0], "tag1")
  TEST_EQUAL(pi1.tags[1], "tag2")

  TEST_EQUAL(pi1.valid_strings.size(), 0)
  TEST_EQUAL(pi1.min_int, -std::numeric_limits<Int>::max())
  TEST_EQUAL(pi1.max_int, std::numeric_limits<Int>::max())
  TEST_EQUAL(pi1.min_float, -std::numeric_limits<double>::max())
  TEST_EQUAL(pi1.max_float, std::numeric_limits<double>::max())
}
END_SECTION

START_SECTION((ParameterInformation& operator=(const ParameterInformation &rhs)))
{
 ParameterInformation pi1("pi1_name", ParameterInformation::STRING, "<STRING>", "def_value", "this is a description", false, true, ListUtils::create<String>("tag1,tag2"));

	TEST_EQUAL(pi1.name, "pi1_name")
  TEST_EQUAL(pi1.type, ParameterInformation::STRING)
  TEST_EQUAL(pi1.default_value, "def_value")
	TEST_EQUAL(pi1.description, "this is a description")
  TEST_EQUAL(pi1.argument, "<STRING>")
  TEST_EQUAL(pi1.required, false)
  TEST_EQUAL(pi1.advanced, true)
  TEST_EQUAL(pi1.tags.size(), 2)
  ABORT_IF(pi1.tags.size() != 2)
  TEST_EQUAL(pi1.tags[0], "tag1")
  TEST_EQUAL(pi1.tags[1], "tag2")

  TEST_EQUAL(pi1.valid_strings.size(), 0)
  TEST_EQUAL(pi1.min_int, -std::numeric_limits<Int>::max())
  TEST_EQUAL(pi1.max_int, std::numeric_limits<Int>::max())
  TEST_EQUAL(pi1.min_float, -std::numeric_limits<double>::max())
  TEST_EQUAL(pi1.max_float, std::numeric_limits<double>::max())

	ParameterInformation pi2;
	pi2 = pi1;

	TEST_EQUAL(pi2.name, "pi1_name")
  TEST_EQUAL(pi2.type, ParameterInformation::STRING)
  TEST_EQUAL(pi2.default_value, "def_value")
	TEST_EQUAL(pi2.description, "this is a description")
  TEST_EQUAL(pi2.argument, "<STRING>")
  TEST_EQUAL(pi2.required, false)
  TEST_EQUAL(pi2.advanced, true)
  TEST_EQUAL(pi2.tags.size(), 2)
  ABORT_IF(pi2.tags.size() != 2)
  TEST_EQUAL(pi2.tags[0], "tag1")
  TEST_EQUAL(pi2.tags[1], "tag2")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
