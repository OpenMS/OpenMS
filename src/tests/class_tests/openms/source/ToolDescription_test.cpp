// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ToolDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolDescription* ptr = 0;
ToolDescription* null_ptr = 0;
START_SECTION(ToolDescription())
{
	ptr = new ToolDescription();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ToolDescription())
{
	delete ptr;
}
END_SECTION

START_SECTION((ToolDescription(const String &p_name, const String &p_category, const StringList &p_types=StringList())))
{
  // TODO
}
END_SECTION

START_SECTION((void addExternalType(const String &type, const ToolExternalDetails &details)))
{
  // TODO
}
END_SECTION

START_SECTION((void append(const ToolDescription &other)))
{
  // TODO
}
END_SECTION

START_SECTION((ToolDescription& operator=(const ToolDescription &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



