// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/MzQuantMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MzQuantMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzQuantMLFile* ptr = nullptr;
MzQuantMLFile* null_ptr = nullptr;
START_SECTION(MzQuantMLFile())
{
	ptr = new MzQuantMLFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MzQuantMLFile())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~MzQuantMLFile()))
{
  // TODO
}
END_SECTION

START_SECTION((void load(const String &filename, MSQuantifications &msq)))
{
  // TODO
}
END_SECTION

START_SECTION((void store(const String &filename, const MSQuantifications &cmsq) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool isSemanticallyValid(const String &filename, StringList &errors, StringList &warnings)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



