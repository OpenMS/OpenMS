// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFileOptions, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFileOptions* ptr = nullptr;
FeatureFileOptions* null_ptr = nullptr;
START_SECTION(FeatureFileOptions())
{
	ptr = new FeatureFileOptions();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFileOptions())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setIntensityRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasIntensityRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getIntensityRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setLoadConvexHull(bool convex)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getLoadConvexHull() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setLoadSubordinates(bool sub)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getLoadSubordinates() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMetadataOnly(bool only)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getMetadataOnly() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setSizeOnly(bool only)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getSizeOnly() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRTRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasRTRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getRTRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMZRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasMZRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getMZRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((~FeatureFileOptions()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



