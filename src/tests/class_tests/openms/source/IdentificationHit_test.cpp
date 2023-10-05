// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/IdentificationHit.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IdentificationHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IdentificationHit* ptr = nullptr;
IdentificationHit* null_ptr = nullptr;
START_SECTION(IdentificationHit())
{
	ptr = new IdentificationHit();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IdentificationHit())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~IdentificationHit()))
{
  // TODO
}
END_SECTION

START_SECTION((IdentificationHit(const IdentificationHit &source)))
{
  // TODO
}
END_SECTION

START_SECTION((IdentificationHit& operator=(const IdentificationHit &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IdentificationHit &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IdentificationHit &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setId(const String &id)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getId() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCharge(Int charge)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCalculatedMassToCharge(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getCalculatedMassToCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setExperimentalMassToCharge(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getExperimentalMassToCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPassThreshold(bool pass)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getPassThreshold() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRank(Int rank)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getRank() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



