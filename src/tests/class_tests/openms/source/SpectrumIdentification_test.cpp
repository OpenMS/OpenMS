// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/SpectrumIdentification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumIdentification* ptr = nullptr;
SpectrumIdentification* null_ptr = nullptr;
START_SECTION(SpectrumIdentification())
{
	ptr = new SpectrumIdentification();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~SpectrumIdentification())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~SpectrumIdentification()))
{
  // TODO
}
END_SECTION

START_SECTION((SpectrumIdentification(const SpectrumIdentification &source)))
{
  // TODO
}
END_SECTION

START_SECTION((SpectrumIdentification& operator=(const SpectrumIdentification &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const SpectrumIdentification &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const SpectrumIdentification &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setHits(const std::vector< IdentificationHit > &hits)))
{
  // TODO
}
END_SECTION

START_SECTION((void addHit(const IdentificationHit &hit)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<IdentificationHit>& getHits() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



