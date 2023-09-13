// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Identification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Identification* ptr = nullptr;
Identification* null_ptr = nullptr;
START_SECTION(Identification())
{
	ptr = new Identification();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Identification())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~Identification()))
{
  // TODO
}
END_SECTION

START_SECTION((Identification(const Identification &source)))
{
  // TODO
}
END_SECTION

START_SECTION((Identification& operator=(const Identification &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const Identification &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const Identification &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCreationDate(const DateTime &date)))
{
  // TODO
}
END_SECTION

START_SECTION((const DateTime& getCreationDate() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setSpectrumIdentifications(const std::vector< SpectrumIdentification > &ids)))
{
  // TODO
}
END_SECTION

START_SECTION((void addSpectrumIdentification(const SpectrumIdentification &id)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<SpectrumIdentification>& getSpectrumIdentifications() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



