// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/PeptideEvidence.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeptideEvidence, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideEvidence* ptr = nullptr;
PeptideEvidence* null_ptr = nullptr;
START_SECTION(PeptideEvidence())
{
	ptr = new PeptideEvidence();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeptideEvidence())
{
	delete ptr;
}
END_SECTION

START_SECTION((PeptideEvidence(const PeptideEvidence &source)))
{
  // TODO
}
END_SECTION

START_SECTION((~PeptideEvidence()))
{
  // TODO
}
END_SECTION

START_SECTION((PeptideEvidence& operator=(const PeptideEvidence &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const PeptideEvidence &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const PeptideEvidence &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getProteinAccession() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProteinAccession(const String &s)))
{
  // TODO
}
END_SECTION

START_SECTION((void setStart(const Int a)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getStart() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setEnd(const Int a)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getEnd() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setAABefore(const char acid)))
{
  // TODO
}
END_SECTION

START_SECTION((char getAABefore() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setAAAfter(const char acid)))
{
  // TODO
}
END_SECTION

START_SECTION((char getAAAfter() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



