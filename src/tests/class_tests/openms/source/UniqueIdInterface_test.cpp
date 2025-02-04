// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(UniqueIdInterface, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

UniqueIdInterface* ptr = nullptr;
UniqueIdInterface* nullPointer = nullptr;
START_SECTION(UniqueIdInterface())
{
	ptr = new UniqueIdInterface();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~UniqueIdInterface())
{
	delete ptr;
}
END_SECTION

START_SECTION((UniqueIdInterface(const UniqueIdInterface &rhs)))
{
  UniqueIdInterface uii1;
  UniqueIdInterface uii2(uii1);
  // to be continued further below
  NOT_TESTABLE
}
END_SECTION

START_SECTION((UniqueIdInterface(UniqueIdInterface &&rhs)))
{
  // Ensure that UniqueIdInterface has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(UniqueIdInterface(std::declval<UniqueIdInterface&&>())), true)
}
END_SECTION

START_SECTION((UniqueIdInterface& operator=(UniqueIdInterface const &rhs)))
{
  UniqueIdInterface uii1;
  UniqueIdInterface uii2;
  uii2 = uii1;
  // to be continued further below
  NOT_TESTABLE
}
END_SECTION

START_SECTION([EXTRA](~UniqueIdInterface()))
{
  {
    UniqueIdInterface uii1;
    UniqueIdInterface uii2;
    uii2.setUniqueId(17);
    // destructor called when the scope is left
  }
  // to be continued further below ;-)
  NOT_TESTABLE
}
END_SECTION

START_SECTION((bool operator==(UniqueIdInterface const &rhs) const ))
{
  UniqueIdInterface uii1;
  UniqueIdInterface uii2;
  TEST_EQUAL(uii1 == uii2,true);
  UniqueIdInterface uii3(uii1);
  TEST_EQUAL(uii1 == uii3,true);
  uii2.setUniqueId(17);
  TEST_EQUAL(uii1 == uii2,false);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1 == uii2,true);
}
END_SECTION

START_SECTION((UInt64 getUniqueId() const ))
{
  UniqueIdInterface uii1;
  TEST_EQUAL(uii1.getUniqueId(),0);
  UniqueIdInterface uii2;
  uii2 = uii1;
  TEST_EQUAL(uii2.getUniqueId(),0);
  UniqueIdInterface uii3(uii1);
  TEST_EQUAL(uii3.getUniqueId(),0);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.getUniqueId(),17);
  uii2 = uii1;
  TEST_EQUAL(uii2.getUniqueId(),17);
  UniqueIdInterface uii4(uii1);
  TEST_EQUAL(uii4.getUniqueId(),17);
}
END_SECTION

START_SECTION((Size clearUniqueId()))
{
  UniqueIdInterface uii1;
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.getUniqueId(),0)
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.getUniqueId(),17);
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.getUniqueId(),0)
}
END_SECTION

START_SECTION((void swap(UniqueIdInterface &from)))
{
	UniqueIdInterface u1;
	u1.setUniqueId(111);
	UniqueIdInterface u2;
	u2.setUniqueId(222);
	u1.swap(u2);
	TEST_EQUAL(u1.getUniqueId(),222);
	TEST_EQUAL(u2.getUniqueId(),111);
	std::swap(u1,u2);
	TEST_EQUAL(u1.getUniqueId(),111);
	TEST_EQUAL(u2.getUniqueId(),222);
}
END_SECTION


START_SECTION((static bool isValid(UInt64 unique_id)))
{
  TEST_EQUAL(UniqueIdInterface::isValid(UniqueIdInterface::INVALID),false);
  TEST_EQUAL(UniqueIdInterface::isValid(1234567890),true);
}
END_SECTION

START_SECTION((Size hasValidUniqueId() const ))
{
  UniqueIdInterface uii1;
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId(0);
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
}
END_SECTION

START_SECTION((Size hasInvalidUniqueId() const ))
{
  UniqueIdInterface uii1;
  TEST_EQUAL(uii1.hasInvalidUniqueId(),true);
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.hasInvalidUniqueId(),true);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.hasInvalidUniqueId(),false);
  uii1.setUniqueId(0);
  TEST_EQUAL(uii1.hasInvalidUniqueId(),true);
  uii1.setUniqueId(17);
  TEST_EQUAL(uii1.hasInvalidUniqueId(),false);
  uii1.clearUniqueId();
  TEST_EQUAL(uii1.hasInvalidUniqueId(),true);
}
END_SECTION

START_SECTION((Size setUniqueId()))
{
  UniqueIdInterface uii1;
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
  uii1.setUniqueId();
  TEST_EQUAL(uii1.hasValidUniqueId(),true);
}
END_SECTION

START_SECTION((Size ensureUniqueId()))
{
  UInt64 the_unique_id;
  UniqueIdInterface uii1;
  the_unique_id = uii1.getUniqueId();
  TEST_EQUAL(the_unique_id,0);
  uii1.ensureUniqueId();
  the_unique_id = uii1.getUniqueId();
  TEST_NOT_EQUAL(the_unique_id,0);
  uii1.ensureUniqueId();
  TEST_EQUAL(uii1.getUniqueId(),the_unique_id);
}
END_SECTION

START_SECTION((void setUniqueId(UInt64 rhs)))
{
  // that one was used throughout the other tests
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((void setUniqueId(const String &rhs)))
{
  UniqueIdInterface uii1;
  TEST_EQUAL(uii1.getUniqueId(),0);
  uii1.setUniqueId("");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId("_");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId("17");
  TEST_EQUAL(uii1.getUniqueId(),17);
  uii1.setUniqueId("18");
  TEST_EQUAL(uii1.getUniqueId(),18);
  uii1.setUniqueId("asdf_19");
  TEST_EQUAL(uii1.getUniqueId(),19);
  uii1.setUniqueId("_20");
  TEST_EQUAL(uii1.getUniqueId(),20);
  uii1.setUniqueId("_021");
  TEST_EQUAL(uii1.getUniqueId(),21);
  uii1.setUniqueId("asdf_19_22");
  TEST_EQUAL(uii1.getUniqueId(),22);
  uii1.setUniqueId("_20_23");
  TEST_EQUAL(uii1.getUniqueId(),23);
  uii1.setUniqueId("_021_024");
  TEST_EQUAL(uii1.getUniqueId(),24);
  uii1.setUniqueId("   _021_025     ");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId("bla");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123_ 456");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 _456");
  TEST_EQUAL(uii1.getUniqueId(),456);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 _ 456");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456_");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456_  ");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456 _ ");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456  _");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("123 456 ff");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("_021bla_");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("_021 bla    ");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("_021 bla");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("_021_bla");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
  uii1.setUniqueId(1000000);
  uii1.setUniqueId("_021_bla_");
  TEST_EQUAL(uii1.hasValidUniqueId(),false);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



