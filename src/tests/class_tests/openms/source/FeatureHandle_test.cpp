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
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap ContainerType;
typedef ContainerType::value_type ElementType;
typedef Feature::PositionType PositionType;

START_TEST(FeatureHandle, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureHandle* ptr = nullptr;
FeatureHandle* nullPointer = nullptr;
START_SECTION((FeatureHandle()))
	ptr = new FeatureHandle();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureHandle()))
	delete ptr;
END_SECTION

START_SECTION((FeatureHandle& operator=(const FeatureHandle &rhs)))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  FeatureHandle it_copy;
  it_copy = it;

  TEST_EQUAL(it.getUniqueId() == it_copy.getUniqueId(), true)
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)
END_SECTION

START_SECTION((FeatureHandle(const FeatureHandle &rhs)))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  FeatureHandle it_copy(it);

  TEST_EQUAL(it.getUniqueId() == it_copy.getUniqueId(), true)
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)
END_SECTION

START_SECTION((void setCharge(ChargeType charge)))
{
  FeatureHandle fh;
  fh.setCharge(-17);
  TEST_EQUAL(fh.getCharge(),-17);
  fh.setCharge(-1717);
  TEST_EQUAL(fh.getCharge(),-1717);
}
END_SECTION

START_SECTION((ChargeType getCharge() const))
{
  NOT_TESTABLE; // see setCharge()
}
END_SECTION


START_SECTION((void setWidth(WidthType width)))
{
    FeatureHandle fh_tmp;
    fh_tmp.setWidth(10.7);
    TEST_REAL_SIMILAR(fh_tmp.getWidth(), 10.7);
    fh_tmp.setWidth(-8.9);
    TEST_REAL_SIMILAR(fh_tmp.getWidth(), -8.9);
}
END_SECTION

START_SECTION((WidthType getWidth() const ))
{
    NOT_TESTABLE;
}
END_SECTION


START_SECTION((FeatureHandle(UInt64 map_index, const Peak2D &point, UInt64 element_index)))
  ElementType e;
  FeatureHandle it(1,e,2);

  TEST_EQUAL(it.getUniqueId() == 2, true)
  TEST_EQUAL(it.getMapIndex() == 1, true)
  TEST_EQUAL(it.getPosition() == e.getPosition(), true)
END_SECTION

START_SECTION((FeatureHandle(UInt64 map_index, const BaseFeature& feature)))

  Feature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setUniqueId(23);
  const Feature& f_cref = f;
  FeatureHandle fh(99,f_cref);

  TEST_EQUAL(fh.getMapIndex(),99);
  TEST_EQUAL(fh.getUniqueId(),23);
  TEST_EQUAL(fh.getRT(),44324.6);
  TEST_EQUAL(fh.getMZ(),867.4);
  TEST_EQUAL(fh.getCharge(),-17);

END_SECTION

START_SECTION((FeatureHandleMutable_ & asMutable() const))
  ConsensusFeature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setUniqueId(23);
 const ConsensusFeature& f_cref = f;
  FeatureHandle fh(99, f_cref);

  const FeatureHandle& fh_cref = fh;
  // fh_cref.setRT(-64544.3); // compile time error
  fh_cref.asMutable().setRT(-64544.3); // ok

  TEST_EQUAL(fh.getMapIndex(),99);
  TEST_EQUAL(fh.getUniqueId(),23);
  TEST_EQUAL(fh.getRT(),-64544.3);
  TEST_EQUAL(fh.getMZ(),867.4);
  TEST_EQUAL(fh.getCharge(),-17);

END_SECTION


START_SECTION((bool operator!=(const FeatureHandle &i) const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it1(1,e);
  FeatureHandle it2(2,e);

  TEST_FALSE(it1 == it2)
END_SECTION

START_SECTION((bool operator==(const FeatureHandle &i) const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it1(2,e);
  FeatureHandle it2(2,e);

  TEST_TRUE(it1 == it2)
END_SECTION

START_SECTION((UInt64 getMapIndex() const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)
END_SECTION

START_SECTION((void setMapIndex(UInt64 i)))
  FeatureHandle it;
  it.setMapIndex(2);
  it.setUniqueId(77);

  TEST_EQUAL(it.getMapIndex() == 2, true)
END_SECTION

START_SECTION(([FeatureHandle::IndexLess] bool operator()(FeatureHandle const &left, FeatureHandle const &right) const))
  FeatureHandle lhs, rhs;
  lhs.setMapIndex(2);
  lhs.setUniqueId(77);
  rhs.setMapIndex(4);
  lhs.setUniqueId(29);

  FeatureHandle::IndexLess il;

  TEST_EQUAL(il(lhs, rhs), 1);
  TEST_EQUAL(il(rhs, lhs), 0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



