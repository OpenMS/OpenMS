// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
///////////////////////////

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

START_TEST(GridFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

GridFeature* gf_ptr = nullptr;
GridFeature* gf_nullPointer = nullptr;

START_SECTION((GridFeature(const BaseFeature& feature, Size map_index, Size feature_index)))
{
  BaseFeature bf;
  gf_ptr = new GridFeature(bf, 0, 0);
  TEST_NOT_EQUAL(gf_ptr, gf_nullPointer);
}
END_SECTION

START_SECTION((~GridFeature()))
{
  delete gf_ptr;
}
END_SECTION

START_SECTION((const BaseFeature& getFeature() const))
{
  BaseFeature bf;
  bf.setRT(1.1);
  bf.setMZ(2.2);
  bf.setCharge(3);
  const BaseFeature bf_const(bf);
  GridFeature gf(bf_const, 0, 0);
  TEST_EQUAL(gf.getFeature() == bf_const, true);
}
END_SECTION

START_SECTION((Size getMapIndex() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 123, 0);
  TEST_EQUAL(gf.getMapIndex(), 123);
}
END_SECTION

START_SECTION((Size getFeatureIndex() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 123);
  TEST_EQUAL(gf.getFeatureIndex(), 123);
}
END_SECTION

START_SECTION((Int getID() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 123);
  TEST_EQUAL(gf.getID(), 123);
}
END_SECTION

START_SECTION((const std::set<AASequence>& getAnnotations() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 0);
  TEST_EQUAL(gf.getAnnotations().size(), 0);
  bf.getPeptideIdentifications().resize(2);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AAA"));
  bf.getPeptideIdentifications()[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("CCC"));
  bf.getPeptideIdentifications()[1].insertHit(hit);
  GridFeature gf2(bf, 0, 0);
  TEST_EQUAL(gf2.getAnnotations().size(), 2);
  TEST_EQUAL(*(gf2.getAnnotations().begin()), AASequence::fromString("AAA"));
  TEST_EQUAL(*(gf2.getAnnotations().rbegin()), AASequence::fromString("CCC"));
}
END_SECTION

START_SECTION((double getRT() const))
{
  BaseFeature bf;
  bf.setRT(4.56);
  GridFeature gf(bf, 0, 123);
  TEST_REAL_SIMILAR(gf.getRT(), 4.56);
}
END_SECTION

START_SECTION((double getMZ() const))
{
  BaseFeature bf;
  bf.setMZ(4.56);
  GridFeature gf(bf, 0, 123);
  TEST_REAL_SIMILAR(gf.getMZ(), 4.56);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
