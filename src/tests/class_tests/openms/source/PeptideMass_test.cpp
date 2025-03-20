// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/PeptideMass.h>

///////////////////////////

START_TEST(PeptideMass, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

PeptideMass* ptr = nullptr;
PeptideMass* nullPointer = nullptr;
START_SECTION(MzCalibration())
ptr = new PeptideMass();
TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~PeptideMass())
delete ptr;
END_SECTION


START_SECTION(void compute(FeatureMap& features))
{
  Feature f;
  PeptideIdentification pi;
  pi.getHits().push_back(PeptideHit(1.0, 1, 3, AASequence::fromString("KKK")));
  pi.setMZ(100.0);
  f.getPeptideIdentifications().push_back(pi);
  FeatureMap fm;
  fm.push_back(f);
  pi.setMZ(200.0);
  pi.getHits().back().setCharge(2);
  f.getPeptideIdentifications().back() = pi;
  fm.push_back(f);
  PeptideMass fw;
  fw.compute(fm);
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mass"), (100.0 - Constants::PROTON_MASS_U) * 3)
  TEST_EQUAL(fm[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mass"), (200.0 - Constants::PROTON_MASS_U) * 2)
}
END_SECTION

START_SECTION(QCBase::Status requirements() const override)
{
  PeptideMass fw;
  TEST_EQUAL(fw.requirements() == (QCBase::Status() | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(const String& getName() const)
{
  TEST_EQUAL(PeptideMass().getName(), "PeptideMass");
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
