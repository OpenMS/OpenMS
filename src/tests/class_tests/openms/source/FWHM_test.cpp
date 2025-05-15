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
#include <OpenMS/QC/FWHM.h>

///////////////////////////

START_TEST(FWHM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

FWHM* ptr = nullptr;
FWHM* nullPointer = nullptr;
START_SECTION(MzCalibration())
ptr = new FWHM();
TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~FWHM())
delete ptr;
END_SECTION


START_SECTION(void compute(FeatureMap& features))
{
  Feature f;
  PeptideIdentification pi;
  pi.getHits().push_back(PeptideHit(1.0, 1, 3, AASequence::fromString("KKK")));
  f.getPeptideIdentifications().push_back(pi);
  f.setMetaValue("FWHM", 123.4);
  FeatureMap fm;
  fm.push_back(f);
  f.clearMetaInfo();
  f.setMetaValue("model_FWHM", 98.1);
  fm.push_back(f);
  FWHM fw;
  fw.compute(fm);
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getMetaValue("FWHM"), 123.4)
  TEST_EQUAL(fm[1].getPeptideIdentifications()[0].getMetaValue("FWHM"), 98.1)
}
END_SECTION

START_SECTION(QCBase::Status requirements() const override)
{
  FWHM fw;
  TEST_EQUAL(fw.requirements() == (QCBase::Status() | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(const String& getName() const)
{
  TEST_EQUAL(FWHM().getName(), "FWHM");
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
