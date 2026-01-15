// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSpectralContrastAngle, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinnedSpectralContrastAngle* ptr = nullptr;
BinnedSpectralContrastAngle* nullPointer = nullptr;
START_SECTION(BinnedSpectralContrastAngle())
{
  ptr = new BinnedSpectralContrastAngle();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~BinnedSpectralContrastAngle())
{
  delete ptr;
}
END_SECTION

ptr = new BinnedSpectralContrastAngle();

START_SECTION((BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle &source)))
{
  BinnedSpectralContrastAngle copy(*ptr);
  TEST_EQUAL(copy.getName(), ptr->getName());
  TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((BinnedSpectralContrastAngle& operator=(const BinnedSpectralContrastAngle &source)))
{
  BinnedSpectralContrastAngle copy;
  copy = *ptr;
  TEST_EQUAL(copy.getName(), ptr->getName());
  TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
{
  PeakSpectrum s1, s2;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
  s2.pop_back();
  BinnedSpectrum bs1 (s1, 1.5, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
  BinnedSpectrum bs2 (s2, 1.5, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
  double score = (*ptr)(bs1, bs2);
  TEST_REAL_SIMILAR(score,0.999985)
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  BinnedSpectrum bs1 (s1, 1.5, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
  double score = (*ptr)(bs1);
  TEST_REAL_SIMILAR(score, 1);
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

