// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzPeakFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzPeakFile* ptr = nullptr;
MzPeakFile* nullPointer = nullptr;
START_SECTION((MzPeakFile()))
ptr = new MzPeakFile;
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzPeakFile()))
delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, MapType& map))
{
  // ----------------------------------------------------------------------
  // Missing file must raise FileNotFound.
  // ----------------------------------------------------------------------
  {
    MSExperiment exp;
    TEST_EXCEPTION(Exception::FileNotFound, MzPeakFile().load("does_not_exist.mzpeak", exp))
  }

  // ----------------------------------------------------------------------
  // Load the bundled point-layout fixture (14 profile MS1 + 34 centroid MS2).
  // ----------------------------------------------------------------------
  MSExperiment exp;
  MzPeakFile().load(OPENMS_GET_TEST_DATA_PATH("small.mzpeak"), exp);

  // 48 spectra total.
  TEST_EQUAL(exp.size(), 48)

  // ms_level split: 14 MS1 + 34 MS2.
  Size n_ms1 = 0;
  Size n_ms2 = 0;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getMSLevel() == 1) ++n_ms1;
    else if (exp[i].getMSLevel() == 2)
      ++n_ms2;
  }
  TEST_EQUAL(n_ms1, 14)
  TEST_EQUAL(n_ms2, 34)

  // updateRanges() must have populated sane, non-default ranges.
  TEST_EQUAL(exp.getMinMZ() < exp.getMaxMZ(), true)
  TEST_EQUAL(exp.getMaxMZ() > 1000.0, true)

  // ----------------------------------------------------------------------
  // Profile spectrum (mzpeak index 0): native id "index=0", MS1, RT 0.004935 s.
  // After null-mz reconstruction it has 13589 points, strictly ascending m/z,
  // NO interior zero m/z; first/last m/z ~ 202.6066 / 1999.8404.
  // ----------------------------------------------------------------------
  const MSSpectrum* prof = nullptr;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getNativeID() == "index=0")
    {
      prof = &exp[i];
      break;
    }
  }
  TEST_NOT_EQUAL(prof, nullptr)
  if (prof)
  {
    TEST_EQUAL(prof->getMSLevel(), 1)
    TEST_REAL_SIMILAR(prof->getRT(), 0.004935)
    TEST_EQUAL(prof->size(), 13589)
    TEST_EQUAL(prof->empty(), false)

    // Strictly ascending m/z and no interior zero m/z (reconstruction worked).
    bool ascending = true;
    bool interior_zero = false;
    for (Size k = 1; k < prof->size(); ++k)
    {
      if ((*prof)[k].getMZ() <= (*prof)[k - 1].getMZ()) ascending = false;
      if ((*prof)[k].getMZ() == 0.0) interior_zero = true;
    }
    if ((*prof)[0].getMZ() == 0.0) interior_zero = true;
    TEST_EQUAL(ascending, true)
    TEST_EQUAL(interior_zero, false)

    TEST_REAL_SIMILAR((*prof)[0].getMZ(), 202.60657495520474)
    TEST_REAL_SIMILAR((*prof)[prof->size() - 1].getMZ(), 1999.8404377599534)

    // Checked intensity: point 1 has intensity 1938.1174 (point 0 intensity 0).
    TOLERANCE_ABSOLUTE(1e-2)
    TEST_REAL_SIMILAR((*prof)[1].getIntensity(), 1938.117431640625)
    TOLERANCE_ABSOLUTE(1e-4)
  }

  // ----------------------------------------------------------------------
  // Centroid spectrum (mzpeak index 2): native id "index=2", MS2, RT 0.011218 s.
  // 485 peaks; first peak (231.38884, 26.5451), last (1560.71985, 22.9731).
  // ----------------------------------------------------------------------
  const MSSpectrum* cent = nullptr;
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp[i].getNativeID() == "index=2")
    {
      cent = &exp[i];
      break;
    }
  }
  TEST_NOT_EQUAL(cent, nullptr)
  if (cent)
  {
    TEST_EQUAL(cent->getMSLevel(), 2)
    TEST_REAL_SIMILAR(cent->getRT(), 0.011218333333)
    TEST_EQUAL(cent->size(), 485)
    TEST_EQUAL(cent->empty(), false)

    // Peaks are sorted by m/z.
    bool sorted = true;
    for (Size k = 1; k < cent->size(); ++k)
    {
      if ((*cent)[k].getMZ() < (*cent)[k - 1].getMZ()) sorted = false;
    }
    TEST_EQUAL(sorted, true)

    TEST_REAL_SIMILAR((*cent)[0].getMZ(), 231.3888397216797)
    TEST_REAL_SIMILAR((*cent)[cent->size() - 1].getMZ(), 1560.7198486328125)

    TOLERANCE_ABSOLUTE(1e-3)
    TEST_REAL_SIMILAR((*cent)[0].getIntensity(), 26.54511260986328)
    TEST_REAL_SIMILAR((*cent)[cent->size() - 1].getIntensity(), 22.973094940185547)
    TOLERANCE_ABSOLUTE(1e-4)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
