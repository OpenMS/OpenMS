// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/SpectrumCount.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <map>

///////////////////////////

START_TEST(SpectrumCount, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumCount* ptr = nullptr;
SpectrumCount* null_ptr = nullptr;

START_SECTION((SpectrumCount()))
{
  ptr = new SpectrumCount();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~SpectrumCount()))
{
  delete ptr;
}
END_SECTION

START_SECTION((std::map<Size, UInt> compute(const MSExperiment& exp)))
{
  SpectrumCount sc;

  // 3 MS1, 5 MS2, 1 MS3 spectra
  MSExperiment exp;
  for (int i = 0; i < 3; ++i) { MSSpectrum s; s.setMSLevel(1); exp.addSpectrum(s); }
  for (int i = 0; i < 5; ++i) { MSSpectrum s; s.setMSLevel(2); exp.addSpectrum(s); }
  { MSSpectrum s; s.setMSLevel(3); exp.addSpectrum(s); }

  std::map<Size, UInt> counts = sc.compute(exp);
  TEST_EQUAL(counts.size(), 3)
  TEST_EQUAL(counts[1], 3)
  TEST_EQUAL(counts[2], 5)
  TEST_EQUAL(counts[3], 1)

  // empty experiment -> empty map
  MSExperiment empty;
  TEST_EQUAL(sc.compute(empty).empty(), true)
}
END_SECTION

START_SECTION((const std::string& getName() const override))
{
  SpectrumCount sc;
  TEST_STRING_EQUAL(sc.getName(), "SpectrumCount")
}
END_SECTION

START_SECTION((QCBase::Status requirements() const override))
{
  SpectrumCount sc;
  TEST_EQUAL(sc.requirements() == QCBase::Status(QCBase::Requires::RAWMZML), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
