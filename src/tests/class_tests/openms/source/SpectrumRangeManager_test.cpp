// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/SpectrumRangeManager.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <set>

///////////////////////////

START_TEST(SpectrumRangeManager, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumRangeManager* ptr = nullptr;
SpectrumRangeManager* null_ptr = nullptr;

START_SECTION((SpectrumRangeManager()))
{
  ptr = new SpectrumRangeManager();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getMSLevels().empty(), true)
}
END_SECTION

START_SECTION((~SpectrumRangeManager()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void extendRT(double rt, UInt ms_level = 0)))
{
  SpectrumRangeManager srm;
  // ms_level 0 -> global (base) range, does not register an MS level
  srm.extendRT(100.0);
  srm.extendRT(200.0);
  TEST_REAL_SIMILAR(srm.getMinRT(), 100.0)
  TEST_REAL_SIMILAR(srm.getMaxRT(), 200.0)
  TEST_EQUAL(srm.getMSLevels().empty(), true)

  // ms_level 2 -> level-specific range, does not touch the global range
  srm.extendRT(10.0, 2);
  srm.extendRT(20.0, 2);
  TEST_REAL_SIMILAR(srm.byMSLevel(2).getMinRT(), 10.0)
  TEST_REAL_SIMILAR(srm.byMSLevel(2).getMaxRT(), 20.0)
  TEST_REAL_SIMILAR(srm.getMinRT(), 100.0) // global unchanged
  TEST_REAL_SIMILAR(srm.getMaxRT(), 200.0)
}
END_SECTION

START_SECTION((void extendMZ(double mz, UInt ms_level = 0)))
{
  SpectrumRangeManager srm;
  srm.extendMZ(500.0);
  srm.extendMZ(600.0);
  TEST_REAL_SIMILAR(srm.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(srm.getMaxMZ(), 600.0)

  srm.extendMZ(300.0, 2);
  srm.extendMZ(400.0, 2);
  TEST_REAL_SIMILAR(srm.byMSLevel(2).getMinMZ(), 300.0)
  TEST_REAL_SIMILAR(srm.byMSLevel(2).getMaxMZ(), 400.0)
  TEST_REAL_SIMILAR(srm.getMinMZ(), 500.0) // global unchanged
}
END_SECTION

START_SECTION((std::set<UInt> getMSLevels() const))
{
  SpectrumRangeManager srm;
  TEST_EQUAL(srm.getMSLevels().empty(), true)

  srm.extendRT(1.0, 2);
  srm.extendRT(1.0, 3);
  srm.extendMZ(1.0, 2); // level 2 already present

  std::set<UInt> levels = srm.getMSLevels();
  TEST_EQUAL(levels.size(), 2)
  TEST_EQUAL(levels.count(2), 1)
  TEST_EQUAL(levels.count(3), 1)

  // global extends (ms_level 0) do not register a level
  srm.extendRT(5.0);
  TEST_EQUAL(srm.getMSLevels().size(), 2)
}
END_SECTION

START_SECTION((const BaseType& byMSLevel(UInt ms_level = 0) const))
{
  SpectrumRangeManager srm;
  srm.extendRT(10.0, 2);
  srm.extendRT(20.0, 2);
  TEST_REAL_SIMILAR(srm.byMSLevel(2).getMinRT(), 10.0)

  // global ranges live in the base class, not in the per-level map:
  // byMSLevel(0) therefore has no entry and throws
  srm.extendRT(100.0); // global
  TEST_EXCEPTION(Exception::InvalidValue, srm.byMSLevel(0))

  // unknown level throws
  TEST_EXCEPTION(Exception::InvalidValue, srm.byMSLevel(99))
}
END_SECTION

START_SECTION((void extend(const BaseType& other, UInt ms_level = 0)))
{
  SpectrumRangeManager srm;
  SpectrumRangeManager::BaseType other;
  other.extendRT(1000.0);
  other.extendRT(2000.0);

  // into a specific level
  srm.extend(other, 5);
  TEST_REAL_SIMILAR(srm.byMSLevel(5).getMinRT(), 1000.0)
  TEST_REAL_SIMILAR(srm.byMSLevel(5).getMaxRT(), 2000.0)

  // into the global range (ms_level 0)
  srm.extendRT(100.0);
  srm.extend(other, 0);
  TEST_REAL_SIMILAR(srm.getMinRT(), 100.0)
  TEST_REAL_SIMILAR(srm.getMaxRT(), 2000.0)
}
END_SECTION

START_SECTION((void extendUnsafe(const MSSpectrum& spectrum, UInt ms_level = 0)))
{
  SpectrumRangeManager srm;
  MSSpectrum s;
  Peak1D p1; p1.setMZ(700.0); p1.setIntensity(5.0); s.push_back(p1);
  Peak1D p2; p2.setMZ(800.0); p2.setIntensity(9.0); s.push_back(p2);
  s.updateRanges();

  srm.extendUnsafe(s, 4);
  TEST_REAL_SIMILAR(srm.byMSLevel(4).getMinMZ(), 700.0)
  TEST_REAL_SIMILAR(srm.byMSLevel(4).getMaxMZ(), 800.0)
}
END_SECTION

START_SECTION((void clearRanges()))
{
  SpectrumRangeManager srm;
  srm.extendRT(100.0);       // global
  srm.extendRT(10.0, 2);     // level 2
  TEST_EQUAL(srm.getMSLevels().empty(), false)

  srm.clearRanges();
  TEST_EQUAL(srm.getMSLevels().empty(), true)
  TEST_EQUAL(srm.getRangeForDim(MSDim::RT).isEmpty(), true)
}
END_SECTION

START_SECTION((SpectrumRangeManager(const SpectrumRangeManager& source)))
{
  SpectrumRangeManager srm;
  srm.extendRT(100.0);       // global
  srm.extendRT(10.0, 2);
  srm.extendRT(20.0, 2);

  SpectrumRangeManager copy(srm);
  TEST_REAL_SIMILAR(copy.getMinRT(), 100.0)           // global copied
  TEST_EQUAL(copy.getMSLevels().count(2), 1)          // level copied
  TEST_REAL_SIMILAR(copy.byMSLevel(2).getMaxRT(), 20.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
