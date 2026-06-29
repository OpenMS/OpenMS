// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <memory>

using namespace OpenMS;
using namespace std;

namespace
{
  std::shared_ptr<MSExperiment> mkExp()
  {
    auto exp = std::make_shared<MSExperiment>();
    MSSpectrum s0; s0.setRT(10.0); s0.setMSLevel(1);
    Peak1D p; p.setMZ(100.0); p.setIntensity(500.0); s0.push_back(p);
    Peak1D p2; p2.setMZ(200.0); p2.setIntensity(800.0); s0.push_back(p2);
    exp->addSpectrum(s0);
    MSSpectrum s1; s1.setRT(20.0); s1.setMSLevel(2); exp->addSpectrum(s1);
    MSChromatogram c; c.setNativeID("chrom_42");
    ChromatogramPeak cp; cp.setRT(5.0); cp.setIntensity(99.0); c.push_back(cp);
    exp->addChromatogram(c);
    return exp;
  }
}

START_TEST(SpectrumAccessOpenMSInMemory, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessOpenMSInMemory* ptr = nullptr;
SpectrumAccessOpenMSInMemory* null_ptr = nullptr;

START_SECTION((explicit SpectrumAccessOpenMSInMemory(OpenSwath::ISpectrumAccess& origin)))
{
  SpectrumAccessOpenMS origin(mkExp());
  ptr = new SpectrumAccessOpenMSInMemory(origin);
  TEST_NOT_EQUAL(ptr, null_ptr)
  // the in-memory copy pre-loads the same spectra/chromatograms
  TEST_EQUAL(ptr->getNrSpectra(), 2)
  TEST_EQUAL(ptr->getNrChromatograms(), 1)
}
END_SECTION

START_SECTION((~SpectrumAccessOpenMSInMemory()))
{
  delete ptr;
}
END_SECTION

START_SECTION((OpenSwath::SpectrumPtr getSpectrumById(int id)))
{
  SpectrumAccessOpenMS origin(mkExp());
  SpectrumAccessOpenMSInMemory sa(origin);
  OpenSwath::SpectrumPtr s = sa.getSpectrumById(0);
  TEST_EQUAL(s->getMZArray()->data.size(), 2)
  TEST_REAL_SIMILAR(s->getMZArray()->data[0], 100.0)
  TEST_REAL_SIMILAR(s->getMZArray()->data[1], 200.0)
  TEST_REAL_SIMILAR(s->getIntensityArray()->data[1], 800.0)
}
END_SECTION

START_SECTION((OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const))
{
  SpectrumAccessOpenMS origin(mkExp());
  SpectrumAccessOpenMSInMemory sa(origin);
  OpenSwath::SpectrumMeta m = sa.getSpectrumMetaById(0);
  TEST_REAL_SIMILAR(m.RT, 10.0)
}
END_SECTION

START_SECTION((OpenSwath::ChromatogramPtr getChromatogramById(int id)))
{
  SpectrumAccessOpenMS origin(mkExp());
  SpectrumAccessOpenMSInMemory sa(origin);
  OpenSwath::ChromatogramPtr c = sa.getChromatogramById(0);
  TEST_EQUAL(c->getTimeArray()->data.size(), 1)
  TEST_REAL_SIMILAR(c->getTimeArray()->data[0], 5.0)
  TEST_REAL_SIMILAR(c->getIntensityArray()->data[0], 99.0)
}
END_SECTION

START_SECTION((std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const))
{
  SpectrumAccessOpenMS origin(mkExp());
  SpectrumAccessOpenMSInMemory sa(origin);
  auto clone = sa.lightClone();
  TEST_NOT_EQUAL(clone.get(), nullptr)
  TEST_EQUAL(clone->getNrSpectra(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
