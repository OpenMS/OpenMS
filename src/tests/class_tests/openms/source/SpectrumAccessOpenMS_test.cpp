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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>

#include <memory>

using namespace OpenMS;
using namespace std;

namespace
{
  // a 2-spectrum, 1-chromatogram experiment
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

START_TEST(SpectrumAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessOpenMS* ptr = nullptr;
SpectrumAccessOpenMS* null_ptr = nullptr;

START_SECTION((explicit SpectrumAccessOpenMS(std::shared_ptr<MSExperimentType> ms_experiment)))
{
  ptr = new SpectrumAccessOpenMS(mkExp());
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((~SpectrumAccessOpenMS()))
{
  delete ptr;
}
END_SECTION

START_SECTION((size_t getNrSpectra() const))
{
  SpectrumAccessOpenMS sa(mkExp());
  TEST_EQUAL(sa.getNrSpectra(), 2)
}
END_SECTION

START_SECTION((size_t getNrChromatograms() const))
{
  SpectrumAccessOpenMS sa(mkExp());
  TEST_EQUAL(sa.getNrChromatograms(), 1)
}
END_SECTION

START_SECTION((OpenSwath::SpectrumPtr getSpectrumById(int id)))
{
  SpectrumAccessOpenMS sa(mkExp());
  OpenSwath::SpectrumPtr s = sa.getSpectrumById(0);
  TEST_EQUAL(s->getMZArray()->data.size(), 2)
  TEST_REAL_SIMILAR(s->getMZArray()->data[0], 100.0)
  TEST_REAL_SIMILAR(s->getMZArray()->data[1], 200.0)
  TEST_REAL_SIMILAR(s->getIntensityArray()->data[0], 500.0)
  TEST_REAL_SIMILAR(s->getIntensityArray()->data[1], 800.0)
}
END_SECTION

START_SECTION((OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const))
{
  SpectrumAccessOpenMS sa(mkExp());
  OpenSwath::SpectrumMeta m = sa.getSpectrumMetaById(0);
  TEST_REAL_SIMILAR(m.RT, 10.0)
  TEST_EQUAL(m.ms_level, 1)
}
END_SECTION

START_SECTION((std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const))
{
  SpectrumAccessOpenMS sa(mkExp());
  // exact RT, zero tolerance -> only spectrum 0 (RT 10)
  auto hits = sa.getSpectraByRT(10.0, 0.0);
  TEST_EQUAL(hits.size(), 1)
  TEST_EQUAL(hits[0], 0)
  // window covering both spectra (RT 10 and 20)
  TEST_EQUAL(sa.getSpectraByRT(15.0, 6.0).size(), 2)
}
END_SECTION

START_SECTION((OpenSwath::ChromatogramPtr getChromatogramById(int id)))
{
  SpectrumAccessOpenMS sa(mkExp());
  OpenSwath::ChromatogramPtr c = sa.getChromatogramById(0);
  TEST_EQUAL(c->getTimeArray()->data.size(), 1)
  TEST_REAL_SIMILAR(c->getTimeArray()->data[0], 5.0)
  TEST_REAL_SIMILAR(c->getIntensityArray()->data[0], 99.0)
}
END_SECTION

START_SECTION((std::string getChromatogramNativeID(int id) const))
{
  SpectrumAccessOpenMS sa(mkExp());
  TEST_STRING_EQUAL(sa.getChromatogramNativeID(0), "chrom_42")
}
END_SECTION

START_SECTION((std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const))
{
  SpectrumAccessOpenMS sa(mkExp());
  auto clone = sa.lightClone();
  TEST_NOT_EQUAL(clone.get(), nullptr)
  TEST_EQUAL(clone->getNrSpectra(), 2)
  TEST_EQUAL(clone->getNrChromatograms(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
