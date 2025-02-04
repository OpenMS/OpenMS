// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <boost/shared_ptr.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessOpenMS* ptr = nullptr;
SpectrumAccessOpenMS* nullPointer = nullptr;

START_SECTION(SpectrumAccessOpenMS())
{
  boost::shared_ptr< PeakMap > exp ( new PeakMap );
  ptr = new SpectrumAccessOpenMS(exp);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SpectrumAccessOpenMS())
{
  delete ptr;
}
END_SECTION

START_SECTION( size_t getNrSpectra() const)
{
  {
    boost::shared_ptr< PeakMap > exp ( new PeakMap );
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 0);
    TEST_EQUAL(spectrum_acc.getNrChromatograms(), 0);
  }

  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;
    MSChromatogram c;
    new_exp->addSpectrum(s);
    new_exp->addSpectrum(s);
    new_exp->addChromatogram(c);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 2);
    TEST_EQUAL(spectrum_acc.getNrChromatograms(), 1);
  }
}
END_SECTION

START_SECTION ( boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const)
{
  PeakMap* new_exp = new PeakMap;
  MSSpectrum s;
  MSChromatogram c;
  new_exp->addSpectrum(s);
  new_exp->addSpectrum(s);
  new_exp->addChromatogram(c);
  boost::shared_ptr< PeakMap > exp (new_exp);
  SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

  TEST_EQUAL(spectrum_acc.getNrSpectra(), 2);
  TEST_EQUAL(spectrum_acc.getNrChromatograms(), 1);

  boost::shared_ptr<OpenSwath::ISpectrumAccess> sa_clone = spectrum_acc.lightClone();
  TEST_EQUAL(sa_clone->getNrSpectra(), 2);
  TEST_EQUAL(sa_clone->getNrChromatograms(), 1);
}
END_SECTION

START_SECTION ( OpenSwath::SpectrumPtr getSpectrumById(int id);)
{
  {
    boost::shared_ptr< PeakMap > exp ( new PeakMap );
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 0);
    // OpenSwath::SpectrumPtr sptr = spectrum_acc.getSpectrumById(0);
  }

  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;

    s.setRT(20);
    Peak1D p;
    p.setMZ(20.0);
    s.push_back(p);

    new_exp->addSpectrum(s);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 1);
    OpenSwath::SpectrumPtr sptr = spectrum_acc.getSpectrumById(0);
    TEST_REAL_SIMILAR (sptr->getMZArray()->data[0], 20.0);
  }

  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;

    s.setRT(20);
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(22.0);
    s.push_back(p);

    OpenMS::DataArrays::FloatDataArray fda;
    fda.push_back(50);
    fda.setName("testName");
    auto fdas = s.getFloatDataArrays();
    fdas.push_back(fda);
    s.setFloatDataArrays(fdas);

    OpenMS::DataArrays::IntegerDataArray ida;
    ida.push_back(51);
    ida.setName("testName_integer");
    auto idas = s.getIntegerDataArrays();
    idas.push_back(ida);
    s.setIntegerDataArrays(idas);

    new_exp->addSpectrum(s);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 1)
    OpenSwath::SpectrumPtr sptr = spectrum_acc.getSpectrumById(0);
    TEST_REAL_SIMILAR (sptr->getMZArray()->data[0], 20.0)
    TEST_REAL_SIMILAR (sptr->getIntensityArray()->data[0], 22.0)
    TEST_EQUAL (sptr->getDataArrays().size(), 4)
    TEST_EQUAL (sptr->getDataArrays()[2]->description, "testName")
    TEST_REAL_SIMILAR (sptr->getDataArrays()[2]->data[0], 50.0)
    TEST_EQUAL (sptr->getDataArrays()[3]->description, "testName_integer")
    TEST_REAL_SIMILAR (sptr->getDataArrays()[3]->data[0], 51.0)
  }
}
END_SECTION

START_SECTION ( OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const)
{
  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;
    s.setRT(20);
    new_exp->addSpectrum(s);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 1);
    OpenSwath::SpectrumMeta spmeta = spectrum_acc.getSpectrumMetaById(0);
    TEST_REAL_SIMILAR(spmeta.RT, 20.0);
  }
}
END_SECTION

START_SECTION ( SpectrumSettings getSpectraMetaInfo(int id) const)
{
  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;
    s.setComment("remember me");
    new_exp->addSpectrum(s);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 1);
    SpectrumSettings spmeta = spectrum_acc.getSpectraMetaInfo(0);
    TEST_EQUAL (spmeta.getComment(), "remember me");
  }
}
END_SECTION

START_SECTION ( std::vector<std::size_t> SpectrumAccessOpenMS::getSpectraByRT(double RT, double deltaRT) const)
{
  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;
    MSChromatogram c;
    s.setRT(20);
    new_exp->addSpectrum(s);
    s.setRT(40);
    new_exp->addSpectrum(s);
    new_exp->addChromatogram(c);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 2);
    TEST_EQUAL(spectrum_acc.getNrChromatograms(), 1);

    TEST_EQUAL(spectrum_acc.getSpectraByRT(20, 5.0).size(),  1);
    TEST_EQUAL(spectrum_acc.getSpectraByRT(20, 25.0).size(), 2);
    TEST_EQUAL(spectrum_acc.getSpectraByRT(40, 5.0).size(),  1);
    TEST_EQUAL(spectrum_acc.getSpectraByRT(40, 25.0).size(), 2);
    TEST_EQUAL(spectrum_acc.getSpectraByRT(50, 5.0).size(),  0);
  }

}
END_SECTION

START_SECTION( size_t getNrChromatograms() const)
{
  NOT_TESTABLE // see getNrSpectra
}
END_SECTION

START_SECTION(OpenSwath::ChromatogramPtr getChromatogramById(int id))
{
  {
    boost::shared_ptr< PeakMap > exp ( new PeakMap );
    SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(spectrum_acc.getNrSpectra(), 0);
    TEST_EQUAL(spectrum_acc.getNrChromatograms(), 0);
  }

  {
    PeakMap* new_exp = new PeakMap;
    MSSpectrum s;
    MSChromatogram c;

    c.setName("chrom_nr_1");
    c.setNativeID("native_id_nr_1");
    ChromatogramPeak p;
    p.setRT(20.0);
    c.push_back(p);

    new_exp->addSpectrum(s);
    new_exp->addSpectrum(s);
    new_exp->addChromatogram(c);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS chrom_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(chrom_acc.getNrSpectra(), 2);
    TEST_EQUAL(chrom_acc.getNrChromatograms(), 1);

    OpenSwath::ChromatogramPtr cptr = chrom_acc.getChromatogramById(0);
    TEST_REAL_SIMILAR (cptr->getTimeArray()->data[0], 20.0);
  }

  {
    PeakMap* new_exp = new PeakMap;
    MSChromatogram chrom;

    ChromatogramPeak p;
    p.setRT(20.0);
    p.setIntensity(22.0);
    chrom.push_back(p);

    OpenMS::DataArrays::FloatDataArray fda;
    fda.push_back(50);
    fda.setName("testName");
    auto fdas = chrom.getFloatDataArrays();
    fdas.push_back(fda);
    chrom.setFloatDataArrays(fdas);

    OpenMS::DataArrays::IntegerDataArray ida;
    ida.push_back(51);
    ida.setName("testName_integer");
    auto idas = chrom.getIntegerDataArrays();
    idas.push_back(ida);
    chrom.setIntegerDataArrays(idas);

    new_exp->addChromatogram(chrom);
    boost::shared_ptr< PeakMap > exp (new_exp);
    SpectrumAccessOpenMS chrom_acc = SpectrumAccessOpenMS(exp);

    TEST_EQUAL(chrom_acc.getNrChromatograms(), 1)
    OpenSwath::ChromatogramPtr cptr = chrom_acc.getChromatogramById(0);
    TEST_REAL_SIMILAR (cptr->getTimeArray()->data[0], 20.0)
    TEST_REAL_SIMILAR (cptr->getIntensityArray()->data[0], 22.0)
    TEST_EQUAL (cptr->getDataArrays().size(), 4)
    TEST_EQUAL (cptr->getDataArrays()[2]->description, "testName")
    TEST_REAL_SIMILAR (cptr->getDataArrays()[2]->data[0], 50.0)
    TEST_EQUAL (cptr->getDataArrays()[3]->description, "testName_integer")
    TEST_REAL_SIMILAR (cptr->getDataArrays()[3]->data[0], 51.0)
  }
}
END_SECTION

START_SECTION(std::string getChromatogramNativeID(int id) const)
{
  PeakMap* new_exp = new PeakMap;
  MSSpectrum s;
  MSChromatogram c;

  c.setName("chrom_nr_1");
  c.setNativeID("native_id_nr_1");
  ChromatogramPeak p;
  p.setRT(20.0);
  c.push_back(p);

  new_exp->addSpectrum(s);
  new_exp->addSpectrum(s);
  new_exp->addChromatogram(c);
  boost::shared_ptr< PeakMap > exp (new_exp);
  SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

  OpenSwath::ChromatogramPtr cptr = spectrum_acc.getChromatogramById(0);
  TEST_EQUAL (spectrum_acc.getChromatogramNativeID(0), "native_id_nr_1")
}
END_SECTION

START_SECTION (ChromatogramSettings getChromatogramMetaInfo(int id) const)
{
  PeakMap* new_exp = new PeakMap;
  MSChromatogram c;
  c.setComment("remember me");
  new_exp->addChromatogram(c);
  boost::shared_ptr< PeakMap > exp (new_exp);
  SpectrumAccessOpenMS spectrum_acc = SpectrumAccessOpenMS(exp);

  TEST_EQUAL(spectrum_acc.getNrChromatograms(), 1);
  ChromatogramSettings cpmeta = spectrum_acc.getChromatogramMetaInfo(0);
  TEST_EQUAL (cpmeta.getComment(), "remember me");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



