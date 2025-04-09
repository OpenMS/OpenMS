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
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>
///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumAccessSqMass, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessSqMass* ptr = nullptr;
SpectrumAccessSqMass* nullPointer = nullptr;

boost::shared_ptr<PeakMap > exp(new PeakMap);
OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);


START_SECTION(SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  ptr = new SpectrumAccessSqMass(handler);
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
}
END_SECTION

    

START_SECTION(SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int> & indices))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  std::vector<int> indices;
  indices.push_back(1);
  ptr = new SpectrumAccessSqMass(handler, indices);
  TEST_NOT_EQUAL(ptr, nullPointer)

  TEST_EQUAL(ptr->getNrSpectra(), 1)
  delete ptr;
}
END_SECTION

START_SECTION(SpectrumAccessSqMass(const SpectrumAccessSqMass& sp, const std::vector<int>& indices))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  ptr = new SpectrumAccessSqMass(handler);
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->getNrSpectra(), 2)
  delete ptr;

  SpectrumAccessSqMass sasm(handler);

  // select subset of the data (all two spectra)
  {
    std::vector<int> indices;
    indices.push_back(0);
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(sasm, indices);
    TEST_NOT_EQUAL(ptr, nullPointer)
    TEST_EQUAL(ptr->getNrSpectra(), 2)
    delete ptr;
  }

  // select subset of the data (only second spectrum)
  {
    std::vector<int> indices;
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(sasm, indices);
    TEST_NOT_EQUAL(ptr, nullPointer)
    TEST_EQUAL(ptr->getNrSpectra(), 1)
  }

  // this should not work:
  // selecting subset of data (only second spectrum) shouldn't work if there is only a single spectrum
  {
    std::vector<int> indices;
    indices.push_back(1);
    TEST_EXCEPTION(Exception::IllegalArgument, new SpectrumAccessSqMass(*ptr, indices))

    std::vector<int> indices2;
    indices2.push_back(50);
    TEST_EXCEPTION(Exception::IllegalArgument, new SpectrumAccessSqMass(*ptr, indices))
  }
}
END_SECTION

START_SECTION(~SpectrumAccessSqMass())
{
  delete ptr;
}
END_SECTION

START_SECTION(size_t getNrSpectra() const)
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  ptr = new SpectrumAccessSqMass(handler);
  TEST_EQUAL(ptr->getNrSpectra(), 2)
  delete ptr;
}
END_SECTION

START_SECTION(boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const)
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  ptr = new SpectrumAccessSqMass(handler);
  TEST_EQUAL(ptr->getNrSpectra(), 2)

  boost::shared_ptr<OpenSwath::ISpectrumAccess> ptr2 = ptr->lightClone();
  TEST_EQUAL(ptr2->getNrSpectra(), 2)
  delete ptr;
}
END_SECTION

START_SECTION(void getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  {
    ptr = new SpectrumAccessSqMass(handler);
    TEST_EQUAL(ptr->getNrSpectra(), 2)

    std::vector< OpenSwath::SpectrumPtr > spectra;
    std::vector< OpenSwath::SpectrumMeta > spectra_meta;
    ptr->getAllSpectra(spectra, spectra_meta);

    TEST_EQUAL(spectra.size(), 2)
    TEST_EQUAL(spectra_meta.size(), 2)

    TEST_EQUAL(spectra[0]->getMZArray()->data.size(), 19914)
    TEST_EQUAL(spectra[0]->getIntensityArray()->data.size(), 19914)

    TEST_EQUAL(spectra[1]->getMZArray()->data.size(), 19800)
    TEST_EQUAL(spectra[1]->getIntensityArray()->data.size(), 19800)
    delete ptr;
  }

  {
    std::vector<int> indices;
    indices.push_back(0);
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(handler, indices);
    TEST_EQUAL(ptr->getNrSpectra(), 2)

    std::vector< OpenSwath::SpectrumPtr > spectra;
    std::vector< OpenSwath::SpectrumMeta > spectra_meta;
    ptr->getAllSpectra(spectra, spectra_meta);

    TEST_EQUAL(spectra.size(), 2)
    TEST_EQUAL(spectra_meta.size(), 2)

    TEST_EQUAL(spectra[0]->getMZArray()->data.size(), 19914)
    TEST_EQUAL(spectra[0]->getIntensityArray()->data.size(), 19914)

    TEST_EQUAL(spectra[1]->getMZArray()->data.size(), 19800)
    TEST_EQUAL(spectra[1]->getIntensityArray()->data.size(), 19800)
    delete ptr;
  }

  // select only 2nd spectrum 
  {
    std::vector<int> indices;
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(handler, indices);
    TEST_EQUAL(ptr->getNrSpectra(), 1)

    std::vector< OpenSwath::SpectrumPtr > spectra;
    std::vector< OpenSwath::SpectrumMeta > spectra_meta;
    ptr->getAllSpectra(spectra, spectra_meta);

    TEST_EQUAL(spectra.size(), 1)
    TEST_EQUAL(spectra_meta.size(), 1)

    TEST_EQUAL(spectra[0]->getMZArray()->data.size(), 19800)
    TEST_EQUAL(spectra[0]->getIntensityArray()->data.size(), 19800)
    delete ptr;
  }

  // select only 2nd spectrum iteratively
  {

    std::vector<int> indices;
    indices.push_back(1);

    SpectrumAccessSqMass* sasm = new SpectrumAccessSqMass(handler, indices);
    TEST_EQUAL(sasm->getNrSpectra(), 1)

    // now we have an interface with a single spectrum in it, so if we select
    // the first spectrum of THAT interface, it should be the 2nd spectrum from
    // the initial dataset
    indices.clear();
    // indices.push_back(1); // this should not work as we now have only a single spectrum (out of bounds access!)
    indices.push_back(0);

    ptr = new SpectrumAccessSqMass(*sasm, indices);
    TEST_EQUAL(ptr->getNrSpectra(), 1)

    std::vector< OpenSwath::SpectrumPtr > spectra;
    std::vector< OpenSwath::SpectrumMeta > spectra_meta;
    ptr->getAllSpectra(spectra, spectra_meta);

    TEST_EQUAL(spectra.size(), 1)
    TEST_EQUAL(spectra_meta.size(), 1)

    TEST_EQUAL(spectra[0]->getMZArray()->data.size(), 19800)
    TEST_EQUAL(spectra[0]->getIntensityArray()->data.size(), 19800)
    delete ptr;
    delete sasm;
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

