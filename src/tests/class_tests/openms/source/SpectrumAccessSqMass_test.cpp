// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

SpectrumAccessSqMass* ptr = 0;
SpectrumAccessSqMass* nullPointer = 0;

boost::shared_ptr<PeakMap > exp(new PeakMap);
OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);


START_SECTION(SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  ptr = new SpectrumAccessSqMass(handler);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(SpectrumAccessSqMass(OpenMS::Internal::MzMLSqliteHandler handler, std::vector<int> indices))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  std::vector<int> indices;
  indices.push_back(1);
  ptr = new SpectrumAccessSqMass(handler, indices);
  TEST_NOT_EQUAL(ptr, nullPointer)

  TEST_EQUAL(ptr->getNrSpectra(), 1)
}
END_SECTION

START_SECTION(SpectrumAccessSqMass(SpectrumAccessSqMass sp, std::vector<int> indices))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  ptr = new SpectrumAccessSqMass(handler);
  TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EQUAL(ptr->getNrSpectra(), 2)

  // select subset of the data (all two spectra)
  ptr = new SpectrumAccessSqMass(handler);
  {
    std::vector<int> indices;
    indices.push_back(0);
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(*ptr, indices);
    TEST_NOT_EQUAL(ptr, nullPointer)
    TEST_EQUAL(ptr->getNrSpectra(), 2)
  }

  // select subset of the data (only second spectrum)
  ptr = new SpectrumAccessSqMass(handler);
  {
    std::vector<int> indices;
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(*ptr, indices);
    TEST_NOT_EQUAL(ptr, nullPointer)
    TEST_EQUAL(ptr->getNrSpectra(), 1)
  }

  // this should not work, ptr has now only a single spectrum
  std::vector<int> indices;
  indices.push_back(1);
  ptr = new SpectrumAccessSqMass(*ptr, indices);

}
END_SECTION

START_SECTION(~SpectrumAccessSqMass())
{
  delete ptr;
}
END_SECTION

START_SECTION(size_t getNrSpectra() const)
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  ptr = new SpectrumAccessSqMass(handler);
  TEST_EQUAL(ptr->getNrSpectra(), 2)
}
END_SECTION

START_SECTION(boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const)
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

  ptr = new SpectrumAccessSqMass(handler);
  TEST_EQUAL(ptr->getNrSpectra(), 2)

  boost::shared_ptr<OpenSwath::ISpectrumAccess> ptr2 = ptr->lightClone();
  TEST_EQUAL(ptr2->getNrSpectra(), 2)
}
END_SECTION

START_SECTION(void getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"));

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
  }

  // select only 2nd spectrum iteratively
  {

    std::vector<int> indices;
    indices.push_back(1);

    ptr = new SpectrumAccessSqMass(handler, indices);
    TEST_EQUAL(ptr->getNrSpectra(), 1)

    // now we have an interface with a single spectrum in it, so if we select
    // the first spectrum of THAT interface, it should be the 2nd spectrum from
    // the initial dataset
    indices.clear();
    // indices.push_back(1); // this should not work as we now have only a single spectrum (out of bounds access!)
    indices.push_back(0);

    ptr = new SpectrumAccessSqMass(*ptr, indices);
    TEST_EQUAL(ptr->getNrSpectra(), 1)

    std::vector< OpenSwath::SpectrumPtr > spectra;
    std::vector< OpenSwath::SpectrumMeta > spectra_meta;
    ptr->getAllSpectra(spectra, spectra_meta);

    TEST_EQUAL(spectra.size(), 1)
    TEST_EQUAL(spectra_meta.size(), 1)

    TEST_EQUAL(spectra[0]->getMZArray()->data.size(), 19800)
    TEST_EQUAL(spectra[0]->getIntensityArray()->data.size(), 19800)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

