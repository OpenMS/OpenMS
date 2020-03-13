// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


boost::shared_ptr<PeakMap > getData()
{
  boost::shared_ptr<PeakMap > exp2(new PeakMap);
  MSSpectrum spec;
  Peak1D p;
  p.setMZ(100);
  p.setIntensity(50);
  spec.push_back(p);
  p.setMZ(500);
  p.setIntensity(150);
  spec.push_back(p);
  exp2->addSpectrum(spec);
  return exp2;
}

START_TEST(SpectrumAccessQuadMZTransforming, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumAccessQuadMZTransforming* ptr = nullptr;
SpectrumAccessQuadMZTransforming* nullPointer = nullptr;

boost::shared_ptr<PeakMap > exp(new PeakMap);
OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);


START_SECTION(SpectrumAccessQuadMZTransforming())
{
  ptr = new SpectrumAccessQuadMZTransforming(expptr, 0, 0, 0, false);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SpectrumAccessQuadMZTransforming())
{
  delete ptr;
}
END_SECTION


START_SECTION(size_t getNrSpectra() const)
{
  boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr(new SpectrumAccessQuadMZTransforming(expptr, 0, 0, 0, false));
  TEST_EQUAL(ptr->getNrSpectra(), 0)

  boost::shared_ptr<PeakMap > exp2 = getData();
  OpenSwath::SpectrumAccessPtr expptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);
  boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr2(new SpectrumAccessQuadMZTransforming(expptr2, 0, 0, 0, false));
  TEST_EQUAL(ptr2->getNrSpectra(), 1)
}
END_SECTION

START_SECTION(OpenSwath::SpectrumPtr getSpectrumById(int id))
{
  {
    boost::shared_ptr<PeakMap > exp2 = getData();
    OpenSwath::SpectrumAccessPtr expptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);
    boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr2(new SpectrumAccessQuadMZTransforming(expptr2, 0, 0, 0, false));
    OpenSwath::SpectrumPtr spec1 = ptr2->getSpectrumById(0);
    TEST_NOT_EQUAL(spec1.get(), nullPointer) // pointer is present
    TEST_EQUAL(bool(spec1), true) // pointer is not null
    TEST_EQUAL(spec1->getMZArray()->data.size(), 2) 
    TEST_EQUAL(spec1->getIntensityArray()->data.size(), 2) 

    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[0], 50) 
    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[1], 150) 
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[0], 0)
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[1], 0)
  }

  {
    boost::shared_ptr<PeakMap > exp2 = getData();
    OpenSwath::SpectrumAccessPtr expptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);
    boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr2(new SpectrumAccessQuadMZTransforming(expptr2, 10, 5, 2, false));
    OpenSwath::SpectrumPtr spec1 = ptr2->getSpectrumById(0);
    TEST_NOT_EQUAL(spec1.get(), nullPointer) // pointer is present
    TEST_EQUAL(bool(spec1), true) // pointer is not null
    TEST_EQUAL(spec1->getMZArray()->data.size(), 2) 
    TEST_EQUAL(spec1->getIntensityArray()->data.size(), 2) 

    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[0], 50) 
    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[1], 150) 
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[0], 10 + 100*5 + 100*100* 2)
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[1], 10 + 500*5 + 500*500* 2)
  }

}
END_SECTION

START_SECTION(boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const)
{
  boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr(new SpectrumAccessQuadMZTransforming(expptr, 0, 0, 0, false));
  boost::shared_ptr<OpenSwath::ISpectrumAccess> clone_ptr_empty = ptr->lightClone();

  TEST_EQUAL(ptr->getNrSpectra(), clone_ptr_empty->getNrSpectra())

  {
    boost::shared_ptr<PeakMap > exp2 = getData();
    OpenSwath::SpectrumAccessPtr expptr2 = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp2);
    boost::shared_ptr<SpectrumAccessQuadMZTransforming> ptr2(new SpectrumAccessQuadMZTransforming(expptr2, 10, 5, 2, false));
    OpenSwath::SpectrumPtr spec1 = ptr2->getSpectrumById(0);
    TEST_NOT_EQUAL(spec1.get(), nullPointer) // pointer is present
    TEST_EQUAL(bool(spec1), true) // pointer is not null
    TEST_EQUAL(spec1->getMZArray()->data.size(), 2) 
    TEST_EQUAL(spec1->getIntensityArray()->data.size(), 2) 

    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[0], 50) 
    TEST_REAL_SIMILAR(spec1->getIntensityArray()->data[1], 150) 
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[0], 10 + 100*5 + 100*100* 2)
    TEST_REAL_SIMILAR(spec1->getMZArray()->data[1], 10 + 500*5 + 500*500* 2)

    boost::shared_ptr<OpenSwath::ISpectrumAccess> clone_ptr = ptr2->lightClone();
    TEST_EQUAL(ptr2->getNrSpectra(), clone_ptr->getNrSpectra())

    OpenSwath::SpectrumPtr spec_clone = ptr2->getSpectrumById(0);
    TEST_NOT_EQUAL(spec_clone.get(), nullPointer) // pointer is present
    TEST_EQUAL(bool(spec_clone), true) // pointer is not null
    TEST_EQUAL(spec_clone->getMZArray()->data.size(), 2) 
    TEST_EQUAL(spec_clone->getIntensityArray()->data.size(), 2) 

    TEST_REAL_SIMILAR(spec_clone->getIntensityArray()->data[0], 50) 
    TEST_REAL_SIMILAR(spec_clone->getIntensityArray()->data[1], 150) 
    TEST_REAL_SIMILAR(spec_clone->getMZArray()->data[0], 10 + 100*5 + 100*100* 2)
    TEST_REAL_SIMILAR(spec_clone->getMZArray()->data[1], 10 + 500*5 + 500*500* 2)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

