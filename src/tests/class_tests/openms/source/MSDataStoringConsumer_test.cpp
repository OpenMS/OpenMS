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

#include <OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>

START_TEST(MSDataStoringConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataStoringConsumer* storing_consumer_ptr = nullptr;
MSDataStoringConsumer* storing_consumer_nullPointer = nullptr;

START_SECTION((MSDataStoringConsumer()))
  storing_consumer_ptr = new MSDataStoringConsumer();
  TEST_NOT_EQUAL(storing_consumer_ptr, storing_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataStoringConsumer()))
    delete storing_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  MSDataStoringConsumer * storing_consumer = new MSDataStoringConsumer();

  MSSpectrum s;
  s.setName("spec1");
  s.setComment("comm1");
  s.setRT(5);
  storing_consumer->consumeSpectrum(s);
  s.setName("spec2");
  s.setComment("comm2");
  s.setRT(15);
  storing_consumer->consumeSpectrum(s);
  s.setName("spec3");
  s.setComment("comm3");
  s.setRT(25);
  storing_consumer->consumeSpectrum(s);

  TEST_EQUAL(storing_consumer->getData().getNrSpectra(), 3);
  TEST_EQUAL(storing_consumer->getData().getNrChromatograms(), 0)


  TEST_EQUAL(storing_consumer->getData().getSpectra()[0].getName(), "spec1")
  TEST_EQUAL(storing_consumer->getData().getSpectra()[1].getName(), "spec2")
  TEST_EQUAL(storing_consumer->getData().getSpectra()[2].getName(), "spec3")

  TEST_EQUAL(storing_consumer->getData().getSpectra()[0].getComment(), "comm1")
  TEST_EQUAL(storing_consumer->getData().getSpectra()[1].getComment(), "comm2")
  TEST_EQUAL(storing_consumer->getData().getSpectra()[2].getComment(), "comm3")

  delete storing_consumer;
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  MSDataStoringConsumer * storing_consumer = new MSDataStoringConsumer();

  MSChromatogram c;
  c.setNativeID("testid");
  storing_consumer->consumeChromatogram(c);

  TEST_EQUAL(storing_consumer->getData().getNrSpectra(), 0)
  TEST_EQUAL(storing_consumer->getData().getNrChromatograms(), 1)
  TEST_EQUAL(storing_consumer->getData().getChromatograms()[0].getNativeID(), "testid")

  delete storing_consumer;
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  MSDataStoringConsumer * storing_consumer = new MSDataStoringConsumer();
  storing_consumer->setExpectedSize(1,1);

  MSChromatogram c;
  c.setNativeID("testid");
  storing_consumer->consumeChromatogram(c);

  MSSpectrum spec;
  spec.setName("spec1");
  spec.setRT(5);
  storing_consumer->consumeSpectrum(spec);

  ExperimentalSettings s;
  s.setComment("mySettings");
  storing_consumer->setExperimentalSettings(s);

  TEST_NOT_EQUAL(storing_consumer, storing_consumer_nullPointer)

  TEST_EQUAL(storing_consumer->getData().getNrSpectra(), 1)
  TEST_EQUAL(storing_consumer->getData().getNrChromatograms(), 1)
  TEST_EQUAL(storing_consumer->getData().getComment(), "mySettings")

  delete storing_consumer;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
