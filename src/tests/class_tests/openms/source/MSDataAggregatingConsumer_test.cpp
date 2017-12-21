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

#include <OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>

///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataStoringConsumer.h>

START_TEST(MSDataAggregatingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataAggregatingConsumer* agg_consumer_ptr = nullptr;
MSDataAggregatingConsumer* agg_consumer_nullPointer = nullptr;

START_SECTION((MSDataAggregatingConsumer()))
  agg_consumer_ptr = new MSDataAggregatingConsumer(agg_consumer_nullPointer); // dont do that ...
  TEST_NOT_EQUAL(agg_consumer_ptr, agg_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataAggregatingConsumer()))
    delete agg_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  // no adding up
  {
    MSDataStoringConsumer * storage = new MSDataStoringConsumer();
    MSDataAggregatingConsumer * agg_consumer = new MSDataAggregatingConsumer(storage);

    MSSpectrum s;
    s.setName("spec1");
    s.setRT(5);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec2");
    s.setRT(15);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec3");
    s.setRT(25);
    agg_consumer->consumeSpectrum(s);

    // note how we can / have to destroy the aggregate consumer to ensure it
    // flushes the data. The storage object will still be around.
    delete agg_consumer;

    TEST_EQUAL(storage->getData().getNrSpectra(), 3)
    TEST_EQUAL(storage->getData().getNrChromatograms(), 0)

    TEST_EQUAL(storage->getData().getSpectra()[0].getName(), "spec1")
    TEST_EQUAL(storage->getData().getSpectra()[1].getName(), "spec2")
    TEST_EQUAL(storage->getData().getSpectra()[2].getName(), "spec3")

    delete storage;
  }

  // adding empty spectra
  {
    MSDataStoringConsumer * storage = new MSDataStoringConsumer();
    MSDataAggregatingConsumer * agg_consumer = new MSDataAggregatingConsumer(storage);

    MSSpectrum s;
    s.setName("spec1");
    s.setComment("comm1");
    s.setRT(5);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec2");
    s.setComment("comm2");
    s.setRT(5);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec3");
    s.setComment("comm3");
    s.setRT(25);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec4");
    s.setComment("comm4");
    s.setRT(25);
    agg_consumer->consumeSpectrum(s);
    s.setName("spec5");
    s.setComment("comm5");
    s.setRT(35);
    agg_consumer->consumeSpectrum(s);

    // note how we can / have to destroy the aggregate consumer to ensure it
    // flushes the data. The storage object will still be around.
    delete agg_consumer;

    TEST_EQUAL(storage->getData().getNrSpectra(), 3)
    TEST_EQUAL(storage->getData().getNrChromatograms(), 0)

    TEST_EQUAL(storage->getData().getSpectra()[0].getName(), "spec1")
    TEST_EQUAL(storage->getData().getSpectra()[1].getName(), "spec3")
    TEST_EQUAL(storage->getData().getSpectra()[2].getName(), "spec5")

    TEST_EQUAL(storage->getData().getSpectra()[0].getComment(), "comm1")
    TEST_EQUAL(storage->getData().getSpectra()[1].getComment(), "comm3")
    TEST_EQUAL(storage->getData().getSpectra()[2].getComment(), "comm5")

    delete storage;
  }

  // adding full spectra
  {
    MSDataStoringConsumer * storage = new MSDataStoringConsumer();
    MSDataAggregatingConsumer * agg_consumer = new MSDataAggregatingConsumer(storage);

    MSSpectrum s;
    s.setName("spec1");
    s.setComment("comm1");
    s.setRT(5);
    s.push_back(Peak1D(5, 7));
    s.push_back(Peak1D(10, 20));
    s.push_back(Peak1D(15, 30));
    agg_consumer->consumeSpectrum(s);
    s.clear(true);

    s.setName("spec2");
    s.setComment("comm2");
    s.setRT(5);
    s.push_back(Peak1D(5, 10));
    s.push_back(Peak1D(10, 100));
    s.push_back(Peak1D(15, 200));
    agg_consumer->consumeSpectrum(s);
    s.clear(true);

    s.setName("spec3");
    s.setComment("comm3");
    s.setRT(25);
    agg_consumer->consumeSpectrum(s);
    s.clear(true);

    s.setName("spec4");
    s.setComment("comm4");
    s.setRT(25);
    agg_consumer->consumeSpectrum(s);
    s.clear(true);

    s.setName("spec5");
    s.setComment("comm5");
    s.setRT(35);
    agg_consumer->consumeSpectrum(s);
    s.clear(true);

    // note how we can / have to destroy the aggregate consumer to ensure it
    // flushes the data. The storage object will still be around.
    delete agg_consumer;

    TEST_EQUAL(storage->getData().getNrSpectra(), 3)
    TEST_EQUAL(storage->getData().getNrChromatograms(), 0)

    TEST_EQUAL(storage->getData().getSpectra()[0].getName(), "spec1")
    TEST_EQUAL(storage->getData().getSpectra()[1].getName(), "spec3")
    TEST_EQUAL(storage->getData().getSpectra()[2].getName(), "spec5")

    TEST_EQUAL(storage->getData().getSpectra()[0].getComment(), "comm1")
    TEST_EQUAL(storage->getData().getSpectra()[1].getComment(), "comm3")
    TEST_EQUAL(storage->getData().getSpectra()[2].getComment(), "comm5")

    MSSpectrum snew = storage->getData().getSpectra()[0];

    TEST_EQUAL(snew.size(), 3)
    TEST_REAL_SIMILAR(snew[0].getMZ(), 5)
    TEST_REAL_SIMILAR(snew[0].getIntensity(), 17)
    TEST_REAL_SIMILAR(snew[1].getMZ(), 10)
    TEST_REAL_SIMILAR(snew[1].getIntensity(), 120)
    TEST_REAL_SIMILAR(snew[2].getMZ(), 15)
    TEST_REAL_SIMILAR(snew[2].getIntensity(), 230)

    delete storage;
  }
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  MSDataStoringConsumer * storage = new MSDataStoringConsumer();
  MSDataAggregatingConsumer * agg_consumer = new MSDataAggregatingConsumer(storage);

  MSChromatogram c;
  c.setNativeID("testid");
  agg_consumer->consumeChromatogram(c);

  delete agg_consumer;

  TEST_EQUAL(storage->getData().getNrSpectra(), 0)
  TEST_EQUAL(storage->getData().getNrChromatograms(), 1)
  TEST_EQUAL(storage->getData().getChromatograms()[0].getNativeID(), "testid")

  delete storage;
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  /*
     MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

     transforming_consumer->setExpectedSize(2,0);
     ExperimentalSettings s;
     transforming_consumer->setExperimentalSettings( s );

     TEST_NOT_EQUAL(transforming_consumer, transforming_consumer_nullPointer)
     delete transforming_consumer;
     */
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
