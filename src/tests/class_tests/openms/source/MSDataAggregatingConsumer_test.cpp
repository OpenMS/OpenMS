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
  agg_consumer_ptr = new MSDataAggregatingConsumer(agg_consumer_nullPointer); // don't do that ...
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
