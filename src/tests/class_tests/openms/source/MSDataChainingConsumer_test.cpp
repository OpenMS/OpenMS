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

#include <OpenMS/FORMAT/DATAACCESS/MSDataChainingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/NoopMSDataConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/MzMLFile.h>


START_TEST(MSDataChainingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataChainingConsumer* chaining_consumer_ptr = nullptr;
MSDataChainingConsumer* chaining_consumer_nullPointer = nullptr;

START_SECTION((MSDataChainingConsumer()))
  chaining_consumer_ptr = new MSDataChainingConsumer();
  TEST_NOT_EQUAL(chaining_consumer_ptr, chaining_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataChainingConsumer()))
    delete chaining_consumer_ptr;
END_SECTION

START_SECTION(( MSDataChainingConsumer(std::vector<IMSDataConsumer*> consumers) ))
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  chaining_consumer_ptr = new MSDataChainingConsumer(consumer_list);
  TEST_NOT_EQUAL(chaining_consumer_ptr, chaining_consumer_nullPointer)
  delete chaining_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), true) // nothing happened

  for (auto& consumer : consumer_list)
  {
    delete consumer;
  }
  delete chaining_consumer;
}
END_SECTION

START_SECTION(([EXTRA] void consumeSpectrum(SpectrumType & s)))
{
   auto f = [](OpenMS::MSSpectrum & s)
  {
    s.sortByIntensity();
  };

  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setSpectraProcessingFunc(f);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(transforming_consumer);
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)


  // note how the transforming consumer still works as deleting the chaining
  // consumer does not take ownership of the consumers
  transforming_consumer->consumeSpectrum(exp.getSpectrum(0) );

  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  for (auto& consumer : consumer_list)
  {
    delete consumer;
  }
  delete chaining_consumer;
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  chaining_consumer->setExpectedSize(0,1);
  chaining_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), true) // nothing happened

  for (auto& consumer : consumer_list)
  {
    delete consumer;
  }
  delete chaining_consumer;
}
END_SECTION

START_SECTION(([EXTRA]void consumeChromatogram(ChromatogramType & c)))
{
  auto f2 = [](OpenMS::MSChromatogram & c)
  {
    c.sortByIntensity();
  };
  
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setChromatogramProcessingFunc(f2);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(transforming_consumer);
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  chaining_consumer->setExpectedSize(0,1);
  chaining_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), false) // something happened
  TEST_EQUAL(first_chromatogram.isSorted(), true)
  TEST_EQUAL(exp.getChromatogram(0).isSorted(), false)

  for (auto& consumer : consumer_list)
  {
    delete consumer;
  }
  delete chaining_consumer;
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer();

  chaining_consumer->setExpectedSize(2,0);
  ExperimentalSettings s;
  chaining_consumer->setExperimentalSettings( s );

  TEST_NOT_EQUAL(chaining_consumer, chaining_consumer_nullPointer)
  delete chaining_consumer;
}
END_SECTION

START_SECTION(( void appendConsumer(IMSDataConsumer * consumer) ))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  auto f = [](OpenMS::MSSpectrum & s)
  {
    s.sortByIntensity();
  };
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setSpectraProcessingFunc(f);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);
  chaining_consumer->appendConsumer(transforming_consumer);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  for (auto& consumer : consumer_list)
  {
    delete consumer;
  }
  delete transforming_consumer;
  delete chaining_consumer;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
