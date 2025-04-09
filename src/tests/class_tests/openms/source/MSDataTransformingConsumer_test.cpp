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
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>


START_TEST(MSDataTransformingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataTransformingConsumer* transforming_consumer_ptr = nullptr;
MSDataTransformingConsumer* transforming_consumer_nullPointer = nullptr;

PeakMap expc;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), expc);


START_SECTION((MSDataTransformingConsumer()))
  transforming_consumer_ptr = new MSDataTransformingConsumer();
  TEST_NOT_EQUAL(transforming_consumer_ptr, transforming_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataTransformingConsumer()))
    delete transforming_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  MSDataTransformingConsumer transforming_consumer;

  PeakMap exp = expc;
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  transforming_consumer.setExpectedSize(2,0);
  transforming_consumer.consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), true) // nothing happened
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  MSDataTransformingConsumer transforming_consumer;

  PeakMap exp = expc;
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  transforming_consumer.setExpectedSize(0,1);
  transforming_consumer.consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), true) // nothing happened
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  MSDataTransformingConsumer transforming_consumer;

  transforming_consumer.setExpectedSize(2,0);
  ExperimentalSettings s;
  transforming_consumer.setExperimentalSettings( s );

  NOT_TESTABLE
}
END_SECTION

START_SECTION(( void MSDataTransformingConsumer::setSpectraProcessingFunc( std::function<void (SpectrumType&)> f_spec ) ))
{
  MSDataTransformingConsumer transforming_consumer;

  PeakMap exp = expc;
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  exp.getSpectrum(0).sortByPosition();
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  auto f = [](OpenMS::MSSpectrum & s)
  {
    s.sortByIntensity();
  };

  transforming_consumer.setExpectedSize(2,0);
  transforming_consumer.setSpectraProcessingFunc(f);
  transforming_consumer.consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

}
END_SECTION

START_SECTION(( void MSDataTransformingConsumer::setChromatogramProcessingFunc( std::function<void (ChromatogramType&)> f_chrom ) ))
{
  MSDataTransformingConsumer transforming_consumer;

  PeakMap exp = expc;
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  exp.getChromatogram(0).sortByPosition();
  MSChromatogram first_chromatogram = exp.getChromatogram(0);
  
  auto f2 = [](OpenMS::MSChromatogram & c)
  {
    c.sortByIntensity();
  };

  transforming_consumer.setExpectedSize(0,1);
  transforming_consumer.setChromatogramProcessingFunc(f2);
  transforming_consumer.consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), false) // something happened
  TEST_EQUAL(first_chromatogram.isSorted(), true)
  TEST_EQUAL(exp.getChromatogram(0).isSorted(), false)

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
