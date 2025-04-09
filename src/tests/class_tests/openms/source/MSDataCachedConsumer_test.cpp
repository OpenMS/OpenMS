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
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/FORMAT/CachedMzML.h>

START_TEST(MSDataCachedConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataCachedConsumer* cached_consumer_ptr = nullptr;
MSDataCachedConsumer* cached_consumer_nullPointer = nullptr;

START_SECTION((MSDataCachedConsumer()))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  cached_consumer_ptr = new MSDataCachedConsumer(tmp_filename);
  TEST_NOT_EQUAL(cached_consumer_ptr, cached_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataCachedConsumer()))
    delete cached_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, false);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)

  cached_consumer->setExpectedSize(2,0);
  cached_consumer->consumeSpectrum(exp.getSpectrum(0));
  cached_consumer->consumeSpectrum(exp.getSpectrum(1));
  delete cached_consumer;

  // Check whether it was written to disk correctly...
  {
    // Create the index from the given file
    Internal::CachedMzMLHandler cache;
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> spectra_index = cache.getSpectraIndex();
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the spectrum
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    ifs_.seekg(spectra_index[0]);
    int ms_level = -1;
    double rt = -1.0;
    Internal::CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

    TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(0).size())
    TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(0).size())

    // retrieve the spectrum
    ifs_.seekg(spectra_index[1]);
    Internal::CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

    TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(1).size())
    TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(1).size())
  }
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, false);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  cached_consumer->setExpectedSize(0,1);
  cached_consumer->consumeChromatogram(exp.getChromatogram(0));
  delete cached_consumer;

  // Check whether it was written to disk correctly...
  {
    // Create the index from the given file
    Internal::CachedMzMLHandler cache;
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> chrom_index = cache.getChromatogramIndex();;
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the chromatogram
    OpenSwath::BinaryDataArrayPtr time_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    ifs_.seekg(chrom_index[0]);
    Internal::CachedMzMLHandler::readChromatogramFast(time_array, intensity_array, ifs_);

    TEST_EQUAL(time_array->data.size(), exp.getChromatogram(0).size())
    TEST_EQUAL(intensity_array->data.size(), exp.getChromatogram(0).size())
  }
}
END_SECTION

START_SECTION((MSDataCachedConsumer(String filename, bool clearData=true)))
{
  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, true);

    PeakMap exp;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
    TEST_EQUAL(exp.getNrSpectra() > 0, true)
    MSSpectrum first_spectrum = exp.getSpectrum(0);

    cached_consumer->setExpectedSize(2,0);

    TEST_EQUAL(!exp.getSpectrum(0).empty(), true)

    cached_consumer->consumeSpectrum(exp.getSpectrum(0));

    TEST_EQUAL(exp.getSpectrum(0).size(), 0)
    TEST_EQUAL(exp.getSpectrum(0) == first_spectrum, false)

    delete cached_consumer;
  }
  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, false);

    PeakMap exp;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
    TEST_EQUAL(exp.getNrSpectra() > 0, true)
    MSSpectrum first_spectrum = exp.getSpectrum(0);

    cached_consumer->setExpectedSize(2,0);

    TEST_EQUAL(!exp.getSpectrum(0).empty(), true)

    cached_consumer->consumeSpectrum(exp.getSpectrum(0));

    TEST_EQUAL(!exp.getSpectrum(0).empty(), true)
    TEST_EQUAL(exp.getSpectrum(0) == first_spectrum, true)

    delete cached_consumer;
  }
}
END_SECTION

START_SECTION((void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION([EXTRA] test empty file)
{
  // try an empty file
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, false);
  delete cached_consumer;

  // Check whether it was written to disk correctly...
  {
    // Create the index from the given file
    Internal::CachedMzMLHandler cache;
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> spectra_index = cache.getSpectraIndex();
    TEST_EQUAL(cache.getSpectraIndex().size(), 0)
    TEST_EQUAL(cache.getChromatogramIndex().size(), 0)
  }
}
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  MSDataCachedConsumer * cached_consumer = new MSDataCachedConsumer(tmp_filename, true);

  cached_consumer->setExpectedSize(2,0);
  ExperimentalSettings s;
  cached_consumer->setExperimentalSettings( s );

  TEST_NOT_EQUAL(cached_consumer, cached_consumer_nullPointer)
  delete cached_consumer;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
