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
#include <OpenMS/FORMAT/CachedMzML.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

START_TEST(CachedmzML, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CachedmzML* ptr = nullptr;
CachedmzML* nullPointer = nullptr;

START_SECTION(CachedmzML())
{
  ptr = new CachedmzML();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CachedmzML())
{
  delete ptr;
}
END_SECTION

// Load experiment
PeakMap exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);

std::string tmpf;
NEW_TMP_FILE(tmpf);

// Cache the experiment to a temporary file
CachedmzML::store(tmpf, exp);
CachedmzML cache_example;
CachedmzML::load(tmpf, cache_example);

// see also MSDataCachedConsumer_test.cpp -> consumeSpectrum
// this is a complete test of the caching object
START_SECTION(( [EXTRA] testCaching))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedmzML::store(tmp_filename, exp);

  // Check whether spectra were written to disk correctly...
  {
    // Create the index from the given file
    CachedmzML cache;
    CachedmzML::load(tmp_filename, cache);

    TEST_EQUAL(cache.getNrSpectra(), 4)

    // retrieve the spectrum (old interface)
    for (int i = 0; i < 4; i++)
    {
      TEST_EQUAL(cache.getSpectrum(i).size(), exp.getSpectrum(i).size())

      // identical except DataProcessing (and extra data arrays -- does not have all fields)
      auto tmp1 = cache.getSpectrum(i);
      auto tmp2 = exp.getSpectrum(i);
      tmp1.getDataProcessing().clear();
      tmp2.getDataProcessing().clear();
      tmp1.getFloatDataArrays().clear(); // clear for now, see test below
      tmp2.getFloatDataArrays().clear(); // clear for now, see test below
      TEST_TRUE(tmp1 == tmp2)
    }

    // test spec 1
    auto scomp = exp.getSpectrum(1);
    TEST_EQUAL(scomp.getFloatDataArrays().size(), 2)
    TEST_EQUAL(scomp.getIntegerDataArrays().size(), 0)
    TEST_EQUAL(scomp.getStringDataArrays().size(), 0)

    // test spec 1
    auto s = cache.getSpectrum(1);
    TEST_EQUAL(s.getFloatDataArrays().size(), 2)
    TEST_EQUAL(s.getIntegerDataArrays().size(), 0)
    TEST_EQUAL(s.getStringDataArrays().size(), 0)

    TEST_EQUAL(s.getFloatDataArrays()[0].getName(), scomp.getFloatDataArrays()[0].getName())
    TEST_EQUAL(s.getFloatDataArrays()[1].getName(), scomp.getFloatDataArrays()[1].getName())
    TEST_EQUAL(s.getFloatDataArrays()[0].getName(), "signal to noise array")
    TEST_EQUAL(s.getFloatDataArrays()[1].getName(), "user-defined name")

    TEST_EQUAL(s.getFloatDataArrays()[0].size(), scomp.getFloatDataArrays()[0].size())
    TEST_EQUAL(s.getFloatDataArrays()[1].size(), scomp.getFloatDataArrays()[1].size())

    for (Size k = 0; k < s.getFloatDataArrays()[0].size(); k++)
    {
      TEST_REAL_SIMILAR(s.getFloatDataArrays()[0][k], scomp.getFloatDataArrays()[0][k])
    }

    for (Size k = 0; k < s.getFloatDataArrays()[1].size(); k++)
    {
      TEST_REAL_SIMILAR(s.getFloatDataArrays()[1][k], scomp.getFloatDataArrays()[1][k])
    }

  }

  // Check whether chromatograms were written to disk correctly...
  {
    // Create the index from the given file
    CachedmzML cache;
    CachedmzML::load(tmp_filename, cache);

    TEST_EQUAL(cache.getNrChromatograms(), 2)

    // retrieve the chromatogram
    for (int i = 0; i < 2; i++)
    {
      TEST_EQUAL(cache.getChromatogram(i).size(), exp.getChromatogram(i).size())
      TEST_EQUAL(cache.getChromatogram(i).getNativeID(), exp.getChromatogram(i).getNativeID())
      TEST_EQUAL(cache.getChromatogram(i).getInstrumentSettings() == exp.getChromatogram(i).getInstrumentSettings(), true)

      // identical except DataProcessing
      auto tmp1 = cache.getChromatogram(i);
      auto tmp2 = exp.getChromatogram(i);
      tmp1.getDataProcessing().clear();
      tmp2.getDataProcessing().clear();
      TEST_TRUE(tmp1 == tmp2)
    }

  }
}
END_SECTION

START_SECTION(( size_t getNrSpectra() const ))
    TEST_EQUAL(cache_example.getNrSpectra(), 4)
END_SECTION

START_SECTION(( size_t getNrChromatograms() const ))
    TEST_EQUAL(cache_example.getNrChromatograms(), 2)
END_SECTION

START_SECTION(( const MSExperiment& getMetaData() const ))
    TEST_EQUAL(cache_example.getMetaData().size(), 4)
    TEST_EQUAL(cache_example.getMetaData().getNrSpectra(), 4)
    TEST_EQUAL(cache_example.getMetaData().getNrChromatograms(), 2)
END_SECTION

START_SECTION(( const MSExperiment& getMetaData() const ))
{
  TEST_EQUAL(cache_example.getNrSpectra(), cache_example.getMetaData().getNrSpectra())
  for (int i = 0; i < 4; i++)
  {
    // identical except DataProcessing
    SpectrumSettings tmp1 = cache_example.getMetaData()[i];
    SpectrumSettings tmp2 = exp.getSpectrum(i);
    tmp1.getDataProcessing().clear();
    tmp2.getDataProcessing().clear();
    TEST_TRUE(tmp1 == tmp2)
  }

  TEST_EQUAL(cache_example.getNrChromatograms(), cache_example.getMetaData().getNrChromatograms())
  for (int i = 0; i < 2; i++)
  {
    // identical except DataProcessing
    ChromatogramSettings tmp1 = cache_example.getMetaData().getChromatograms()[i];
    ChromatogramSettings tmp2 = exp.getChromatogram(i);
    tmp1.getDataProcessing().clear();
    tmp2.getDataProcessing().clear();
    TEST_TRUE(tmp1 == tmp2)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


#pragma clang diagnostic pop

