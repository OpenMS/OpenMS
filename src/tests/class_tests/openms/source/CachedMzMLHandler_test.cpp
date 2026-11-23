// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshadow"
#endif

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

void cacheFile(CachedMzMLHandler& cache, std::string & tmp_filename, PeakMap& exp)
{
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  cache.writeMemdump(exp, tmp_filename);
  // Create the index from the given file
  cache.createMemdumpIndex(tmp_filename);
}

START_TEST(CachedMzMLHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CachedMzMLHandler* ptr = nullptr;
CachedMzMLHandler* nullPointer = nullptr;

START_SECTION(CachedMzMLHandler())
{
  ptr = new CachedMzMLHandler();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CachedMzMLHandler())
{
  delete ptr;
}
END_SECTION

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
  CachedMzMLHandler cache;
  cache.writeMemdump(exp, tmp_filename);

  // Check whether spectra were written to disk correctly...
  {
    // Create the index from the given file
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> spectra_index = cache.getSpectraIndex();
    TEST_EQUAL(spectra_index.size(), 4)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the spectrum (old interface)
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    int ms_level = -1;
    double rt = -1.0;
    for (int i = 0; i < 4; i++)
    {
      ifs_.seekg(spectra_index[i]);
      CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);
      TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(i).size())
    }

    // retrieve the spectrum (new interface)
    ms_level = -1;
    rt = -1.0;
    for (int i = 0; i < 4; i++)
    {
      ifs_.seekg(spectra_index[i]);
      std::vector<OpenSwath::BinaryDataArrayPtr> darray = CachedMzMLHandler::readSpectrumFast(ifs_, ms_level, rt);
      TEST_EQUAL(darray.size() >= 2, true)
      mz_array = darray[0];
      intensity_array = darray[1];

      TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(i).size())
    }

    // test spec 1
    ifs_.seekg(spectra_index[1]);
    std::vector<OpenSwath::BinaryDataArrayPtr> darray = CachedMzMLHandler::readSpectrumFast(ifs_, ms_level, rt);
    TEST_EQUAL(darray.size(), 4)
    TEST_EQUAL(darray[0]->description, "") // mz
    TEST_EQUAL(darray[1]->description, "") // intensity
    TEST_EQUAL(darray[2]->description, "signal to noise array")
    TEST_EQUAL(darray[3]->description, "user-defined name")

  }

  // Check whether chromatograms were written to disk correctly...
  {
    // Create the index from the given file
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> chrom_index = cache.getChromatogramIndex();;
    TEST_EQUAL(chrom_index.size(), 2)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the chromatogram
    OpenSwath::BinaryDataArrayPtr time_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    for (int i = 0; i < 2; i++)
    {
      ifs_.seekg(chrom_index[i]);
      CachedMzMLHandler::readChromatogramFast(time_array, intensity_array, ifs_);

      TEST_EQUAL(time_array->data.size(), exp.getChromatogram(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getChromatogram(i).size())
    }

    // retrieve the chromatogram (new interface)
    for (int i = 0; i < 2; i++)
    {
      ifs_.seekg(chrom_index[i]);
      std::vector<OpenSwath::BinaryDataArrayPtr> darray = CachedMzMLHandler::readChromatogramFast(ifs_);
      TEST_EQUAL(darray.size() >= 2, true)
      time_array = darray[0];
      intensity_array = darray[1];

      TEST_EQUAL(time_array->data.size(), exp.getChromatogram(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getChromatogram(i).size())
    }
  }
}
END_SECTION

START_SECTION(( void writeMemdump(MapType& exp, std::string out )))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedMzMLHandler cache;
  cache.writeMemdump(exp, tmp_filename);

  NOT_TESTABLE // not testable independently, see testCaching
}
END_SECTION

START_SECTION(( void writeMetadata(MapType exp, std::string out_meta, bool addCacheMetaValue=false) ))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedMzMLHandler cache;

  PeakMap meta_exp;
  // without adding the cache value, the meta data of the two experiments should be equal
  cache.writeMetadata(exp, tmp_filename, false);
  MzMLFile().load(tmp_filename, meta_exp);
  TEST_EQUAL( (ExperimentalSettings)(meta_exp), (ExperimentalSettings)(exp) )
  TEST_EQUAL( (SpectrumSettings)(meta_exp.getSpectrum(0)), (SpectrumSettings)(exp.getSpectrum(0)) )
  TEST_EQUAL( (ChromatogramSettings)(meta_exp.getChromatogram(0)), (ChromatogramSettings)(exp.getChromatogram(0)) )

  // without adding the cache value, the meta data except the "cache" meta value should be equal
  cache.writeMetadata(exp, tmp_filename, true);
  MzMLFile().load(tmp_filename, meta_exp);
  TEST_EQUAL( (ExperimentalSettings)(meta_exp), (ExperimentalSettings)(exp) )

}
END_SECTION

// Create a single CachedMzML file and use it for the following computations
// (may be somewhat faster)
std::string tmp_filename;
PeakMap exp;
CachedMzMLHandler cache_;
cacheFile(cache_, tmp_filename, exp);

START_SECTION(( void readMemdump(MapType& exp_reading, std::string filename) const ))
{

  std::string tmp_filename;
  PeakMap exp;
  CachedMzMLHandler cache;
  cacheFile(cache, tmp_filename, exp);

  PeakMap exp_new;
  cache.readMemdump(exp_new, tmp_filename);

  TEST_EQUAL(exp_new.size(), exp.size())
  TEST_EQUAL(exp_new.getChromatograms().size(), exp.getChromatograms().size())

  std::string unused_tmp_filename;
  NEW_TMP_FILE(unused_tmp_filename);
  TEST_EXCEPTION(Exception::FileNotFound, cache.readMemdump(exp_new, unused_tmp_filename) )
  TEST_EXCEPTION(Exception::ParseError, cache.readMemdump(exp_new, OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML") ) )
}
END_SECTION

START_SECTION(( void createMemdumpIndex(std::string filename) ))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  // Cache the experiment to a temporary file
  CachedMzMLHandler cache;
  cache.writeMemdump(exp, tmp_filename);

  // create the memory dump
  cache.createMemdumpIndex(tmp_filename);

  // check whether we actually did read something
  TEST_EQUAL( cache.getSpectraIndex().size(), 4);
  TEST_EQUAL( cache.getChromatogramIndex().size(), 2);

  // Test error conditions
  std::string unused_tmp_filename;
  NEW_TMP_FILE(unused_tmp_filename);
  TEST_EXCEPTION(Exception::FileNotFound, cache.createMemdumpIndex(unused_tmp_filename) )
  TEST_EXCEPTION(Exception::ParseError, cache.createMemdumpIndex(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML") ) )
}
END_SECTION

START_SECTION(([EXTRA] spectra and chromatograms without peaks keep their data arrays))
{
  // A spectrum without peaks may still have data arrays, e.g. the empty ion mobility array of a Bruker frame
  // without peaks. They have to be written, or what follows them is read from the wrong position.
  PeakMap exp;
  MSSpectrum empty_spectrum;
  empty_spectrum.setRT(1.0);
  empty_spectrum.getFloatDataArrays().resize(1);
  empty_spectrum.getFloatDataArrays()[0].setName("Ion Mobility");
  empty_spectrum.getIntegerDataArrays().resize(1);
  exp.addSpectrum(empty_spectrum);
  MSSpectrum spectrum;
  spectrum.setRT(2.0);
  spectrum.setMSLevel(2);
  spectrum.push_back(Peak1D(100.0, 10.0f));
  spectrum.push_back(Peak1D(200.0, 20.0f));
  spectrum.getFloatDataArrays().resize(1);
  spectrum.getFloatDataArrays()[0].setName("Ion Mobility");
  spectrum.getFloatDataArrays()[0].push_back(0.8f);
  spectrum.getFloatDataArrays()[0].push_back(0.9f);
  exp.addSpectrum(spectrum);
  MSChromatogram empty_chromatogram;
  empty_chromatogram.getFloatDataArrays().resize(1);
  MSChromatogram chromatogram;
  chromatogram.push_back(ChromatogramPeak(5.0, 50.0));
  chromatogram.getFloatDataArrays().resize(1);
  chromatogram.getFloatDataArrays()[0].setName("noise");
  chromatogram.getFloatDataArrays()[0].push_back(2.5f);
  exp.setChromatograms({empty_chromatogram, chromatogram});

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  CachedMzMLHandler cache;
  cache.writeMemdump(exp, tmp_filename);

  PeakMap read;
  cache.readMemdump(read, tmp_filename);
  TEST_EQUAL(read.size(), 2)
  ABORT_IF(read.size() != 2)
  TEST_EQUAL(read[0].size(), 0)
  TEST_REAL_SIMILAR(read[0].getRT(), 1.0)
  TEST_EQUAL(read[0].getFloatDataArrays().size(), 2) // the integer array is read back as a float array
  TEST_EQUAL(read[1].size(), 2)
  TEST_REAL_SIMILAR(read[1].getRT(), 2.0)
  TEST_EQUAL(read[1].getMSLevel(), 2)
  TEST_REAL_SIMILAR(read[1][1].getMZ(), 200.0)
  TEST_REAL_SIMILAR(read[1][1].getIntensity(), 20.0)
  TEST_EQUAL(read[1].getFloatDataArrays().size(), 1)
  ABORT_IF(read[1].getFloatDataArrays().size() != 1)
  TEST_STRING_EQUAL(read[1].getFloatDataArrays()[0].getName(), "Ion Mobility")
  TEST_EQUAL(read[1].getFloatDataArrays()[0].size(), 2)
  TEST_EQUAL(read.getChromatograms().size(), 2)
  ABORT_IF(read.getChromatograms().size() != 2)
  TEST_EQUAL(read.getChromatograms()[0].size(), 0)
  TEST_EQUAL(read.getChromatograms()[1].size(), 1)
  TEST_REAL_SIMILAR(read.getChromatograms()[1][0].getRT(), 5.0)
  // the values of a chromatogram's data arrays are read back as well
  TEST_EQUAL(read.getChromatograms()[1].getFloatDataArrays().size(), 1)
  ABORT_IF(read.getChromatograms()[1].getFloatDataArrays().size() != 1)
  TEST_STRING_EQUAL(read.getChromatograms()[1].getFloatDataArrays()[0].getName(), "noise")
  TEST_EQUAL(read.getChromatograms()[1].getFloatDataArrays()[0].size(), 1)
  ABORT_IF(read.getChromatograms()[1].getFloatDataArrays()[0].size() != 1)
  TEST_REAL_SIMILAR(read.getChromatograms()[1].getFloatDataArrays()[0][0], 2.5)

  // the index finds every spectrum and chromatogram where it begins
  cache.createMemdumpIndex(tmp_filename);
  TEST_EQUAL(cache.getSpectraIndex().size(), 2)
  TEST_EQUAL(cache.getChromatogramIndex().size(), 2)
  ABORT_IF(cache.getSpectraIndex().size() != 2 || cache.getChromatogramIndex().size() != 2)
  std::ifstream ifs(tmp_filename.c_str(), std::ios::binary);
  ifs.seekg(cache.getSpectraIndex()[1]);
  int ms_level = -1;
  double rt = -1.0;
  std::vector<OpenSwath::BinaryDataArrayPtr> data = CachedMzMLHandler::readSpectrumFast(ifs, ms_level, rt);
  TEST_EQUAL(ms_level, 2)
  TEST_REAL_SIMILAR(rt, 2.0)
  TEST_EQUAL(data.size(), 3)
  ABORT_IF(data.size() != 3)
  TEST_EQUAL(data[0]->data.size(), 2)
  ifs.seekg(cache.getChromatogramIndex()[1]);
  data = CachedMzMLHandler::readChromatogramFast(ifs);
  TEST_EQUAL(data[0]->data.size(), 1)
}
END_SECTION

START_SECTION(( const std::vector<std::streampos>& getSpectraIndex() const ))
{
  TEST_EQUAL( cache_.getSpectraIndex().size(), 4);
}
END_SECTION

START_SECTION(( const std::vector<std::streampos>& getChromatogramIndex() const ))
{
  TEST_EQUAL( cache_.getChromatogramIndex().size(), 2);
}
END_SECTION

START_SECTION(static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs, int& ms_level, double& rt))
{

  // Check whether spectra were written to disk correctly...
  {
    // Create the index from the given file
    std::vector<std::streampos> spectra_index = cache_.getSpectraIndex();
    TEST_EQUAL(spectra_index.size(), 4)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the spectrum
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    ifs_.seekg(spectra_index[0]);
    int ms_level = -1;
    double rt = -1.0;
    CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

    TEST_EQUAL(!mz_array->data.empty(), true)
    TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(0).size())
    TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(0).size())

    TEST_EQUAL(ms_level, 1)
    TEST_REAL_SIMILAR(rt, 5.1)

    for (Size i = 0; i < mz_array->data.size(); i++)
    {
      TEST_REAL_SIMILAR(mz_array->data[i], exp.getSpectrum(0)[i].getMZ())
      TEST_REAL_SIMILAR(intensity_array->data[i], exp.getSpectrum(0)[i].getIntensity())
    }
  }

  // Check error conditions
  {
    // Create the index from the given file
    std::vector<std::streampos> spectra_index = cache_.getSpectraIndex();
    TEST_EQUAL(spectra_index.size(), 4)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the spectrum
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    int ms_level = -1;
    double rt = -1.0;

    // should not read before the file starts
    ifs_.seekg( -1 );
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt),
      "filestream in: Read an invalid spectrum length, something is wrong here. Aborting.")

    // should not read after the file ends
    ifs_.seekg(spectra_index.back() * 20);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedMzMLHandler::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt),
      "filestream in: Read an invalid spectrum length, something is wrong here. Aborting.")
  }

}
END_SECTION

START_SECTION( static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1, OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs) )
{
  // Check whether chromatograms were written to disk correctly...
  {
    // Create the index from the given file
    std::vector<std::streampos> chrom_index = cache_.getChromatogramIndex();;
    TEST_EQUAL(chrom_index.size(), 2)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the chromatogram
    OpenSwath::BinaryDataArrayPtr time_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);

    ifs_.seekg(chrom_index[0]);
    CachedMzMLHandler::readChromatogramFast(time_array, intensity_array, ifs_);

    TEST_EQUAL(!time_array->data.empty(), true)
    TEST_EQUAL(time_array->data.size(), exp.getChromatogram(0).size())
    TEST_EQUAL(intensity_array->data.size(), exp.getChromatogram(0).size())

    for (Size i = 0; i < time_array->data.size(); i++)
    {
      TEST_REAL_SIMILAR(time_array->data[i], exp.getChromatogram(0)[i].getRT())
      TEST_REAL_SIMILAR(intensity_array->data[i], exp.getChromatogram(0)[i].getIntensity())
    }
  }

  // Check error conditions
  {
    // Create the index from the given file
    std::vector<std::streampos> chrom_index = cache_.getChromatogramIndex();;
    TEST_EQUAL(chrom_index.size(), 2)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the chromatogram
    OpenSwath::BinaryDataArrayPtr time_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);

    // should not read before the file starts
    ifs_.seekg( -1 );
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedMzMLHandler::readChromatogramFast(time_array, intensity_array, ifs_),
      "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.")

    // should not read after the file ends
    ifs_.seekg(chrom_index.back() * 20);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedMzMLHandler::readChromatogramFast(time_array, intensity_array, ifs_),
      "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST

#ifdef __clang__
  #pragma clang diagnostic pop
#endif
