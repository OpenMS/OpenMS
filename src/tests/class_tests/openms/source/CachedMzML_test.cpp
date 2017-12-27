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
#include <OpenMS/FORMAT/CachedMzML.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;


CachedmzML cacheFile(std::string & tmp_filename, PeakMap& exp)
{
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedmzML cache;
  cache.writeMemdump(exp, tmp_filename);
  // Create the index from the given file
  cache.createMemdumpIndex(tmp_filename);
  return cache;
}

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
  CachedmzML cache;
  cache.writeMemdump(exp, tmp_filename);

  // Check whether spectra were written to disk correctly...
  {
    // Create the index from the given file
    cache.createMemdumpIndex(tmp_filename);
    std::vector<std::streampos> spectra_index = cache.getSpectraIndex();
    TEST_EQUAL(spectra_index.size(), 4)
    std::ifstream ifs_(tmp_filename.c_str(), std::ios::binary);

    // retrieve the spectrum
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    int ms_level = -1;
    double rt = -1.0;
    for (int i = 0; i  < 4; i++)
    {
      ifs_.seekg(spectra_index[i]);
      CachedmzML::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

      TEST_EQUAL(mz_array->data.size(), exp.getSpectrum(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getSpectrum(i).size())
    }
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
    for (int i = 0; i  < 2; i++)
    {
      ifs_.seekg(chrom_index[i]);
      CachedmzML::readChromatogramFast(time_array, intensity_array, ifs_);

      TEST_EQUAL(time_array->data.size(), exp.getChromatogram(i).size())
      TEST_EQUAL(intensity_array->data.size(), exp.getChromatogram(i).size())
    }
  }
}
END_SECTION

START_SECTION(( void writeMemdump(MapType& exp, String out )))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedmzML cache;
  cache.writeMemdump(exp, tmp_filename);

  NOT_TESTABLE // not testable independently, see testCaching
}
END_SECTION

START_SECTION(( void writeMetadata(MapType exp, String out_meta, bool addCacheMetaValue=false) ))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)

  // Cache the experiment to a temporary file
  CachedmzML cache;

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
CachedmzML cache_ = cacheFile(tmp_filename, exp);

START_SECTION(( void readMemdump(MapType& exp_reading, String filename) const ))
{

  std::string tmp_filename;
  PeakMap exp;
  CachedmzML cache = cacheFile(tmp_filename, exp);

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

START_SECTION(( void createMemdumpIndex(String filename) ))
{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  // Load experiment
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  // Cache the experiment to a temporary file
  CachedmzML cache;
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
    CachedmzML::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

    TEST_EQUAL(mz_array->data.size() > 0, true)
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
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedmzML::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt),
      "filestream in: Read an invalid spectrum length, something is wrong here. Aborting.")

    // should not read after the file ends
    ifs_.seekg(spectra_index.back() * 20);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedmzML::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt),
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
    CachedmzML::readChromatogramFast(time_array, intensity_array, ifs_);

    TEST_EQUAL(time_array->data.size() > 0, true)
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
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedmzML::readChromatogramFast(time_array, intensity_array, ifs_),
      "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.")

    // should not read after the file ends
    ifs_.seekg(chrom_index.back() * 20);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, CachedmzML::readChromatogramFast(time_array, intensity_array, ifs_),
      "filestream in: Read an invalid chromatogram length, something is wrong here. Aborting.")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

