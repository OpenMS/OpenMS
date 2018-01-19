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

#include <OpenMS/FORMAT/IndexedMzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>

// for comparison
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IndexedMzMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexedMzMLFile* ptr = nullptr;
IndexedMzMLFile* nullPointer = nullptr;
START_SECTION((IndexedMzMLFile(String filename) ))
	ptr = new IndexedMzMLFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IndexedMzMLFile()))
	delete ptr;
END_SECTION

START_SECTION((IndexedMzMLFile() ))
	ptr = new IndexedMzMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION

START_SECTION((IndexedMzMLFile(const IndexedMzMLFile &source)))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  IndexedMzMLFile file2(file);

  TEST_EQUAL(file.getParsingSuccess(), file2.getParsingSuccess())
  TEST_EQUAL(file.getNrSpectra(), file2.getNrSpectra())
  TEST_EQUAL(file.getNrChromatograms(), file2.getNrChromatograms())

  ABORT_IF(file.getNrSpectra() != 2)
  TEST_EQUAL(file.getSpectrumById(0)->getMZArray()->data == file2.getSpectrumById(0)->getMZArray()->data, true)
  TEST_EQUAL(file.getSpectrumById(0)->getIntensityArray()->data == file2.getSpectrumById(0)->getIntensityArray()->data, true)
  TEST_EQUAL(file.getSpectrumById(1)->getMZArray()->data == file2.getSpectrumById(1)->getMZArray()->data, true)
  TEST_EQUAL(file.getSpectrumById(1)->getIntensityArray()->data == file2.getSpectrumById(1)->getIntensityArray()->data, true)
  ABORT_IF(file.getNrChromatograms() != 1)
  TEST_EQUAL(file.getChromatogramById(0)->getTimeArray()->data == file2.getChromatogramById(0)->getTimeArray()->data, true)
  TEST_EQUAL(file.getChromatogramById(0)->getIntensityArray()->data == file2.getChromatogramById(0)->getIntensityArray()->data, true)
  /*
  TEST_EQUAL(file.getChromatogramById(0) == file2.getChromatogramById(0), true)
  TEST_EQUAL(file.getSpectrumById(1), file2.getSpectrumById(1))
  */
}
END_SECTION

START_SECTION(( bool getParsingSuccess() const))
{
  {
    IndexedMzMLFile file;
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    file.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML") );
    TEST_EQUAL(file.getParsingSuccess(), true)
  }

  {
    IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), false)
  }

  {
    IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), true)
  }
}
END_SECTION

START_SECTION(( void openFile(String filename) ))
{
  IndexedMzMLFile file;
  TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")))
  TEST_EQUAL(file.getParsingSuccess(), false)
  file.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(file.getParsingSuccess(), false)
  file.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getParsingSuccess(), true)
}
END_SECTION

START_SECTION(( size_t getNrSpectra() const ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrSpectra(), 2)
}
END_SECTION

START_SECTION(( size_t getNrChromatograms() const ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrChromatograms(), 1)
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id)  ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrSpectra(), exp.getSpectra().size())

  OpenMS::Interfaces::SpectrumPtr spec = file.getSpectrumById(0);
  TEST_EQUAL(spec->getMZArray()->data.size(), exp.getSpectra()[0].size() )
  TEST_EQUAL(spec->getIntensityArray()->data.size(), exp.getSpectra()[0].size() )

  // Test Exceptions
  TEST_EXCEPTION(Exception::IllegalArgument,file.getSpectrumById(-1));
  TEST_EXCEPTION(Exception::IllegalArgument,file.getSpectrumById( file.getNrSpectra()+1));

  {
    IndexedMzMLFile file;
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::ParseError,file.getSpectrumById( 0 ));
  }
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int id) ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrChromatograms(), exp.getChromatograms().size())

  OpenMS::Interfaces::ChromatogramPtr chrom = file.getChromatogramById(0);
  TEST_EQUAL(chrom->getTimeArray()->data.size(), exp.getChromatograms()[0].size() )
  TEST_EQUAL(chrom->getIntensityArray()->data.size(), exp.getChromatograms()[0].size() )
}
END_SECTION

START_SECTION(([EXTRA] load broken file))
{

  // Contains an unparseable value (2^64) in the indexListOffset field that
  // will not fit into a long long.
  // NOTE: this will not be true on all systems, if long long is larger than 64
  // bits it will fit, however parsing will fail since the file is not actually
  // 2^64 bit long...
  if ( sizeof(long long)*8 <= 64 )
  {
    TEST_EXCEPTION(Exception::ConversionError, new IndexedMzMLFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_2_broken.mzML")))
  }
  else
  {
    IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_2_broken.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), false)
  }
}
END_SECTION

START_SECTION(([EXTRA] load broken file))
{

  // Contains an value (2^63-1) in the indexListOffset field that should not
  // trigger an exception - however parsing will fail since the file is
  // actually shorter.
  if (sizeof(std::streampos)*8 > 32 )
  {
    IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_3_broken.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), false)
  }
  else
  {
    // 
    // On systems that use 32 bit or less to represent std::streampos, we
    // cannot fit our value in the indexListOffset (2^63-1) into std::streampos
    // -> this should throw an exception in the constructor which can test here
    // instead when loading the file.
    //
    // This code path is hard to test on most machines since almost all modern
    // compilers and filesystems support file access for files > 2 GB 
    // Manually, one can cast the indexoffset variable to int to trigger this
    // behavior in IndexedMzMLDecoder.cpp
    // 
    TEST_EXCEPTION_WITH_MESSAGE (Exception::ConversionError, 
      new IndexedMzMLFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_3_broken.mzML")), 
      "Could not convert string '9223372036854775807' to an integer on your system." )
  }
}
END_SECTION
    
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

