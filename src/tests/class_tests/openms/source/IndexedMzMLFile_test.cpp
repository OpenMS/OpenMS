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
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
///////////////////////////

#include <OpenMS/FORMAT/FileTypes.h>

// for comparison
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

///////////////////////////

START_TEST(IndexedMzMLHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IndexedMzMLHandler* ptr = nullptr;
IndexedMzMLHandler* nullPointer = nullptr;
START_SECTION((IndexedMzMLHandler(String filename) ))
  ptr = new IndexedMzMLHandler(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IndexedMzMLHandler()))
  delete ptr;
END_SECTION

START_SECTION((IndexedMzMLHandler() ))
  ptr = new IndexedMzMLHandler();
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
END_SECTION

START_SECTION((IndexedMzMLHandler(const IndexedMzMLHandler &source)))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  IndexedMzMLHandler file2(file);

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
    IndexedMzMLHandler file;
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    file.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML") );
    TEST_EQUAL(file.getParsingSuccess(), true)
  }

  {
    IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), false)
  }

  {
    IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
    TEST_EQUAL(file.getParsingSuccess(), true)
  }
}
END_SECTION

START_SECTION(( void openFile(String filename) ))
{
  IndexedMzMLHandler file;
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
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrSpectra(), 2)
}
END_SECTION

START_SECTION(( size_t getNrChromatograms() const ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrChromatograms(), 1)
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id)  ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

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
    IndexedMzMLHandler file;
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::ParseError,file.getSpectrumById( 0 ));
  }
}
END_SECTION

START_SECTION(( OpenMS::MSSpectrum getMSSpectrumById(int id)  ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrSpectra(), exp.getSpectra().size())

  OpenMS::MSSpectrum spec = file.getMSSpectrumById(0);
  TEST_EQUAL(spec.size(), exp.getSpectra()[0].size() )
  // TEST_EQUAL(spec.getNativeID(), exp.getSpectra()[0].getNativeID() )

  // Test Exceptions
  TEST_EXCEPTION(Exception::IllegalArgument,file.getMSSpectrumById(-1));
  TEST_EXCEPTION(Exception::IllegalArgument,file.getMSSpectrumById( file.getNrSpectra()+1));

  {
    IndexedMzMLHandler file;
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::ParseError,file.getMSSpectrumById( 0 ));
  }
}
END_SECTION

START_SECTION(( void getMSSpectrumByNativeId(std::string id, OpenMS::MSSpectrum& s) ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrSpectra(), exp.getSpectra().size())

  OpenMS::MSSpectrum spec;
  file.getMSSpectrumByNativeId("controllerType=0 controllerNumber=1 scan=1", spec);
  TEST_EQUAL(spec.size(), exp.getSpectra()[0].size() )
  TEST_EQUAL(spec.getNativeID(), exp.getSpectra()[0].getNativeID() )

  // Test Exceptions
  TEST_EXCEPTION(Exception::IllegalArgument,file.getMSSpectrumByNativeId("TEST", spec));

  {
    IndexedMzMLHandler file;
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::IllegalArgument, file.getMSSpectrumByNativeId( "TEST", spec ));
  }
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int id) ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrChromatograms(), exp.getChromatograms().size())

  OpenMS::Interfaces::ChromatogramPtr chrom = file.getChromatogramById(0);
  TEST_EQUAL(chrom->getTimeArray()->data.size(), exp.getChromatograms()[0].size() )
  TEST_EQUAL(chrom->getIntensityArray()->data.size(), exp.getChromatograms()[0].size() )
}
END_SECTION

START_SECTION(( OpenMS::MSChromatogram getMSChromatogramById(int id) ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrChromatograms(), exp.getChromatograms().size())

  OpenMS::MSChromatogram chrom = file.getMSChromatogramById(0);
  TEST_EQUAL(chrom.size(), exp.getChromatograms()[0].size() )
  TEST_EQUAL(chrom.getNativeID(), exp.getChromatograms()[0].getNativeID() )
}
END_SECTION

START_SECTION(( void getMSChromatogramByNativeId(std::string id, OpenMS::MSChromatogram& c) ))
{
  IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrChromatograms(), exp.getChromatograms().size())

  OpenMS::MSChromatogram chrom;
  file.getMSChromatogramByNativeId("TIC", chrom);
  TEST_EQUAL(chrom.size(), exp.getChromatograms()[0].size() )
  TEST_EQUAL(chrom.getNativeID(), exp.getChromatograms()[0].getNativeID() )

  TEST_EXCEPTION(Exception::IllegalArgument,file.getMSChromatogramByNativeId("TEST", chrom));

  {
    IndexedMzMLHandler file;
    TEST_EXCEPTION(Exception::FileNotFound, file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist")));
    TEST_EQUAL(file.getParsingSuccess(), false)
    TEST_EXCEPTION(Exception::IllegalArgument, file.getMSChromatogramByNativeId( "TEST", chrom ));
  }
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
    TEST_EXCEPTION(Exception::ConversionError, new IndexedMzMLHandler(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_2_broken.mzML")))
  }
  else
  {
    IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_2_broken.mzML"));
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
    IndexedMzMLHandler file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_3_broken.mzML"));
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
    TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError,
      new IndexedMzMLHandler(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_3_broken.mzML")), 
      "Could not convert string '9223372036854775807' to an integer on your system." )
  }
}
END_SECTION
    
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

