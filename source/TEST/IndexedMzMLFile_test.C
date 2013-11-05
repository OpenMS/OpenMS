// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

IndexedMzMLFile* ptr = 0;
IndexedMzMLFile* nullPointer = 0;
START_SECTION((IndexedMzMLFile(String filename) ))
	ptr = new IndexedMzMLFile(OPENMS_GET_TEST_DATA_PATH("small.pwiz.1.1.test.mzML"));
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

START_SECTION(( bool getParsingSuccess() ))
{
  {
    IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist"));
    TEST_EQUAL(file.getParsingSuccess(), false)
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
  file.openFile(OPENMS_GET_TEST_DATA_PATH("fileDoesNotExist"));
  TEST_EQUAL(file.getParsingSuccess(), false)
  file.openFile(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(file.getParsingSuccess(), false)
  file.openFile(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getParsingSuccess(), true)
}
END_SECTION

START_SECTION(( bool getNrSpectra() ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrSpectra(), 2)
}
END_SECTION

START_SECTION(( bool getNrChromatograms() ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(file.getNrChromatograms(), 1)
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id)  ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

	MSExperiment<> exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrSpectra(), exp.getSpectra().size())

  OpenMS::Interfaces::SpectrumPtr spec = file.getSpectrumById(0);
  TEST_EQUAL(spec->getMZArray()->data.size(), exp.getSpectra()[0].size() )
  TEST_EQUAL(spec->getIntensityArray()->data.size(), exp.getSpectra()[0].size() )
}
END_SECTION

START_SECTION(( OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int id) ))
{
  IndexedMzMLFile file(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));

	MSExperiment<> exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"),exp);

  TEST_EQUAL(file.getNrChromatograms(), exp.getChromatograms().size())

  OpenMS::Interfaces::ChromatogramPtr chrom = file.getChromatogramById(0);
  TEST_EQUAL(chrom->getTimeArray()->data.size(), exp.getChromatograms()[0].size() )
  TEST_EQUAL(chrom->getIntensityArray()->data.size(), exp.getChromatograms()[0].size() )
}
END_SECTION
    
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

