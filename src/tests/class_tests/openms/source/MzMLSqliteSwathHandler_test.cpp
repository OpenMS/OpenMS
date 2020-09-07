// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>
///////////////////////////

#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QFile>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

///////////////////////////

START_TEST(SqMassFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLSqliteSwathHandler* ptr = nullptr;
MzMLSqliteSwathHandler* nullPointer = nullptr;

std::string tmp_filename;
NEW_TMP_FILE(tmp_filename);
MSExperiment exp_orig;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SwathFile.mzML"), exp_orig);
SqMassFile file;
file.store(tmp_filename, exp_orig);

START_SECTION((MzMLSqliteSwathHandler()))
  ptr = new MzMLSqliteSwathHandler(tmp_filename);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzMLSqliteSwathHandler()))
  delete ptr;
END_SECTION

TOLERANCE_RELATIVE(1.0005)

START_SECTION(std::vector<OpenSwath::SwathMap> readSwathWindows())
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  std::vector<OpenSwath::SwathMap> maps = handler.readSwathWindows();
  TEST_EQUAL(maps.size(), 5)

  TEST_EQUAL(maps[0].ms1, false)
  TEST_REAL_SIMILAR(maps[0].lower, 400.0)
  TEST_REAL_SIMILAR(maps[0].center, 412.5)
  TEST_REAL_SIMILAR(maps[0].upper, 425.0)
  TEST_REAL_SIMILAR(maps[1].lower, 425.0)
  TEST_REAL_SIMILAR(maps[1].upper, 450.0)
  TEST_REAL_SIMILAR(maps[4].lower, 500.0)
  TEST_REAL_SIMILAR(maps[4].upper, 525.0)
}
END_SECTION

START_SECTION(std::vector<int> readMS1Spectra())
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  TEST_EQUAL(handler.readMS1Spectra().size(), 19)
  TEST_EQUAL(handler.readMS1Spectra()[0], 0)
  TEST_EQUAL(handler.readMS1Spectra()[18], 108)
}
END_SECTION

START_SECTION(std::vector<int> readSpectraForWindow(const OpenSwath::SwathMap& swath_map))
{
  MzMLSqliteSwathHandler handler(tmp_filename);

  std::vector<OpenSwath::SwathMap> maps = handler.readSwathWindows();
  TEST_EQUAL(maps.size(), 5)

  TEST_EQUAL(handler.readSpectraForWindow(maps[0]).size(), 19)
  TEST_EQUAL(handler.readSpectraForWindow(maps[0])[0], 1)
  TEST_EQUAL(handler.readSpectraForWindow(maps[0])[18], 109)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1]).size(), 19)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1])[0], 2)
  TEST_EQUAL(handler.readSpectraForWindow(maps[1])[18], 110)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

