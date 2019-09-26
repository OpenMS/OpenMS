// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/OMSFile.h>

#include <OpenMS/FORMAT/SqliteConnector.h>
#include <sqlite3.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(OMSFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

String oms_path;
IdentificationData ids;

START_SECTION(void store(const String& filename, const IdentificationData& id_data))
{
  // empty input:
  oms_path = OPENMS_GET_TEST_DATA_PATH("OMSFile_test_0.oms");
  OMSFile::store(oms_path, ids);

  sqlite3_stmt* result;
  SqliteConnector con(oms_path);
  con.prepareStatement(&result, "SELECT * FROM version;");
  sqlite3_step(result);
  TEST_EQUAL(sqlite3_column_count(result), 4);
  TEST_STRING_EQUAL(sqlite3_column_name(result, 0), "OMSFile");
  TEST_EQUAL(sqlite3_column_int(result, 0), 1);
  TEST_STRING_EQUAL(sqlite3_column_name(result, 1), "date");
  TEST_STRING_EQUAL(sqlite3_column_name(result, 2), "OpenMS");
  TEST_STRING_EQUAL(sqlite3_column_name(result, 3), "build_date");
  sqlite3_step(result);
  TEST_EQUAL(sqlite3_column_type(result, 0), SQLITE_NULL);
  sqlite3_finalize(result); // free memory

  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), proteins_in, peptides_in);
  // IdentificationData doesn't allow score types with the same name, but different orientations:
  peptides_in[0].setHigherScoreBetter(true);

  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  oms_path = OPENMS_GET_TEST_DATA_PATH("OMSFile_test_1.oms");
  OMSFile::store(oms_path, ids);

}
END_SECTION

START_SECTION(void load(const String& filename, IdentificationData& id_data))
{
  IdentificationData out;
  OMSFile::load(oms_path, out);

  TEST_EQUAL(ids.getParentMolecules().size(), out.getParentMolecules().size());
  TEST_EQUAL(ids.getScoreTypes().size(), out.getScoreTypes().size());
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
