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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>

///////////////////////////

START_TEST(IdentificationData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IdentificationData* ptr = 0;
IdentificationData* null = 0;
START_SECTION((IdentificationDataConverter()))
  ptr = new IdentificationData();
  TEST_NOT_EQUAL(ptr, null);
END_SECTION

START_SECTION((void importIDs(IdentificationData&, const vector<ProteinIdentification>&, const vector<PeptideIdentification>&)))
{
  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), proteins_in, peptides_in);
  // IdentificationData doesn't allow score types with the same name, but different orientations:
  peptides_in[0].setHigherScoreBetter(true);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  vector<ProteinIdentification> proteins_out;
  vector<PeptideIdentification> peptides_out;
  IdentificationDataConverter::exportIDs(ids, proteins_out, peptides_out);

  TEST_EQUAL(peptides_in.size(), peptides_out.size());
  TEST_EQUAL(proteins_in.size(), proteins_out.size());

  String filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out.idXML");
  // NEW_TMP_FILE(filename);
  // cout << "Test file:" << filename << endl;
  IdXMLFile().store(filename, proteins_out, peptides_out);

  MzTab mztab = IdentificationDataConverter::exportMzTab(ids);
  filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out.mzTab");
  // NEW_TMP_FILE(filename);
  // cout << "Test file:" << filename << endl;
  MzTabFile().store(filename, mztab);
}
END_SECTION

START_SECTION((void importSequences(IdentificationData&, const vector<FASTAFile::FASTAEntry>&, IdentificationData::MoleculeType, const String&)))
{
  vector<FASTAFile::FASTAEntry> fasta;
  FASTAFile::load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"), fasta);
  IdentificationData ids;
  IdentificationDataConverter::importSequences(ids, fasta);
  TEST_EQUAL(ids.getParentMolecules().size(), 5);
}
END_SECTION

START_SECTION((void exportIDs(const IdentificationData&, vector<ProteinIdentification>&, vector<PeptideIdentification>&)))
{
  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  PepXMLFile().load(OPENMS_GET_TEST_DATA_PATH("PepXMLFile_test_extended.pepxml"), proteins_in, peptides_in, "PepXMLFile_test");

  TEST_EQUAL(proteins_in.size(), 1);
  TEST_EQUAL(proteins_in[0].getHits().size(), 4);
  TEST_EQUAL(peptides_in.size(), 2);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  vector<ProteinIdentification> proteins_out;
  vector<PeptideIdentification> peptides_out;
  IdentificationDataConverter::exportIDs(ids, proteins_out, peptides_out);

  // the additional Peptide-/InterProphet scores cause multiple ID runs upon
  // export:
  TEST_EQUAL(proteins_out.size(), 3);
  TEST_EQUAL(peptides_out.size(), 4);

  String filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_pepXML_out.idXML");
  // NEW_TMP_FILE(filename);
  // cout << "Test file:" << filename << endl;
  IdXMLFile().store(filename, proteins_out, peptides_out);
}
END_SECTION

/*
// performance test on a large file:
START_SECTION(([[EXTRA]] void importIDs(IdentificationData&, const vector<ProteinIdentification>&, const vector<PeptideIdentification>&)))
{
  SysInfo::MemUsage mem_usage;
  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("large_test.idXML"), proteins_in, peptides_in);
  STATUS(mem_usage.delta("PeptideIdentification/ProteinIdentification"));

  TEST_EQUAL(proteins_in.size(), 1);
  TEST_EQUAL(proteins_in[0].getHits().size(), 11098);
  TEST_EQUAL(peptides_in.size(), 328591);
  TEST_EQUAL(proteins_in[0].getIndistinguishableProteins().size(), 10853);
  TEST_EQUAL(proteins_in[0].getProteinGroups().size(), 9092);

  mem_usage.reset();
  mem_usage.before();
  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);
  STATUS(mem_usage.delta("IdentificationData"));

  TEST_EQUAL(ids.getParentMolecules().size(), 11098);
  // problem: input data comes from multiple files, spectra with matching names
  // in different files get merged together -> lower number of data queries:
  TEST_EQUAL(ids.getDataQueries().size(), 55522);
  TEST_EQUAL(ids.getIdentifiedPeptides().size(), 73950);
  // according to "grep" on the input file, there should be 335250 peptide hits
  // in total - maybe some duplicates?:
  TEST_EQUAL(ids.getMoleculeQueryMatches().size(), 332778);

  TEST_EQUAL(ids.getParentMoleculeGroupings().size(), 2);
  TEST_EQUAL(ids.getParentMoleculeGroupings()[0].groups.size(), 10853);
  TEST_EQUAL(ids.getParentMoleculeGroupings()[1].groups.size(), 9092);
}
END_SECTION
*/

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
