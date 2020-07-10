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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

struct ComparePIdSize
{
      bool operator()(const ProteinIdentification& lhs, const ProteinIdentification& rhs) const
      {
        return lhs.getHits().size() < rhs.getHits().size();
      }
};

START_TEST(IdentificationDataConverter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

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
  vector<PeptideHit> hits_in, hits_out;
  for (const auto& pep : peptides_in)
  {
    hits_in.insert(hits_in.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (const auto& pep : peptides_out)
  {
    hits_out.insert(hits_out.end(), pep.getHits().begin(), pep.getHits().end());
  }
  TEST_EQUAL(hits_in.size(), hits_out.size());
  // order of hits is different, check that every output one is in the input:
  for (const auto& hit : hits_out)
  {
    TEST_EQUAL(find(hits_in.begin(), hits_in.end(), hit) != hits_in.end(),
               true);
  }

  std::sort(proteins_in.begin(), proteins_in.end(), ComparePIdSize());
  std::sort(proteins_out.begin(), proteins_out.end(), ComparePIdSize());
  TEST_EQUAL(proteins_in.size(), proteins_out.size());
  TEST_EQUAL(proteins_in[0].getHits().size(), 1) // is sorted
  TEST_EQUAL(proteins_in[1].getHits().size(), 2) // is sorted

  // the exporter adds target/decoy information (default: target):
  for (auto& hit : proteins_in[0].getHits()) hit.setMetaValue("target_decoy", "target");
  for (auto& hit : proteins_in[1].getHits()) hit.setMetaValue("target_decoy", "target");

  // TEST_EQUAL(proteins_in[0].getIdentifier(), proteins_out[0].getIdentifier() ) // identifiers are not equal
  // TEST_EQUAL(proteins_in[1].getIdentifier(), proteins_out[1].getIdentifier() ) // identifiers are not equal

  TEST_EQUAL(proteins_in[0].getHits().size(), proteins_out[0].getHits().size());
  TEST_EQUAL(proteins_in[1].getHits().size(), proteins_out[1].getHits().size());
  TEST_EQUAL(proteins_in[0].getHits() == proteins_out[0].getHits(), true);
  TEST_EQUAL(proteins_in[1].getHits() == proteins_out[1].getHits(), true);

  TEST_EQUAL(proteins_in[0].getDateTime().get(),
             proteins_out[0].getDateTime().get());
  TEST_EQUAL(proteins_in[1].getDateTime().get(),
             proteins_out[1].getDateTime().get());

  TEST_EQUAL(proteins_in[0].getSearchParameters() == proteins_out[0].getSearchParameters(), true)
  TEST_EQUAL(proteins_in[1].getSearchParameters() == proteins_out[1].getSearchParameters(), true)

  // String filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out.idXML");
  // IdXMLFile().store(filename, proteins_out, peptides_out);
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

  String filename = OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/FidoAdapter_4_output.idXML");
  IdXMLFile().load(filename, proteins_in, peptides_in);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  vector<ProteinIdentification> proteins_out;
  vector<PeptideIdentification> peptides_out;
  IdentificationDataConverter::exportIDs(ids, proteins_out, peptides_out);

  TEST_EQUAL(proteins_in.size(), proteins_out.size());
  TEST_EQUAL(proteins_in[0].getHits().size(),
             proteins_out[0].getHits().size());
  TEST_EQUAL(proteins_in[0].getHits() == proteins_out[0].getHits(), true);

  TEST_EQUAL(proteins_in[0].getIndistinguishableProteins() ==
             proteins_out[0].getIndistinguishableProteins(), true);
  TEST_EQUAL(proteins_in[0].getProteinGroups() ==
             proteins_out[0].getProteinGroups(), true);

  TEST_EQUAL(peptides_in.size(), peptides_out.size());
  // no "operator<" for PeptideHit, otherwise we could use a set:
  vector<PeptideHit> hits_in, hits_out;
  for (const auto& pep : peptides_in)
  {
    hits_in.insert(hits_in.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (const auto& pep : peptides_out)
  {
    hits_out.insert(hits_out.end(), pep.getHits().begin(), pep.getHits().end());
  }
  for (auto& hit : hits_in)
  {
    // "target+decoy" is counted as "target" in IdentificationData:
    if (hit.getMetaValue("target_decoy") == "target+decoy")
    {
      hit.setMetaValue("target_decoy", "target");
    }
  }
  TEST_EQUAL(hits_in.size(), hits_out.size());
  // order of hits is different, check that every output one is in the input:
  for (const auto& hit : hits_out)
  {
    TEST_EQUAL(find(hits_in.begin(), hits_in.end(), hit) != hits_in.end(),
               true);
  }

  // filename = OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out2.idXML");
  // IdXMLFile().store(filename, proteins_out, peptides_out);
}
END_SECTION

START_SECTION((MzTab exportMzTab(const IdentificationData& id_data)))
{
  vector<ProteinIdentification> proteins_in;
  vector<PeptideIdentification> peptides_in;
  String filename = OPENMS_GET_TEST_DATA_PATH("../../../topp/THIRDPARTY/FidoAdapter_4_output.idXML");
  IdXMLFile().load(filename, proteins_in, peptides_in);

  IdentificationData ids;
  IdentificationDataConverter::importIDs(ids, proteins_in, peptides_in);

  MzTab mztab = IdentificationDataConverter::exportMzTab(ids);
  NEW_TMP_FILE(filename);
  MzTabFile().store(filename, mztab);

  TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out1.mzTab"));

  // RNA data, oligonucleotide that matches several times in the same RNA:
  IdentificationData rna_ids;
  IdentificationData::ParentMolecule rna("test", IdentificationData::MoleculeType::RNA, "AUCGAUCG");
  IdentificationData::ParentMoleculeRef ref = rna_ids.registerParentMolecule(rna);
  IdentificationData::IdentifiedOligo oli(NASequence::fromString("AUCG"));
  IdentificationData::MoleculeParentMatch match1(0, 3), match2(4, 7);
  oli.parent_matches[ref].insert(match1);
  oli.parent_matches[ref].insert(match2);
  rna_ids.registerIdentifiedOligo(oli);

  mztab = IdentificationDataConverter::exportMzTab(rna_ids);
  NEW_TMP_FILE(filename);
  MzTabFile().store(filename, mztab);

  TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("IdentificationDataConverter_out2.mzTab"));
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
