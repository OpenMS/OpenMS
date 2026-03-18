// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(RNaseDigestion, "$Id$")

/////////////////////////////////////////////////////////////
RNaseDigestion* rd_ptr = 0;
RNaseDigestion* rd_null = 0;

START_SECTION(([EXTRA] RNaseDigestion()))
{
  rd_ptr = new RNaseDigestion;
  TEST_NOT_EQUAL(rd_ptr, rd_null);
}
END_SECTION

START_SECTION([EXTRA] ~RNaseDigestion())
{
  delete rd_ptr;
}
END_SECTION

START_SECTION((void setEnzyme(const String& enzyme_name)))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");
  TEST_EQUAL(rd.getEnzymeName(), "RNase_T1");
  rd.setEnzyme("cusativin");
  TEST_EQUAL(rd.getEnzymeName(), "cusativin");
  rd.setEnzyme("mazF");
  TEST_EQUAL(rd.getEnzymeName(), "mazF");
  rd.setEnzyme("colicin_E5");
  TEST_EQUAL(rd.getEnzymeName(), "colicin_E5");
  TEST_EXCEPTION(Exception::ElementNotFound, rd.setEnzyme("NoSuchEnzyme"));
}
END_SECTION

START_SECTION((void digest(const NASequence& rna, vector<NASequence>& output, Size min_length, Size max_length) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1"); // cuts after G and leaves a 3'-phosphate
  vector<NASequence> out;

  rd.digest(NASequence::fromString("AUC"), out);
  TEST_EQUAL(out.size(), 1);
  TEST_STRING_EQUAL(out[0].toString(), "AUC");
  out.clear();

  rd.digest(NASequence::fromString("AGUC"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "AGp");
  TEST_STRING_EQUAL(out[1].toString(), "UC");
  out.clear();

  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "pAUGp");
  TEST_STRING_EQUAL(out[1].toString(), "UCGp");
  TEST_STRING_EQUAL(out[2].toString(), "CAG");
  out.clear();

  // RNase T1 should cut after G, m2G and I, but not after m1G, m7G or Gm:
  rd.digest(NASequence::fromString("GI[m2G][m1G][m7G][Gm]A"), out);
  TEST_EQUAL(out.size(), 4);
  TEST_STRING_EQUAL(out[0].toString(), "Gp");
  TEST_STRING_EQUAL(out[1].toString(), "Ip");
  TEST_STRING_EQUAL(out[2].toString(), "[m2G]p");
  TEST_STRING_EQUAL(out[3].toString(), "[m1G][m7G][Gm]A");
  out.clear();

  // RNase T1 should not cut after m22G:
  rd.digest(NASequence::fromString("G[m2,2G]A"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "Gp");
  TEST_STRING_EQUAL(out[1].toString(), "[m2,2G]A");
  out.clear();

  rd.setMissedCleavages(2);
  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 6);
  TEST_STRING_EQUAL(out[0].toString(), "pAUGp");
  TEST_STRING_EQUAL(out[1].toString(), "pAUGUCGp");
  TEST_STRING_EQUAL(out[2].toString(), "pAUGUCGCAG");
  TEST_STRING_EQUAL(out[3].toString(), "UCGp");
  TEST_STRING_EQUAL(out[4].toString(), "UCGCAG");
  TEST_STRING_EQUAL(out[5].toString(), "CAG");
  out.clear();

  rd.setEnzyme("cusativin");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "CCCp");
  TEST_STRING_EQUAL(out[1].toString(), "AUCCp");
  TEST_STRING_EQUAL(out[2].toString(), "G");
  out.clear();

  rd.setEnzyme("cusativin");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "CCCp");
  TEST_STRING_EQUAL(out[1].toString(), "AUCCp");
  TEST_STRING_EQUAL(out[2].toString(), "G");
  out.clear();

  rd.setEnzyme("mazF");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("A[m6A]CA[m5C]AGGACGACAAAG"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "A[m6A]CA[m5C]AGGACGp");
  TEST_STRING_EQUAL(out[1].toString(), "ACAAAG");
  out.clear();


  rd.setEnzyme("colicin_E5");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("GGAUGUAAA"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "GGAUGp");  
  TEST_STRING_EQUAL(out[1].toString(), "UAAA");  
  out.clear();

  rd.setEnzyme("no cleavage");
  rd.setMissedCleavages(3);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 1);
  TEST_STRING_EQUAL(out[0].toString(), "CCCAUCCG");

  rd.setEnzyme("RNase_4c");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("UCCGUACGG"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "UCCGUc");    
  TEST_STRING_EQUAL(out[1].toString(), "ACGG");

  rd.setEnzyme("unspecific cleavage");
  rd.setMissedCleavages(0); // shouldn't matter for the result
  rd.digest(NASequence::fromString("ACGU"), out);
  TEST_EQUAL(out.size(), 10);
  TEST_STRING_EQUAL(out[0].toString(), "Ap");
  TEST_STRING_EQUAL(out[1].toString(), "ACp");
  TEST_STRING_EQUAL(out[2].toString(), "ACGp");
  TEST_STRING_EQUAL(out[3].toString(), "ACGU");
  TEST_STRING_EQUAL(out[4].toString(), "Cp");
  TEST_STRING_EQUAL(out[5].toString(), "CGp");
  TEST_STRING_EQUAL(out[6].toString(), "CGU");
  TEST_STRING_EQUAL(out[7].toString(), "Gp");
  TEST_STRING_EQUAL(out[8].toString(), "GU");
  TEST_STRING_EQUAL(out[9].toString(), "U");
}
END_SECTION

START_SECTION((void digest(IdentificationData& id_data, Size min_length = 0,
                Size max_length = 0) const))
{
  IdentificationData id_data;
  IdentificationData::ParentSequence rna("test", IdentificationData::MoleculeType::RNA, "pAUGUCGCAG");
  id_data.registerParentSequence(rna);

  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1"); // cuts after G and leaves a 3'-phosphate
  rd.digest(id_data);

  TEST_EQUAL(id_data.getIdentifiedOligos().size(), 3);

  /// multiple occurrences of the same oligo:
  IdentificationData id_data2;
  rna.sequence = "ACUGACUGG";
  id_data2.registerParentSequence(rna);

  rd.digest(id_data2, 2);

  TEST_EQUAL(id_data2.getIdentifiedOligos().size(), 1);
  ABORT_IF(id_data2.getIdentifiedOligos().empty());
  IdentificationData::IdentifiedOligoRef ref = id_data2.getIdentifiedOligos().begin();
  TEST_EQUAL(ref->parent_matches.size(), 1);
  ABORT_IF(ref->parent_matches.empty());
  // oligo sequence matches in two locations:
  const set<IdentificationData::ParentMatch>& matches =
    ref->parent_matches.begin()->second;
  TEST_EQUAL(matches.size(), 2);
  ABORT_IF(matches.size() < 2);
  auto match_it = matches.begin();
  TEST_EQUAL(match_it->start_pos, 0);
  ++match_it;
  TEST_EQUAL(match_it->start_pos, 4);

  // cleavage behavior with m2,2G should match sequence-level digest:
  IdentificationData id_data3;
  rna.sequence = "G[m2,2G]A";
  id_data3.registerParentSequence(rna);
  rd.digest(id_data3);
  TEST_EQUAL(id_data3.getIdentifiedOligos().size(), 2);
}
END_SECTION

START_SECTION((std::vector<std::pair<Size, Size>> getFragmentPositions(const NASequence& rna, Size min_length = 0, Size max_length = 0) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");

  vector<pair<Size, Size>> pos =
    rd.getFragmentPositions(NASequence::fromString("G[m2,2G]A"));
  TEST_EQUAL(pos.size(), 2);
  ABORT_IF(pos.size() < 2);
  TEST_EQUAL(pos[0].first, 0);
  TEST_EQUAL(pos[0].second, 1);
  TEST_EQUAL(pos[1].first, 1);
  TEST_EQUAL(pos[1].second, 2);
}
END_SECTION

START_SECTION((void digest(const NASequence& rna, std::vector<DigestionProduct>& output,
                Size min_length = 0, Size max_length = 0) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");

  vector<RNaseDigestion::DigestionProduct> out;
  rd.digest(NASequence::fromString("AGUC"), out);

  TEST_EQUAL(out.size(), 2);
  ABORT_IF(out.size() < 2);
  TEST_STRING_EQUAL(out[0].fragment.toString(), "AGp");
  TEST_EQUAL(out[0].position.first, 0);
  TEST_EQUAL(out[0].position.second, 2);
  TEST_STRING_EQUAL(out[1].fragment.toString(), "UC");
  TEST_EQUAL(out[1].position.first, 2);
  TEST_EQUAL(out[1].position.second, 2);
}
END_SECTION

START_SECTION((void digestWithCleavageSensitiveMods(const NASequence& rna,
                const CleavageSensitiveModGroups& cleavage_sensitive_mods,
                Size max_sensitive_mods_per_fragment,
                std::vector<DigestionProduct>& output,
                Size min_length = 0,
                Size max_length = 0) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");

  set<const Ribonucleotide*> variable_mods;
  variable_mods.insert(RibonucleotideDB::getInstance()->getRibonucleotide("m2,2G"));
  RNaseDigestion::CleavageSensitiveModGroups groups =
    rd.inferCleavageSensitiveMods(variable_mods);

  vector<RNaseDigestion::DigestionProduct> out;
  rd.digestWithCleavageSensitiveMods(NASequence::fromString("GGA"), groups, 1,
                                     out);

  set<String> sequences;
  for (const auto& product : out)
  {
    sequences.insert(product.fragment.toString());
  }

  TEST_EQUAL(sequences.count("Gp"), 1);
  TEST_EQUAL(sequences.count("A"), 1);
  TEST_EQUAL(sequences.count("[m2,2G]Gp"), 1);
  TEST_EQUAL(sequences.count("[m2,2G]A"), 1);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
