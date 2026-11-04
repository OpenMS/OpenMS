// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <tuple>
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

START_SECTION((void setEnzyme(const std::string& enzyme_name)))
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

  // sequence shorter than the minimum length: nothing to report (used to wrap around and run past the sequence)
  rd.digest(NASequence::fromString("ACGU"), out, 5);
  TEST_EQUAL(out.size(), 0);

  // min_length above max_length: nothing to report (used to wrap around and throw from reserve())
  rd.digest(NASequence::fromString("ACGUACGUACGUACGUACGU"), out, 10, 5);
  TEST_EQUAL(out.size(), 0);
}
END_SECTION

START_SECTION((owning digestion retains repeated evidence and target / decoy membership))
{
  IdentificationData data;
  auto& run = data.addRun("RNA", IdentificationData::MoleculeKind::OLIGONUCLEOTIDE);
  IdentificationData::Database database;
  database.path = "db";
  IdentificationData::DatabaseSequence parent;
  parent.database = run.addDatabase(database);
  parent.accession = "test";
  parent.sequence = "pAUGUCGCAG";
  parent.target_decoy = IdentificationData::TargetDecoy::TARGET;
  run.setDatabaseSequences(std::vector {parent});
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");
  auto output = rd.digest(run);
  TEST_EQUAL(output.size(), 3);
  parent.sequence = "ACUGACUGG";
  run.setDatabaseSequences(std::vector {parent});
  output = rd.digest(run, 2);
  TEST_EQUAL(output.size(), 1);
  ABORT_IF(output.empty());
  TEST_EQUAL(output[0].sequence.toString(), "ACUGp");
  TEST_EQUAL(output[0].sequence_evidence.size(), 2);
  TEST_EQUAL(output[0].sequence_evidence[0].start.value(), 0);
  TEST_EQUAL(output[0].sequence_evidence[1].start.value(), 4);
  TEST_EQUAL(output[0].sequence_evidence[0].accession, "test");
  TEST_EQUAL(output[0].target_decoy == IdentificationData::TargetDecoy::TARGET, true);
  auto decoy = parent;
  decoy.accession = "decoy";
  decoy.target_decoy = IdentificationData::TargetDecoy::DECOY;
  run.setDatabaseSequences(std::vector {parent, decoy});
  output = rd.digest(run, 2);
  TEST_EQUAL(output[0].sequence_evidence.size(), 4);
  TEST_EQUAL(output[0].target_decoy == IdentificationData::TargetDecoy::BOTH, true);

  // cleavage behavior with m2,2G should match sequence-level digest:
  parent.sequence = "G[m2,2G]A";
  run.setDatabaseSequences(std::vector {parent});
  output = rd.digest(run);
  TEST_EQUAL(output.size(), 2);
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

  set<string> sequences;
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


START_SECTION(([EXTRA] digestWithCleavageSensitiveMods agrees with explicit modified-parent digestion))
{
  RNaseDigestion digestion;
  digestion.setEnzyme("RNase_T1");
  const auto* mod = RibonucleotideDB::getInstance()->getRibonucleotide("m2,2G");
  const auto groups = digestion.inferCleavageSensitiveMods({mod});
  using Product = std::tuple<std::string, Size, Size>;
  for (const std::string parent_text : {"AG", "GGA", "AAGGAA", "GGGG", "AG[m2,2G]GAA"})
  {
    const auto original = NASequence::fromString(parent_text);
    for (Size missed = 0; missed <= 2; ++missed)
    {
      digestion.setMissedCleavages(missed);
      for (Size budget = 0; budget <= 2; ++budget)
      {
        for (Size maximum : {Size(0), Size(4)})
        {
          std::set<Product> expected, actual;
          const auto enumerate = [&](auto&& self, NASequence& parent, Size position, Size used) -> void
          {
            if (position == parent.size())
            {
              std::vector<RNaseDigestion::DigestionProduct> products;
              digestion.digest(parent, products, 1, maximum);
              for (const auto& p : products) expected.emplace(p.fragment.toString(), p.position.first, p.position.second);
              return;
            }
            self(self, parent, position + 1, used);
            if (used < budget && !parent[position]->isModified() && parent[position]->getCode() == "G")
            {
              const auto* residue = parent[position];
              parent.set(position, mod);
              self(self, parent, position + 1, used + 1);
              parent.set(position, residue);
            }
          };
          NASequence parent = original;
          enumerate(enumerate, parent, 0, 0);
          std::vector<RNaseDigestion::DigestionProduct> products;
          digestion.digestWithCleavageSensitiveMods(original, groups, budget, products, 1, maximum);
          for (const auto& p : products) actual.emplace(p.fragment.toString(), p.position.first, p.position.second);
          TEST_EQUAL(actual.size(), products.size())
          TEST_TRUE(actual == expected)
        }
      }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
