// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <string>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>


///////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

// must be defined up here or it won't compile:
struct IsEven
{
  typedef int argument_type;

  bool operator()(int i) const
  {
    return (i % 2 == 0);
  }
} is_even;

START_TEST(IDFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// load input data
// @TODO: use an example with more than one peptide ID
vector<ProteinIdentification> global_proteins;
vector<PeptideIdentification> global_peptides;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test.idXML"),
                 global_proteins, global_peptides);
global_peptides[0].sort(); // makes it easier to compare results

IDFilter* ptr = nullptr;
IDFilter* nullPointer = nullptr;

START_SECTION((IDFilter()))
  ptr = new IDFilter();
  TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION((~IDFilter()))
  delete ptr;
END_SECTION

START_SECTION((template <class Container, class Predicate> static void removeMatchingItems(Container& items, const Predicate& pred)))
{
  vector<int> numbers(6);
  for (Size i = 0; i < 6; ++i)
  {
    numbers[i] = i;
  }
  IDFilter::removeMatchingItems(numbers, is_even);
  TEST_EQUAL(numbers.size(), 3);
  TEST_EQUAL(numbers[0], 1);
  TEST_EQUAL(numbers[1], 3);
  TEST_EQUAL(numbers[2], 5);
}
END_SECTION

START_SECTION((template <class Container, class Predicate> static void keepMatchingItems(Container& items, const Predicate& pred)))
{
  vector<int> numbers(6);
  for (Size i = 0; i < 6; ++i)
  {
    numbers[i] = i;
  }
  IDFilter::keepMatchingItems(numbers, is_even);
  TEST_EQUAL(numbers.size(), 3);
  TEST_EQUAL(numbers[0], 0);
  TEST_EQUAL(numbers[1], 2);
  TEST_EQUAL(numbers[2], 4);
}
END_SECTION

START_SECTION((template <class IdentificationType> static Size countHits(const vector<IdentificationType>& ids)))
{
  vector<PeptideIdentification> peptides(4);
  peptides[0].getHits().resize(1);
  peptides[1].getHits().resize(3);
  // no hits in peptides[2]
  peptides[3].getHits().resize(2);

  TEST_EQUAL(IDFilter::countHits(peptides), 6);
}
END_SECTION

START_SECTION((template <class IdentificationType> static bool getBestHit(const vector<IdentificationType>& identifications, bool assume_sorted, typename IdentificationType::HitType& best_hit)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  PeptideHit best_hit;
  IDFilter::getBestHit(peptides, true, best_hit);
  TEST_REAL_SIMILAR(best_hit.getScore(), 40);
  TEST_EQUAL(best_hit.getSequence().toString(), "FINFGVNVEVLSRFQTK");

  peptides[0].setHigherScoreBetter(false);
  IDFilter::getBestHit(peptides, false, best_hit);
  TEST_REAL_SIMILAR(best_hit.getScore(), 10);
  TEST_EQUAL(best_hit.getSequence().toString(),
                    "MSLLSNM(Oxidation)ISIVKVGYNAR");
  ProteinHit best_hit2;
  IDFilter::getBestHit(global_proteins, false, best_hit2);
  TEST_REAL_SIMILAR(best_hit2.getScore(), 32.3);
  TEST_EQUAL(best_hit2.getAccession(), "Q824A5");
}
END_SECTION

START_SECTION((static void extractPeptideSequences(const vector<PeptideIdentification>& peptides, set<String>& sequences, bool ignore_mods = false)))
{
  set<String> seqs;
  IDFilter::extractPeptideSequences(global_peptides, seqs);
  TEST_EQUAL(seqs.size(), 11);
  vector<String> expected = ListUtils::create<String>("AITSDFANQAKTVLQNFK,DLEPGTDYEVTVSTLFGR,EGASTDFAALRTFLAEDGK,FINFGVNVEVLSRFQTK,LHASGITVTEIPVTATNFK,MRSLGYVAVISAVATDTDK,MSLLSNM(Oxidation)ISIVKVGYNAR,MSLLSNMISIVKVGYNAR,TGCDTWGQGTLVTVSSASTK,THPYGHAIVAGIERYPSK,TLCHHDATFDNLVWTPK");
  vector<String> expected_unmodified = ListUtils::create<String>("AITSDFANQAKTVLQNFK,DLEPGTDYEVTVSTLFGR,EGASTDFAALRTFLAEDGK,FINFGVNVEVLSRFQTK,LHASGITVTEIPVTATNFK,MRSLGYVAVISAVATDTDK,MSLLSNMISIVKVGYNAR,MSLLSNMISIVKVGYNAR,TGCDTWGQGTLVTVSSASTK,THPYGHAIVAGIERYPSK,TLCHHDATFDNLVWTPK");
  Size counter = 0;
  for (set<String>::iterator it = seqs.begin(); it != seqs.end(); ++it,
         ++counter)
  {
    TEST_EQUAL(*it, expected[counter]);
  }

  seqs.clear();
  IDFilter::extractPeptideSequences(global_peptides, seqs, true);
  TEST_EQUAL(seqs.size(), 10);
  counter = 0;
  for (set<String>::iterator it = seqs.begin(); it != seqs.end(); ++it,
         ++counter)
  {
    if (counter == 6) counter++; // skip the modified sequence
    TEST_EQUAL(*it, expected_unmodified[counter]);
  }
}
END_SECTION

START_SECTION((class PeptideDigestionFilter::operator(PeptideHit& hit)))
{
  ProteaseDigestion digestion;
  digestion.setEnzyme("Trypsin");
  
  IDFilter::PeptideDigestionFilter filter(digestion, 0, 1);
  vector<PeptideHit>hits, test_hits;

  
  // No cleavage
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("(MOD:00051)DFPIANGER")));
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("DFPIANGER")));
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("DFPIAN(Deamidated)GER")));

  // 1 - missed cleavage exception K before P
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("DFKPIARN(Deamidated)GER")));
  
  
  // 2 missed cleavages
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("(MOD:00051)DFPKIARNGER")));
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("DFPKIARNGER")));

  test_hits = hits;

  filter.filterPeptideSequences(test_hits);
  
  TEST_EQUAL(test_hits.size(), 4);
  for (UInt i = 0; i < test_hits.size(); i++)
  {
    TEST_EQUAL(test_hits[i].getSequence(), hits[i].getSequence());
  }

  IDFilter::PeptideDigestionFilter filter2(digestion, 0, 2);
  
  test_hits = hits;
  filter2.filterPeptideSequences(test_hits);
  
  TEST_EQUAL(test_hits.size(), hits.size());
  for (UInt i = 0; i < test_hits.size(); i++)
  {
    TEST_EQUAL(test_hits[i].getSequence(), hits[i].getSequence());
  }


  // Removing sequences
  hits.clear();
  hits.push_back(PeptideHit(0, 0, 0, AASequence::fromString("K(Dimethyl)FPIAUGR")));

  test_hits = hits;
  digestion.setEnzyme("Asp-N_ambic");
  
  //Should have exactly zero missed cleavages
  IDFilter::PeptideDigestionFilter filter3(digestion, 0, 0);

  filter3.filterPeptideSequences(test_hits);
  TEST_EQUAL(test_hits.size(), hits.size());
  for (UInt i = 0; i < test_hits.size(); i++)
  {
    TEST_EQUAL(test_hits[i].getSequence(), hits[i].getSequence());
  }

}
END_SECTION


START_SECTION((template <class IdentificationType> static void updateHitRanks(vector<IdentificationType>& ids)))
{
  TEST_EQUAL(global_peptides[0].getHits()[0].getRank(), 0);
  TEST_EQUAL(global_peptides[0].getHits()[1].getRank(), 0);
  TEST_EQUAL(global_peptides[0].getHits()[2].getRank(), 0);
  IDFilter::updateHitRanks(global_peptides);
  TEST_EQUAL(global_peptides[0].getHits()[0].getRank(), 1);
  TEST_EQUAL(global_peptides[0].getHits()[1].getRank(), 1);
  TEST_EQUAL(global_peptides[0].getHits()[2].getRank(), 2);

  TEST_EQUAL(global_proteins[0].getHits()[0].getRank(), 0);
  TEST_EQUAL(global_proteins[0].getHits()[1].getRank(), 0);
  TEST_EQUAL(global_proteins[0].getHits()[2].getRank(), 0);
  IDFilter::updateHitRanks(global_proteins);
  TEST_EQUAL(global_proteins[0].getHits()[0].getRank(), 1);
  TEST_EQUAL(global_proteins[0].getHits()[1].getRank(), 2);
  TEST_EQUAL(global_proteins[0].getHits()[2].getRank(), 3);
}
END_SECTION

START_SECTION((static void removeUnreferencedProteins(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides)))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test4.idXML"),
                   proteins, peptides);
  IDFilter::removeUnreferencedProteins(proteins, peptides);
  vector<ProteinHit>& hits = proteins[0].getHits();

  TEST_EQUAL(hits.size(), 3);
  TEST_EQUAL(hits[0].getAccession(), "Q824A5");
  TEST_EQUAL(hits[1].getAccession(), "S53854");
  TEST_EQUAL(hits[2].getAccession(), "Q872T5");
}
END_SECTION

START_SECTION((static void updateProteinReferences(vector<PeptideIdentification>& peptides, const vector<ProteinIdentification>& proteins, bool remove_peptides_without_reference = false)))
{
  vector<ProteinIdentification> proteins = global_proteins;
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();
  // create a peptide hit that matches to two proteins:
  peptide_hits[3].addPeptideEvidence(peptide_hits[4].getPeptideEvidences()[0]);
  TEST_EQUAL(peptide_hits[3].getPeptideEvidences().size(), 2);
  TEST_EQUAL(peptide_hits[4].getPeptideEvidences().size(), 1);
  proteins[0].getHits().resize(2);

  IDFilter::updateProteinReferences(peptides, proteins);
  TEST_EQUAL(peptide_hits.size(), 11);
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    if ((i == 3) || (i == 4))
    {
      TEST_EQUAL(peptide_hits[i].getPeptideEvidences().size(), 1);
      TEST_EQUAL(peptide_hits[i].getPeptideEvidences()[0].
                        getProteinAccession(), "Q824A5");
    }
    else
    {
      TEST_EQUAL(peptide_hits[i].getPeptideEvidences().size(), 0);
    }
  }

  // remove peptide hits without any reference to an existing proteins:
  IDFilter::updateProteinReferences(peptides, proteins, true);
  TEST_EQUAL(peptide_hits.size(), 2);
}
END_SECTION

START_SECTION((bool updateProteinGroups(vector<ProteinIdentification::ProteinGroup>& groups, const vector<ProteinHit>& hits)))
{
  vector<ProteinIdentification::ProteinGroup> groups(2);
  groups[0].accessions.push_back("A");
  groups[0].probability = 0.1;
  groups[1].accessions.push_back("B");
  groups[1].accessions.push_back("C");
  groups[1].probability = 0.2;

  vector<ProteinHit> hits(3);
  hits[0].setAccession("C");
  hits[1].setAccession("B");
  hits[2].setAccession("A");

  vector<ProteinIdentification::ProteinGroup> groups_copy = groups;

  // no protein to remove:
  bool valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, true);
  TEST_EQUAL(groups_copy.size(), 2);
  TEST_TRUE(groups_copy == groups);

  // remove full protein group:
  hits.pop_back();
  valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, true);
  TEST_EQUAL(groups_copy.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions.size(), 2);
  TEST_EQUAL(groups_copy[0].accessions[0], "B");
  TEST_EQUAL(groups_copy[0].accessions[1], "C");
  TEST_EQUAL(groups_copy[0].probability, 0.2);

  // remove part of a protein group:
  hits.pop_back();
  valid = IDFilter::updateProteinGroups(groups_copy, hits);
  TEST_EQUAL(valid, false);
  TEST_EQUAL(groups_copy.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions.size(), 1);
  TEST_EQUAL(groups_copy[0].accessions[0], "C");
  TEST_EQUAL(groups_copy[0].probability, 0.2);
}
END_SECTION

START_SECTION((template <class IdentificationType> static void removeEmptyIdentifications(vector<IdentificationType>& ids)))
{
  vector<ProteinIdentification> proteins(2);
  proteins[1].getHits().resize(1);
  IDFilter::removeEmptyIdentifications(proteins);
  TEST_EQUAL(proteins.size(), 1);
  TEST_EQUAL(proteins[0].getHits().size(), 1);

  vector<PeptideIdentification> peptides(2);
  peptides[0].getHits().resize(1);
  IDFilter::removeEmptyIdentifications(peptides);
  TEST_EQUAL(peptides.size(), 1);
  TEST_EQUAL(peptides[0].getHits().size(), 1);
}
END_SECTION

START_SECTION((template <class IdentificationType> static void filterHitsByScore(vector<IdentificationType>& ids, double threshold_score)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 11);

  IDFilter::filterHitsByScore(peptides, 33);
  TEST_EQUAL(peptide_hits.size(), 5);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(), 
                    "FINFGVNVEVLSRFQTK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MSLLSNMISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "THPYGHAIVAGIERYPSK");
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[3].getSequence().toString(),
                    "LHASGITVTEIPVTATNFK");
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[4].getSequence().toString(),
                    "MRSLGYVAVISAVATDTDK");

  IDFilter::filterHitsByScore(peptides, 41);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 0);
}
END_SECTION

START_SECTION((template <class IdentificationType> static void keepNBestHits(vector<IdentificationType>& ids, Size n)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();

  IDFilter::keepNBestHits(peptides, 3);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");

  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(), 
                    "FINFGVNVEVLSRFQTK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MSLLSNMISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "THPYGHAIVAGIERYPSK");
}
END_SECTION

START_SECTION((template <class IdentificationType> static void filterHitsByRank(vector<IdentificationType>& ids, Size min_rank, Size max_rank)))
{
  vector<ProteinIdentification> proteins = global_proteins;
  vector<PeptideIdentification> peptides = global_peptides;

  IDFilter::filterHitsByRank(peptides, 1, 5);
  TEST_EQUAL(peptides[0].getHits().size(), 6); // two rank 1 hits (same score)

  IDFilter::filterHitsByRank(proteins, 3, 10);
  TEST_EQUAL(proteins[0].getHits().size(), 2);
}
END_SECTION

START_SECTION((template <class IdentificationType> static void removeDecoyHits(vector<IdentificationType>& ids)))
{
  vector<ProteinIdentification> proteins(1);
  proteins[0].getHits().resize(5);
  proteins[0].getHits()[0].setMetaValue("target_decoy", "target");
  proteins[0].getHits()[1].setMetaValue("target_decoy", "decoy");
  // no meta value on hit 2
  proteins[0].getHits()[3].setMetaValue("isDecoy", "true");
  proteins[0].getHits()[4].setMetaValue("isDecoy", "false");
  IDFilter::removeDecoyHits(proteins);
  TEST_EQUAL(proteins[0].getHits().size(), 3);
  TEST_EQUAL(proteins[0].getHits()[0].getMetaValue("target_decoy"),
                    "target");
  TEST_EQUAL(proteins[0].getHits()[1].metaValueExists("target_decoy"), false);
  TEST_EQUAL(proteins[0].getHits()[1].metaValueExists("isDecoy"), false);
  TEST_EQUAL(proteins[0].getHits()[2].getMetaValue("isDecoy"), "false");

  vector<PeptideIdentification> peptides(1);
  peptides[0].getHits().resize(6);
  peptides[0].getHits()[0].setMetaValue("target_decoy", "target");
  peptides[0].getHits()[1].setMetaValue("target_decoy", "decoy");
  peptides[0].getHits()[2].setMetaValue("target_decoy", "target+decoy");
  // no meta value on hit 3
  peptides[0].getHits()[4].setMetaValue("isDecoy", "true");
  peptides[0].getHits()[5].setMetaValue("isDecoy", "false");
  IDFilter::removeDecoyHits(peptides);
  TEST_EQUAL(peptides[0].getHits().size(), 4);
  TEST_EQUAL(peptides[0].getHits()[0].getMetaValue("target_decoy"),
                    "target");
  TEST_EQUAL(peptides[0].getHits()[1].getMetaValue("target_decoy"),
                    "target+decoy");
  TEST_EQUAL(peptides[0].getHits()[2].metaValueExists("target_decoy"), false);
  TEST_EQUAL(peptides[0].getHits()[2].metaValueExists("isDecoy"), false);
  TEST_EQUAL(peptides[0].getHits()[3].getMetaValue("isDecoy"), "false");
}
END_SECTION

START_SECTION((template <class IdentificationType> static void removeHitsMatchingProteins(vector<IdentificationType>& ids, const set<String> accessions)))
{
  set<String> accessions;
  accessions.insert("Q824A5");
  accessions.insert("Q872T5");

  vector<ProteinIdentification> proteins = global_proteins;
  IDFilter::removeHitsMatchingProteins(proteins, accessions);

  TEST_EQUAL(proteins[0].getScoreType(), "Mascot");
  TEST_EQUAL(proteins[0].getHits().size(), 2);
  TEST_EQUAL(proteins[0].getHits()[0].getAccession(), "AAD30739");
  TEST_EQUAL(proteins[0].getHits()[1].getAccession(), "S53854");

  vector<PeptideIdentification> peptides = global_peptides;
  IDFilter::removeHitsMatchingProteins(peptides, accessions);

  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptides[0].getHits().size(), 9);
  // check some examples:
  TEST_EQUAL(peptides[0].getHits()[0].getSequence().toString(),
                    "FINFGVNVEVLSRFQTK");
  TEST_EQUAL(peptides[0].getHits()[3].getSequence().toString(),
                    "EGASTDFAALRTFLAEDGK");
  TEST_EQUAL(peptides[0].getHits()[8].getSequence().toString(),
                    "MSLLSNM(Oxidation)ISIVKVGYNAR");
}
END_SECTION

START_SECTION((template <class IdentificationType> static void keepHitsMatchingProteins(vector<IdentificationType>& ids, const set<String> accessions)))
{
  set<String> accessions;
  accessions.insert("Q824A5");
  accessions.insert("Q872T5");

  vector<ProteinIdentification> proteins = global_proteins;
  IDFilter::keepHitsMatchingProteins(proteins, accessions);

  TEST_EQUAL(proteins[0].getScoreType(), "Mascot");
  TEST_EQUAL(proteins[0].getHits().size(), 2);
  TEST_EQUAL(proteins[0].getHits()[0].getAccession(), "Q824A5");
  TEST_EQUAL(proteins[0].getHits()[1].getAccession(), "Q872T5");

  vector<PeptideIdentification> peptides = global_peptides;
  IDFilter::keepHitsMatchingProteins(peptides, accessions);

  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptides[0].getHits().size(), 2);
  TEST_EQUAL(peptides[0].getHits()[0].getSequence().toString(),
                    "LHASGITVTEIPVTATNFK");
  TEST_EQUAL(peptides[0].getHits()[1].getSequence().toString(),
                    "MRSLGYVAVISAVATDTDK");
}
END_SECTION

START_SECTION((static void keepBestPeptideHits(vector<PeptideIdentification>& peptides, bool strict = false)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();

  // not strict:
  IDFilter::keepBestPeptideHits(peptides);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(),
                    "FINFGVNVEVLSRFQTK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MSLLSNMISIVKVGYNAR");

  // strict:
  IDFilter::keepBestPeptideHits(peptides, true);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 0);
}
END_SECTION

START_SECTION((static void filterPeptidesByLength(vector<PeptideIdentification>& peptides, Size min_length, Size max_length = UINT_MAX)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  AASequence eighter = AASequence::fromString("OKTAMERR");
  AASequence niner = AASequence::fromString("NONAMERRR");
  AASequence tener = AASequence::fromString("DECAMERRRR");
  peptides[0].insertHit(PeptideHit(99.99, 1, 2, eighter));
  peptides[0].insertHit(PeptideHit(99.99, 1, 2, niner));
  peptides[0].insertHit(PeptideHit(99.99, 1, 2, tener));
  TEST_EQUAL(peptides[0].getHits().size(), 14);

  vector<PeptideIdentification> peptides2 = peptides;
  vector<PeptideHit>& peptide_hits = peptides2[0].getHits();
  IDFilter::filterPeptidesByLength(peptides2, 10);
  TEST_EQUAL(peptide_hits.size(), 12)
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 10, true);
  }

  peptides2 = peptides;
  IDFilter::filterPeptidesByLength(peptides2, 9, 10);
  TEST_EQUAL(peptide_hits.size(), 2);
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 9, true);
    TEST_EQUAL(peptide_hits[i].getSequence().size() <= 10, true);
  }

  peptides2 = peptides;
  IDFilter::filterPeptidesByLength(peptides2, 9, 8);
  TEST_EQUAL(peptide_hits.size(), 13)
  for (Size i = 0; i < peptide_hits.size(); ++i)
  {
    TEST_EQUAL(peptide_hits[i].getSequence().size() >= 9, true);
  }
}
END_SECTION

START_SECTION((static void filterPeptidesByCharge(vector<PeptideIdentification>& peptides, Size min_charge, Size max_charge)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& hits = peptides[0].getHits();
  hits[3].setCharge(3);
  hits[4].setCharge(4);
  hits[6].setCharge(3);
  hits[8].setCharge(1);
  hits[10].setCharge(5);

  IDFilter::filterPeptidesByCharge(peptides, 3, 4);
  TEST_EQUAL(hits.size(), 3);
  TEST_EQUAL(hits[0].getCharge(), 3);
  TEST_EQUAL(hits[1].getCharge(), 4);
  TEST_EQUAL(hits[2].getCharge(), 3);
}
END_SECTION

START_SECTION((static void filterPeptidesByRT(vector<PeptideIdentification>& peptides, double min_rt, double max_rt)))
{
  vector<PeptideIdentification> peptides(5);
  peptides[1].setRT(1);
  peptides[2].setRT(2);
  peptides[3].setRT(2.5);
  peptides[4].setRT(1.5);

  IDFilter::filterPeptidesByRT(peptides, 1.0, 1.9);
  TEST_EQUAL(peptides.size(), 2);
  TEST_EQUAL(peptides[0].getRT(), 1.0);
  TEST_EQUAL(peptides[1].getRT(), 1.5);
}
END_SECTION

START_SECTION((static void filterPeptidesByMZ(vector<PeptideIdentification>& peptides, double min_mz, double max_mz)))
{
  vector<PeptideIdentification> peptides(5);
  peptides[1].setMZ(111.1);
  peptides[2].setMZ(222.2);
  peptides[3].setMZ(225.5);
  peptides[4].setMZ(115.5);

  IDFilter::filterPeptidesByMZ(peptides, 112.0, 223.3);
  TEST_EQUAL(peptides.size(), 2);
  TEST_EQUAL(peptides[0].getMZ(), 222.2);
  TEST_EQUAL(peptides[1].getMZ(), 115.5);
}
END_SECTION

START_SECTION((static void filterPeptidesByMZError(vector<PeptideIdentification>& peptides, double mass_error, bool unit_ppm)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  peptides[0].setMZ(1000.0);
  IDFilter::filterPeptidesByMZError(peptides, 1, false); // in Da
  TEST_EQUAL(peptides[0].getHits().size(), 7);
  for (vector<PeptideHit>::iterator it = peptides[0].getHits().begin(); 
       it != peptides[0].getHits().end(); ++it)
  {
    double mz = it->getSequence().getMonoWeight(Residue::Full, 2) / 2.0;
    TEST_EQUAL((mz >= 999.0) && (mz <= 1001.0), true);
  }

  IDFilter::filterPeptidesByMZError(peptides, 100.0, true); // in PPM
  TEST_EQUAL(peptides[0].getHits().size(), 4);
}
END_SECTION

START_SECTION((static void filterPeptidesByRTPredictPValue(vector<PeptideIdentification>& peptides, const String& metavalue_key, double threshold = 0.05)))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;

  { // RT prediction:
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test2.idXML"), 
                     proteins, peptides);
    IDFilter::filterPeptidesByRTPredictPValue(peptides, "predicted_RT_p_value",
                                              0.08);
    vector<PeptideHit>& hits = peptides[0].getHits();

    TEST_EQUAL(hits.size(), 4);
    TEST_EQUAL(hits[0].getSequence().toString(), "LHASGITVTEIPVTATNFK");
    TEST_EQUAL(hits[1].getSequence().toString(), "DLEPGTDYEVTVSTLFGR");
    TEST_EQUAL(hits[2].getSequence().toString(), "FINFGVNVEVLSRFQTK");
    TEST_EQUAL(hits[3].getSequence().toString(), "MSLLSNMISIVKVGYNAR");
  }
  { // first dim. RT prediction:
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test3.idXML"),
                     proteins, peptides);
    IDFilter::filterPeptidesByRTPredictPValue(peptides,
                                              "predicted_RT_p_value_first_dim",
                                              0.08);
    vector<PeptideHit>& hits = peptides[0].getHits();

    TEST_EQUAL(hits.size(), 4);
    TEST_EQUAL(hits[0].getSequence().toString(), "LHASGITVTEIPVTATNFK");
    TEST_EQUAL(hits[1].getSequence().toString(), "DLEPGTDYEVTVSTLFGR");
    TEST_EQUAL(hits[2].getSequence().toString(), "FINFGVNVEVLSRFQTK");
    TEST_EQUAL(hits[3].getSequence().toString(), "MSLLSNMISIVKVGYNAR");
  }
}
END_SECTION

START_SECTION((static void removePeptidesWithMatchingModifications(vector<PeptideIdentification>& peptides, const set<String>& modifications)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  set<String> mods;
  mods.insert("Carbamidomethyl (C)"); // not present in the data
  IDFilter::removePeptidesWithMatchingModifications(peptides, mods);
  TEST_TRUE(peptides == global_peptides); // no changes

  mods.clear(); // filter any mod.
  IDFilter::removePeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 10);
  for (vector<PeptideHit>::iterator it = peptides[0].getHits().begin();
       it != peptides[0].getHits().end(); ++it)
  {
    TEST_EQUAL(it->getSequence().isModified(), false);
  }

  peptides = global_peptides;
  mods.insert("Oxidation (M)"); // present in the data
  IDFilter::removePeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 10);
  for (vector<PeptideHit>::iterator it = peptides[0].getHits().begin();
       it != peptides[0].getHits().end(); ++it)
  {
    TEST_EQUAL(it->getSequence().isModified(), false);
  }
}
END_SECTION

START_SECTION((static void removePeptidesWithMatchingRegEx(vector<PeptideIdentification>& peptides, const String& regex)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  String re{"[BJXZ]"};

  IDFilter::removePeptidesWithMatchingRegEx(peptides, re);
  TEST_TRUE(peptides == global_peptides); // no changes

  PeptideHit aaa_hit1;
  aaa_hit1.setSequence(AASequence::fromString("BBBBB"));
  PeptideHit aaa_hit2;
  aaa_hit2.setSequence(AASequence::fromString("JJJJJ"));
  PeptideHit aaa_hit3;
  aaa_hit3.setSequence(AASequence::fromString("XXXXX"));
  peptides[0].getHits().push_back(aaa_hit1);
  peptides[0].getHits().push_back(aaa_hit2);
  peptides[0].getHits().push_back(aaa_hit3);

  TEST_EQUAL(peptides == global_peptides, false); // added aaa peptides
  TEST_EQUAL(peptides[0].getHits().size(), 14);

  IDFilter::removePeptidesWithMatchingRegEx(peptides, re);
  /// aaa peptides should now be removed
  TEST_TRUE(peptides == global_peptides);
  TEST_EQUAL(peptides[0].getHits().size(), 11);
}
END_SECTION

START_SECTION((static void keepPeptidesWithMatchingModifications(vector<PeptideIdentification>& peptides, const set<String>& modifications)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  set<String> mods;
  mods.insert("Oxidation (M)");
  IDFilter::keepPeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 1);
  TEST_EQUAL(peptides[0].getHits()[0].getSequence().toString(),
                    "MSLLSNM(Oxidation)ISIVKVGYNAR");

  // terminal mods:
  AASequence seq = AASequence::fromString("(Acetyl)PEPTIDER.(Arg-loss)");
  peptides[0].getHits().resize(2);
  peptides[0].getHits()[1].setSequence(seq);
  mods.insert("Acetyl (N-term)");
  IDFilter::keepPeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 2);

  mods.clear();
  mods.insert("Arg-loss (C-term R)");
  IDFilter::keepPeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 1);

  // mod. not present in the data:
  mods.clear();
  mods.insert("Carbamidomethyl (C)");
  IDFilter::keepPeptidesWithMatchingModifications(peptides, mods);
  TEST_EQUAL(peptides[0].getHits().size(), 0);
}
END_SECTION

START_SECTION((static void removePeptidesWithMatchingSequences(vector<PeptideIdentification>& peptides, const vector<PeptideIdentification>& bad_peptides, bool ignore_mods = false)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();
  vector<PeptideIdentification> bad_peptides(1);
  vector<PeptideHit>& bad_hits = bad_peptides[0].getHits();
  bad_hits.resize(8);
  bad_hits[0].setSequence(AASequence::fromString("LHASGITVTEIPVTATNFK"));
  bad_hits[1].setSequence(AASequence::fromString("MRSLGYVAVISAVATDTDK"));
  bad_hits[2].setSequence(AASequence::fromString("EGASTDFAALRTFLAEDGK"));
  bad_hits[3].setSequence(AASequence::fromString("DLEPGTDYEVTVSTLFGR"));
  bad_hits[4].setSequence(AASequence::fromString("FINFGVNVEVLSRFQTK"));
  bad_hits[5].setSequence(AASequence::fromString("MSLLSNMISIVKVGYNAR"));
  bad_hits[6].setSequence(AASequence::fromString("THPYGHAIVAGIERYPSK"));
  bad_hits[7].setSequence(AASequence::fromString("AITSDFANQAKTVLQNFK"));

  // modification-aware filtering:
  IDFilter::removePeptidesWithMatchingSequences(peptides, bad_peptides, false);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_EQUAL(peptide_hits[0].getSequence(),
             AASequence::fromString("TGCDTWGQGTLVTVSSASTK"));
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[1].getSequence(),
             AASequence::fromString("TLCHHDATFDNLVWTPK"));
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.37);
  TEST_EQUAL(peptide_hits[2].getSequence(),
             AASequence::fromString("MSLLSNM(Oxidation)ISIVKVGYNAR"));
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 10);

  // modification-unaware filtering:
  IDFilter::removePeptidesWithMatchingSequences(peptides, bad_peptides, true);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(), 
                    "TGCDTWGQGTLVTVSSASTK");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(), 
                    "TLCHHDATFDNLVWTPK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.37);
}
END_SECTION

START_SECTION((static void keepPeptidesWithMatchingSequences(vector<PeptideIdentification>& peptides, const vector<PeptideIdentification>& good_peptides, bool ignore_mods = false)))
{
  vector<PeptideIdentification> peptides = global_peptides;
  vector<PeptideHit>& peptide_hits = peptides[0].getHits();
  vector<PeptideIdentification> good_peptides(1);
  vector<PeptideHit>& good_hits = good_peptides[0].getHits();
  good_hits.resize(3);
  good_hits[0].setSequence(AASequence::fromString("TGCDTWGQGTLVTVSSASTK"));
  good_hits[1].setSequence(AASequence::fromString("TLCHHDATFDNLVWTPK"));
  good_hits[2].setSequence(AASequence::fromString("MSLLSNM(Oxidation)ISIVKVGYNAR"));

  // modification-unaware filtering:
  IDFilter::keepPeptidesWithMatchingSequences(peptides, good_peptides, true);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 4);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(), 
                    "MSLLSNMISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "TGCDTWGQGTLVTVSSASTK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "TLCHHDATFDNLVWTPK");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 10.37);
  TEST_EQUAL(peptide_hits[3].getSequence().toString(),
                    "MSLLSNM(Oxidation)ISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 10);

  // modification-aware filtering:
  IDFilter::keepPeptidesWithMatchingSequences(peptides, good_peptides, false);
  TEST_EQUAL(peptides[0].getScoreType(), "Mascot");
  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(),
                    "TGCDTWGQGTLVTVSSASTK");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 10.93);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "TLCHHDATFDNLVWTPK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 10.37);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "MSLLSNM(Oxidation)ISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 10);
}
END_SECTION

START_SECTION((static void keepUniquePeptidesPerProtein(vector<PeptideIdentification>& peptides)))
{
  vector<PeptideIdentification> peptides(1);
  vector<PeptideHit>& hits = peptides[0].getHits();
  hits.resize(4);
  hits[0].setMetaValue("protein_references", "non-unique");
  hits[1].setMetaValue("protein_references", "unmatched");
  // no meta value for hit 2
  hits[3].setMetaValue("protein_references", "unique");
  IDFilter::keepUniquePeptidesPerProtein(peptides);
  TEST_EQUAL(hits.size(), 1);
  TEST_EQUAL(hits[0].getMetaValue("protein_references"), "unique");
}
END_SECTION

START_SECTION((static void removeDuplicatePeptideHits(vector<PeptideIdentification>& peptides, bool seq_only)))
{
  vector<PeptideIdentification> peptides(1, global_peptides[0]);
  vector<PeptideHit>& hits = peptides[0].getHits();
  hits.clear();
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("DFPIANGER"));
  hit.setCharge(1);
  hit.setScore(0.3);
  hits.push_back(hit);
  hit.setCharge(2);
  hits.push_back(hit);
  hit.setScore(0.5);
  hits.push_back(hit);
  hit.setSequence(AASequence::fromString("DFPIANGEK"));
  hits.push_back(hit);
  hits.push_back(hit);
  hits.push_back(hit);
  hit.setCharge(5);
  hits.push_back(hit);
  TEST_EQUAL(hits.size(), 7);

  IDFilter::removeDuplicatePeptideHits(peptides);
  TEST_EQUAL(hits.size(), 5);
  TEST_EQUAL(hits[3].getSequence().toString(), "DFPIANGEK");
  TEST_EQUAL(hits[3].getCharge(), 2);
  TEST_EQUAL(hits[4].getSequence().toString(), "DFPIANGEK");
  TEST_EQUAL(hits[4].getCharge(), 5);

  IDFilter::removeDuplicatePeptideHits(peptides, true);
  TEST_EQUAL(hits.size(), 2);
  TEST_EQUAL(hits[0].getSequence().toString(), "DFPIANGER");
  TEST_EQUAL(hits[0].getScore(), 0.3);
  TEST_EQUAL(hits[1].getSequence().toString(), "DFPIANGEK");
}
END_SECTION

START_SECTION((template <class PeakT> static void filterHitsByScore(MSExperiment<PeakT>& experiment, double peptide_threshold_score, double protein_threshold_score)))
{
  PeakMap experiment;
  vector<PeptideIdentification> ids(1, global_peptides[0]);

  ids[0].assignRanks();

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(ids);

  IDFilter::filterHitsByScore(experiment, 31.8621, 0);
  PeptideIdentification& identification = experiment[3].getPeptideIdentifications()[0];
  TEST_EQUAL(identification.getScoreType(), "Mascot");

  vector<PeptideHit>& peptide_hits = identification.getHits();
  TEST_EQUAL(peptide_hits.size(), 5);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(),
                    "FINFGVNVEVLSRFQTK");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getRank(), 1);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MSLLSNMISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getRank(), 1);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "THPYGHAIVAGIERYPSK");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getRank(), 2);
  TEST_EQUAL(peptide_hits[3].getSequence().toString(),
                    "LHASGITVTEIPVTATNFK");
  TEST_REAL_SIMILAR(peptide_hits[3].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[3].getRank(), 3);
  TEST_EQUAL(peptide_hits[4].getSequence().toString(),
                    "MRSLGYVAVISAVATDTDK");
  TEST_REAL_SIMILAR(peptide_hits[4].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[4].getRank(), 4);
}
END_SECTION

START_SECTION((template <class PeakT> static void keepNBestHits(MSExperiment<PeakT>& experiment, Size n)))
{
  PeakMap experiment;
  vector<PeptideIdentification> ids(1, global_peptides[0]);

  ids[0].assignRanks();

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(ids);

  IDFilter::keepNBestHits(experiment, 3);
  PeptideIdentification& identification = experiment[3].getPeptideIdentifications()[0];
  TEST_EQUAL(identification.getScoreType(), "Mascot");

  vector<PeptideHit>& peptide_hits = identification.getHits();
  TEST_EQUAL(peptide_hits.size(), 3);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(),
                    "FINFGVNVEVLSRFQTK");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 40);
  TEST_EQUAL(peptide_hits[0].getRank(), 1);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MSLLSNMISIVKVGYNAR");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 40);
  TEST_EQUAL(peptide_hits[1].getRank(), 1);
  TEST_EQUAL(peptide_hits[2].getSequence().toString(),
                    "THPYGHAIVAGIERYPSK");
  TEST_REAL_SIMILAR(peptide_hits[2].getScore(), 39);
  TEST_EQUAL(peptide_hits[2].getRank(), 2);
}
END_SECTION

START_SECTION((static void keepNBestSpectra(std::vector<PeptideIdentification>& peptides, Size n)))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDFilter_test5.idXML"),
                   proteins, peptides);

  cout << peptides[0].getHits()[0].getSequence().toString() << endl;
  cout << peptides[1].getHits()[0].getSequence().toString() << endl;

  IDFilter::keepNBestSpectra(peptides, 2); // keep best two spectra (those with best hits)

  TEST_EQUAL(peptides.size(), 2);

  vector<PeptideHit> peptide_hits = peptides[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 2);

  peptide_hits = peptides[1].getHits();
  TEST_EQUAL(peptide_hits.size(), 2);

  cout << peptides[0].getHits()[0].getSequence().toString() << endl;
  cout << peptides[1].getHits()[0].getSequence().toString() << endl;
  TEST_REAL_SIMILAR(peptides[0].getHits()[0].getScore(), 1000);
  TEST_REAL_SIMILAR(peptides[1].getHits()[0].getScore(), 40);
}
END_SECTION

START_SECTION((template<class PeakT> static void keepHitsMatchingProteins(MSExperiment<PeakT>& experiment, const vector<FASTAFile::FASTAEntry>& proteins)))
{
  PeakMap experiment;
  vector<FASTAFile::FASTAEntry> proteins;
  vector<PeptideIdentification> peptides = global_peptides;

  proteins.push_back(FASTAFile::FASTAEntry("Q824A5", "first desription",
                                           "LHASGITVTEIPVTATNFK"));
  proteins.push_back(FASTAFile::FASTAEntry("Q872T5", "second description",
                                           "THPYGHAIVAGIERYPSK"));

  for (Size i = 0; i < 5; ++i)
  {
    experiment.addSpectrum(MSSpectrum());
  }
  experiment[3].setMSLevel(2);
  experiment[3].setPeptideIdentifications(peptides);

  IDFilter::keepHitsMatchingProteins(experiment, proteins);
  TEST_EQUAL(experiment[3].getPeptideIdentifications()[0].getScoreType(),
             "Mascot");

  vector<PeptideHit>& peptide_hits =
    experiment[3].getPeptideIdentifications()[0].getHits();
  TEST_EQUAL(peptide_hits.size(), 2);
  TEST_EQUAL(peptide_hits[0].getSequence().toString(), 
                    "LHASGITVTEIPVTATNFK");
  TEST_REAL_SIMILAR(peptide_hits[0].getScore(), 34.85);
  TEST_EQUAL(peptide_hits[0].getRank(), 1);
  TEST_EQUAL(peptide_hits[1].getSequence().toString(),
                    "MRSLGYVAVISAVATDTDK");
  TEST_REAL_SIMILAR(peptide_hits[1].getScore(), 33.85);
  TEST_EQUAL(peptide_hits[1].getRank(), 2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#pragma clang diagnostic pop
