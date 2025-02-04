// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/DBSuitability.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/regex.hpp>

int countAS(const vector<FASTAFile::FASTAEntry>& fasta)
{
  int counter = 0;

  for (const auto& entry : fasta)
  {
    counter += entry.sequence.size();
  }

  return counter;
}

START_TEST(Suitability, "$Id$")

/////////////////////////////////////////////////////////////
////////////////////// CREATE DATA //////////////////////////
/////////////////////////////////////////////////////////////

PeptideEvidence decoy_protein("DECOY_PROT", 0, 0, 'A', 'A');
PeptideEvidence target_protein("DB_PROT", 0, 0, 'A', 'A');
PeptideEvidence novo_protein(Constants::UserParam::CONCAT_PEPTIDE, 0, 0, 'A', 'A');

// target-db hits with different q-values
PeptideHit target_db_hit1;
target_db_hit1.setSequence(AASequence::fromString("PEP"));
target_db_hit1.setPeptideEvidences({ target_protein });
target_db_hit1.setMetaValue("target_decoy", "target+decoy");
target_db_hit1.setMetaValue("MS:1002252", 0.8);
target_db_hit1.setScore(0.002);

PeptideHit target_db_hit2;
target_db_hit2.setSequence(AASequence::fromString("PEP"));
target_db_hit2.setPeptideEvidences({ target_protein });
target_db_hit2.setMetaValue("target_decoy", "target");
target_db_hit2.setMetaValue("MS:1002252", 0.8);
target_db_hit2.setScore(0.011);

// target-novo hits with different xcorr scores
PeptideHit target_novo_hit1;
target_novo_hit1.setSequence(AASequence::fromString("PEP"));
target_novo_hit1.setPeptideEvidences({ novo_protein });
target_novo_hit1.setMetaValue("target_decoy", "target");
target_novo_hit1.setMetaValue("MS:1002252", 0.85); // diff to db 0.05
target_novo_hit1.setScore(0.001);

PeptideHit target_novo_hit2;
target_novo_hit2.setSequence(AASequence::fromString("PEP"));
target_novo_hit2.setPeptideEvidences({ novo_protein });
target_novo_hit2.setMetaValue("target_decoy", "target");
target_novo_hit2.setMetaValue("MS:1002252", 0.93); // diff to db 0.13
target_novo_hit2.setScore(0.001);

// decoy hits with different xcorr scores -> resulting cut-offs: 0.15, 0.1, 0.05 (devided by 377.3 - weight of "PEP")
PeptideHit decoy1;
decoy1.setSequence(AASequence::fromString("PEP"));
decoy1.setPeptideEvidences({ decoy_protein });
decoy1.setMetaValue("target_decoy", "decoy");
decoy1.setMetaValue("MS:1002252", 0.7);
decoy1.setScore(1);

PeptideHit decoy2;
decoy2.setSequence(AASequence::fromString("PEP"));
decoy2.setPeptideEvidences({ decoy_protein });
decoy2.setMetaValue("target_decoy", "decoy");
decoy2.setMetaValue("MS:1002252", 0.6);
decoy2.setScore(1);

PeptideHit decoy3;
decoy3.setSequence(AASequence::fromString("PEP"));
decoy3.setMetaValue("target_decoy", "decoy");
decoy3.setMetaValue("MS:1002252", 0.55);
decoy3.setScore(1);

PeptideHit high_decoy;
high_decoy.setSequence(AASequence::fromString("PEP"));
high_decoy.setMetaValue("target_decoy", "decoy");
high_decoy.setMetaValue("MS:1002252", 0.55);
high_decoy.setScore(0);

// some error throwing hits

PeptideHit no_xcorr_hit;
no_xcorr_hit.setSequence(AASequence::fromString("PEP"));
no_xcorr_hit.setPeptideEvidences({ decoy_protein });
no_xcorr_hit.setMetaValue("target_decoy", "decoy");
no_xcorr_hit.setScore(1);

// build identifications
vector<PeptideIdentification> pep_ids;
vector<PeptideIdentification> top_decoy;
vector<PeptideIdentification> few_decoys;
vector<PeptideIdentification> no_xcorr_ids;
PeptideIdentification pep_id;
pep_id.setScoreType("some_score");
pep_id.setHigherScoreBetter(false);
pep_id.setHits({ target_novo_hit1, decoy1, decoy2 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ target_db_hit1, decoy1, decoy3 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ target_db_hit2 });
pep_ids.push_back(pep_id);
few_decoys.push_back(pep_id);

pep_id.setHits({ high_decoy, target_db_hit2 });
top_decoy.push_back(pep_id);

pep_id.setHits({ target_novo_hit1, target_db_hit1, decoy2, decoy3});
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);
no_xcorr_ids.push_back(pep_id);

pep_id.setHits({ target_novo_hit2, target_db_hit1 });
pep_ids.push_back(pep_id);
top_decoy.push_back(pep_id);

pep_id.setHits({ no_xcorr_hit });
no_xcorr_ids.push_back(pep_id);

vector<PeptideIdentification> pep_ids_2(pep_ids);
vector<PeptideIdentification> pep_ids_3(pep_ids);

vector<PeptideIdentification> FDR_id;
pep_id.setScoreType("q-value");
pep_id.setHits({ decoy1 });
FDR_id.push_back(pep_id);

vector<FASTAFile::FASTAEntry> empty_fasta;
MSExperiment empty_exp;
ProteinIdentification::SearchParameters empty_params;

/////////////////////////////////////////////////////////////
///////////////////// START TESTING /////////////////////////
/////////////////////////////////////////////////////////////

DBSuitability* ptr = nullptr;
DBSuitability* nulpt = nullptr;
START_SECTION(DBSuitability())
{
  ptr = new DBSuitability();
  TEST_NOT_EQUAL(ptr, nulpt)
}
END_SECTION

START_SECTION(~DBSuitability())
{
  delete ptr;
}
END_SECTION

START_SECTION(void compute(std::vector<PeptideIdentification>&& pep_ids, const MSExperiment& exp, const std::vector<FASTAFile::FASTAEntry>& original_fasta, const std::vector<FASTAFile::FASTAEntry>& novo_fasta, const ProteinIdentification::SearchParameters& search_params))
{
  // Test normal suitability (without correction)
  DBSuitability s;
  Param p;
  p.setValue("disable_correction", "true");
  p.setValue("reranking_cutoff_percentile", 1.);
  s.setParameters(p);
  s.compute(move(pep_ids), empty_exp, empty_fasta, empty_fasta, empty_params);

  p.setValue("reranking_cutoff_percentile", 1./3);
  p.setValue("FDR", 0.);
  s.setParameters(p);
  s.compute(move(pep_ids_2), empty_exp, empty_fasta, empty_fasta, empty_params);
  s.compute(move(top_decoy), empty_exp, empty_fasta, empty_fasta, empty_params);

  p.setValue("reranking_cutoff_percentile", 0.);
  s.setParameters(p);
  s.compute(move(pep_ids_3), empty_exp, empty_fasta, empty_fasta, empty_params);
  vector<DBSuitability::SuitabilityData> d = s.getResults();
  DBSuitability::SuitabilityData data_fract_1 = d[0];
  DBSuitability::SuitabilityData data_fract_05 = d[1];
  DBSuitability::SuitabilityData data_decoy_top = d[2];
  DBSuitability::SuitabilityData data_small_percentile = d[3];
  TEST_REAL_SIMILAR(data_fract_1.cut_off, 0.00044);
  TEST_REAL_SIMILAR(data_fract_05.cut_off, 0.00029);
  TEST_REAL_SIMILAR(data_decoy_top.cut_off, 0.00029);
  TEST_REAL_SIMILAR(data_small_percentile.cut_off, 0.00014);
  TEST_EQUAL(data_fract_1.num_interest, 2);
  TEST_EQUAL(data_fract_05.num_interest, 2);
  TEST_EQUAL(data_decoy_top.num_interest, 0);
  TEST_EQUAL(data_small_percentile.num_interest, 2);
  TEST_EQUAL(data_fract_1.num_re_ranked, 2);
  TEST_EQUAL(data_fract_05.num_re_ranked, 1);
  TEST_EQUAL(data_decoy_top.num_re_ranked, 0);
  TEST_EQUAL(data_small_percentile.num_re_ranked, 0);
  TEST_EQUAL(data_fract_1.num_top_db, 4);
  TEST_EQUAL(data_fract_05.num_top_db, 3);
  TEST_EQUAL(data_decoy_top.num_top_db, 0);
  TEST_EQUAL(data_small_percentile.num_top_db, 2);
  TEST_EQUAL(data_fract_1.num_top_novo, 1);
  TEST_EQUAL(data_fract_05.num_top_novo, 2);
  TEST_EQUAL(data_decoy_top.num_top_novo, 0);
  TEST_EQUAL(data_small_percentile.num_top_novo, 3);
  TEST_REAL_SIMILAR(data_fract_1.suitability, 4./5);
  TEST_REAL_SIMILAR(data_fract_05.suitability, 3./5);
  TEST_REAL_SIMILAR(data_small_percentile.suitability, 2./5);
  TEST_EQUAL(data_decoy_top.suitability, DBL_MAX);

  TEST_EXCEPTION_WITH_MESSAGE(Exception::Precondition, s.compute(move(FDR_id), empty_exp, empty_fasta, empty_fasta, empty_params), "q-value found at PeptideIdentifications. That is not allowed! Please make sure FDR did not run previously.");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, s.compute(move(few_decoys), empty_exp, empty_fasta, empty_fasta, empty_params), "Under 20 % of peptide identifications have two decoy hits. This is not enough for re-ranking. Use the 'no_rerank' flag to still compute a suitability score.");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, s.compute(move(no_xcorr_ids), empty_exp, empty_fasta, empty_fasta, empty_params), "No cross correlation score found at peptide hit. Only Comet search engine is supported for re-ranking. Set 'force' flag to use the default score for this. This may result in undefined behaviour and is not advised.");

  // Corrected Suitability is to complicated to be tested here.
  // The tests for the DatabaseSuitability TOPP tool have to suffice.
}
END_SECTION

START_SECTION(getResults())
{
  NOT_TESTABLE;
}
END_SECTION

DBSuitability_friend private_suit;

START_SECTION(std::vector<FASTAFile::FASTAEntry> getSubsampledFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_data, double subsampling_rate) const)
{
  vector<FASTAFile::FASTAEntry> fasta;
  FASTAFile::FASTAEntry entry;
  entry.sequence = "AAAAAAA";// 7
  fasta.push_back(entry);
  entry.sequence = "PP";// 2
  fasta.push_back(entry);
  entry.sequence = "EEE";// 3
  fasta.push_back(entry);
  entry.sequence = "I";// 1
  fasta.push_back(entry);
  entry.sequence = "KKKKKK";// 6
  fasta.push_back(entry);
  entry.sequence = "LLLLL";// 5
  fasta.push_back(entry);
  entry.sequence = "QQQQ";//4
  fasta.push_back(entry);
  entry.sequence = "YYY";// 3
  fasta.push_back(entry);
  entry.sequence = "GGGG";// 4
  fasta.push_back(entry);
  // 35 AS in fasta

  vector<FASTAFile::FASTAEntry> subsampled_fasta = private_suit.getSubsampledFasta(fasta, 0.3); // 35 * 0.3 = 10.5 --> at least 11 AS should be written (& at max. 17)

  TEST_EQUAL((countAS(subsampled_fasta) >= 11 && countAS(subsampled_fasta) < 17), 1);
  TEST_EXCEPTION(Exception::IllegalArgument, private_suit.getSubsampledFasta(fasta, 2));
  TEST_EXCEPTION(Exception::IllegalArgument, private_suit.getSubsampledFasta(fasta, -1));
}
END_SECTION

START_SECTION(void appendDecoys_(std::vector<FASTAFile::FASTAEntry>& fasta) const)
{
  vector<FASTAFile::FASTAEntry> fasta;
  FASTAFile::FASTAEntry entry;
  entry.sequence = "LIEQKPABIM";
  entry.identifier = "PROTEIN";
  fasta.push_back(entry);

  private_suit.appendDecoys(fasta);

  TEST_STRING_EQUAL(fasta[1].sequence, "LIBAPKQEIM");
  TEST_STRING_EQUAL(fasta[1].identifier, "DECOY_PROTEIN");
}
END_SECTION

START_SECTION(double calculateCorrectionFactor_(const DBSuitability::SuitabilityData& data, const DBSuitability::SuitabilityData& data_sampled, double sampling_rate) const)
{
  DBSuitability::SuitabilityData full_data;
  DBSuitability::SuitabilityData subsampled_data;

  full_data.num_top_db = 100;
  subsampled_data.num_top_db = 50;
  // delta 50

  full_data.num_top_novo = 10;
  subsampled_data.num_top_novo = 30;
  // delta 20

  double factor = private_suit.calculateCorrectionFactor(full_data, subsampled_data, 0.6);
  // rate 0.6 --> db_slope = -50 / -0.4 = 125, novo_slope = 20 / -0.4 = -50
  // factor = - (125) / (-50) = 2.5

  TEST_EQUAL(factor, 2.5);
  TEST_EXCEPTION(Exception::Precondition, private_suit.calculateCorrectionFactor(full_data, subsampled_data, 2));
  TEST_EXCEPTION(Exception::Precondition, private_suit.calculateCorrectionFactor(full_data, subsampled_data, -1));
}
END_SECTION

START_SECTION(UInt numberOfUniqueProteins_(const std::vector<PeptideIdentification>& peps, UInt number_of_hits = 1) const)
{
  PeptideEvidence ev1("PROTEIN_1", 0, 0, '[', ']');
  PeptideEvidence ev2("PROTEIN_2", 0, 0, '[', ']');
  PeptideEvidence ev3("PROTEIN_3", 0, 0, '[', ']');
  PeptideEvidence ev4("PROTEIN_4", 0, 0, '[', ']');
  PeptideEvidence ev5("DECOY_PROTEIN", 0, 0, '[', ']');

  PeptideHit hit1;
  hit1.setPeptideEvidences({ev1, ev1, ev2});
  hit1.setMetaValue("target_decoy", "target");
  PeptideHit hit2;
  hit2.setPeptideEvidences({ev4, ev3, ev5});
  hit2.setMetaValue("target_decoy", "target+decoy");
  PeptideHit hit3;
  hit3.setPeptideEvidences({ev3, ev2, ev3});
  hit3.setMetaValue("target_decoy", "target");
  PeptideHit hit4;
  hit4.setPeptideEvidences({ev5});
  hit4.setMetaValue("target_decoy", "decoy");
  PeptideHit empty_hit;

  PeptideIdentification id1;
  id1.setHits({hit1, hit2});
  PeptideIdentification id2;
  id2.setHits({hit3});
  PeptideIdentification id3;
  id3.setHits({hit4});
  PeptideIdentification empty_id;
  PeptideIdentification id_hit_without_info;
  id_hit_without_info.setHits({empty_hit});

  vector<PeptideIdentification> ids({id1, id2, empty_id, id3});

  TEST_EQUAL(private_suit.numberOfUniqueProteins(ids), 3);
  TEST_EQUAL(private_suit.numberOfUniqueProteins(ids, 5), 4);
  TEST_EXCEPTION(Exception::MissingInformation, private_suit.numberOfUniqueProteins({id_hit_without_info}));
}
END_SECTION

START_SECTION(Size getIndexWithMedianNovoHits_(const std::vector<SuitabilityData>& data) const)
{
  DBSuitability::SuitabilityData d1;
  d1.num_top_novo = 10;
  DBSuitability::SuitabilityData d2;
  d2.num_top_novo = 20;
  DBSuitability::SuitabilityData d3;
  d3.num_top_novo = 15;
  DBSuitability::SuitabilityData d4;
  d4.num_top_novo = 40;

  TEST_EQUAL(private_suit.getIndexWithMedianNovoHits({d1, d2, d3}), 2);
  TEST_EQUAL(private_suit.getIndexWithMedianNovoHits({d1, d2, d3, d4}), 1);
  TEST_EXCEPTION(Exception::IllegalArgument, private_suit.getIndexWithMedianNovoHits({}));
}
END_SECTION

START_SECTION(double getScoreMatchingFDR_(const std::vector<PeptideIdentification>& pep_ids, double FDR, String score_name, bool higher_score_better) const)
{
  PeptideHit hit1;
  hit1.setScore(0.01);
  hit1.setMetaValue("some_score", 120);
  PeptideHit hit2;
  hit2.setScore(0.04);
  hit2.setMetaValue("some_score", 80);
  PeptideHit hit3;
  hit3.setScore(0.5);
  hit3.setMetaValue("some_score", 5);
  PeptideHit hit4;
  hit4.setScore(0.05);
  hit4.setMetaValue("some_score", 75);

  PeptideIdentification id1;
  id1.setScoreType("q-value");
  id1.setHits({hit1});
  PeptideIdentification id2;
  id2.setScoreType("q-value");
  id2.setHits({hit2});
  PeptideIdentification id3;
  id3.setScoreType("q-value");
  id3.setHits({hit3});
  PeptideIdentification id4;
  id4.setScoreType("q-value");
  id4.setHits({hit4});

  TEST_EQUAL(private_suit.getScoreMatchingFDR({id1, id2, id3, id4}, 0.05, "some_score", true), 75);
  TEST_EQUAL(private_suit.getScoreMatchingFDR({id1, id2, id3, id4}, 0.05, "some", false), 120);
  TEST_EXCEPTION(Exception::IllegalArgument, private_suit.getScoreMatchingFDR({id1}, 0.05, "e-value", false));
  id1.setScoreType("e-value");
  TEST_EXCEPTION(Exception::Precondition, private_suit.getScoreMatchingFDR({id1}, 0.05, "some_score", false));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST