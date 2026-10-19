// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/EntrapmentEvaluation.h>

#include <limits>
#include <unordered_set>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  PeptideIdentification makeId(const std::string& seq, double q, const std::vector<std::string>& accessions,
                               bool decoy = false, const std::string& cls = "")
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setScore(q);
    for (const auto& acc : accessions)
    {
      PeptideEvidence pe;
      pe.setProteinAccession(acc);
      hit.addPeptideEvidence(pe);
    }
    hit.setTargetDecoyType(decoy ? PeptideHit::TargetDecoyType::DECOY : PeptideHit::TargetDecoyType::TARGET);
    if (!cls.empty()) hit.setMetaValue("mod_class", cls);
    PeptideIdentification id;
    id.setScoreType("q-value");
    id.setHigherScoreBetter(false);
    id.insertHit(hit);
    return id;
  }
}

START_TEST(EntrapmentEvaluation, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static Kind classify(const PeptideHit& hit, const std::string& entrapment_prefix)))
{
  using K = EntrapmentEvaluation::Kind;
  TEST_EQUAL(EntrapmentEvaluation::classify(makeId("PEPTIDEK", 0, {"P1"}).getHits()[0], "ENTRAP_") == K::TARGET, true)
  TEST_EQUAL(EntrapmentEvaluation::classify(makeId("PEPTIDEK", 0, {"ENTRAP_P1"}).getHits()[0], "ENTRAP_") == K::ENTRAPMENT, true)
  // shared with a target protein: target
  TEST_EQUAL(EntrapmentEvaluation::classify(makeId("PEPTIDEK", 0, {"ENTRAP_P1", "P2"}).getHits()[0], "ENTRAP_") == K::TARGET, true)
  TEST_EQUAL(EntrapmentEvaluation::classify(makeId("PEPTIDEK", 0, {"DECOY_P1"}, true).getHits()[0], "ENTRAP_") == K::DECOY, true)
}
END_SECTION

START_SECTION((static double qValue(const PeptideIdentification& pep_id, const PeptideHit& hit)))
{
  PeptideIdentification id = makeId("PEPTIDEK", 0.02, {"P1"});
  TEST_REAL_SIMILAR(EntrapmentEvaluation::qValue(id, id.getHits()[0]), 0.02)
  id.setScoreType("hyperscore");
  TEST_EXCEPTION(Exception::MissingInformation, EntrapmentEvaluation::qValue(id, id.getHits()[0]))
  PeptideHit h = id.getHits()[0];
  h.setMetaValue("q-value", 0.005);
  TEST_REAL_SIMILAR(EntrapmentEvaluation::qValue(id, h), 0.005)
}
END_SECTION

START_SECTION((static std::map<std::string, Result> evaluate(const PeptideIdentificationList& pep_ids, double q_threshold, const std::string& entrapment_prefix, double te_ratio, const std::string& class_meta_value = "", const std::unordered_set<std::string>* target_peptides = nullptr)))
{
  PeptideIdentificationList ids;
  // 8 target PSMs (7 distinct peptides), 2 entrapment PSMs (2 peptides), 1 decoy, 1 above threshold
  for (const auto& s : {"AAAK", "CCCK", "DDDK", "EEEK", "FFFK", "GGGK", "HHHK"}) ids.push_back(makeId(s, 0.001, {"P1"}, false, "unmod"));
  ids.push_back(makeId("AAAK", 0.002, {"P1"}, false, "unmod"));
  ids.push_back(makeId("KLLLK", 0.005, {"ENTRAP_P9"}, false, "mod"));
  ids.push_back(makeId("KMMMK", 0.009, {"ENTRAP_P9"}, false, "mod"));
  ids.push_back(makeId("NNNK", 0.001, {"DECOY_P1"}, true));
  ids.push_back(makeId("PPPK", 0.5, {"ENTRAP_P2"}));

  auto res = EntrapmentEvaluation::evaluate(ids, 0.01, "ENTRAP_", 1.0, "mod_class");
  const auto& all = res["all"];
  TEST_EQUAL(all.accepted_psms, 10)
  TEST_EQUAL(all.entrapment_psms, 2)
  TEST_REAL_SIMILAR(all.psmFDP(), 2.0 / 10.0 * 2.0)
  TEST_EQUAL(all.accepted_peptides, 9)
  TEST_EQUAL(all.entrapment_peptides, 2)
  TEST_REAL_SIMILAR(all.peptideFDP(), 2.0 / 9.0 * 2.0)
  TEST_EQUAL(res["mod"].accepted_psms, 2)
  TEST_REAL_SIMILAR(res["mod"].psmFDP(), 2.0)
  TEST_EQUAL(res["unmod"].entrapment_psms, 0)
  TEST_EQUAL(res.count("NA"), 0)

  // a class named "all" does not add to the overall result
  PeptideIdentificationList ids_all = ids;
  for (auto& id : ids_all) if (id.getHits()[0].getMetaValue("mod_class") == "mod") id.getHits()[0].setMetaValue("mod_class", "all");
  auto res_all = EntrapmentEvaluation::evaluate(ids_all, 0.01, "ENTRAP_", 1.0, "mod_class");
  TEST_EQUAL(res_all["all"].accepted_psms, 10)
  TEST_EQUAL(res_all["class:all"].accepted_psms, 2)

  // T/E = 3: entrapment hits weigh less
  res = EntrapmentEvaluation::evaluate(ids, 0.01, "ENTRAP_", 3.0);
  TEST_REAL_SIMILAR(res["all"].psmFDP(), 2.0 / 10.0 * 4.0)
  TEST_EQUAL(res.size(), 1)

  // nothing accepted
  res = EntrapmentEvaluation::evaluate(ids, 0.0, "ENTRAP_", 1.0);
  TEST_EQUAL(res["all"].accepted_psms, 0)
  TEST_REAL_SIMILAR(res["all"].psmFDP(), 0.0)

  // a NaN or negative q-value is never accepted
  ids.push_back(makeId("QQQK", std::numeric_limits<double>::quiet_NaN(), {"P1"}));
  ids.push_back(makeId("RRRK", -0.1, {"P1"}));
  res = EntrapmentEvaluation::evaluate(ids, 0.01, "ENTRAP_", 1.0);
  TEST_EQUAL(res["all"].accepted_psms, 10)

  // with the target peptides of the database: KILLK (entrapment) is KLLLK with I = L, which a target protein holds,
  // so it counts as target; KMMMK stays entrapment
  ids.push_back(makeId("KILLK", 0.003, {"ENTRAP_P8"}));
  const std::unordered_set<std::string> target_peptides = {"KLLLK", "AAAK"};
  res = EntrapmentEvaluation::evaluate(ids, 0.01, "ENTRAP_", 1.0, "", &target_peptides);
  TEST_EQUAL(res["all"].accepted_psms, 11)
  TEST_EQUAL(res["all"].entrapment_psms, 1)
  TEST_EQUAL(res["all"].entrapment_peptides, 1)
  res = EntrapmentEvaluation::evaluate(ids, 0.01, "ENTRAP_", 1.0);
  TEST_EQUAL(res["all"].entrapment_psms, 3)
}
END_SECTION

START_SECTION((static double estimateTERatio(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& entrapment_prefix, const std::string& decoy_string, bool decoy_is_prefix, const std::string& enzyme = "Trypsin", Size missed_cleavages = 2, Size min_length = 7, Size max_length = 40)))
{
  std::vector<FASTAFile::FASTAEntry> db;
  db.emplace_back("P1", "", "AAAAAAAKCCCCCCCKDDDDDDDK");
  db.emplace_back("ENTRAP_P1", "", "AAAAAAAKEEEEEEEK"); // AAAAAAAK is shared with the target: not entrapment
  db.emplace_back("DECOY_P1", "", "GGGGGGGKHHHHHHHKIIIIIIIK");
  // targets (0 missed cleavages): 3; entrapment-only: 1
  TEST_REAL_SIMILAR(EntrapmentEvaluation::estimateTERatio(db, "ENTRAP_", "DECOY_", true, "Trypsin", 0, 7, 40), 3.0)
  // suffix decoys, also numbered as DecoyDatabase writes them with shuffle_decoy_ratio > 1
  db[2].identifier = "P1_DECOY";
  db.emplace_back("P1_DECOY_2", "", "MMMMMMMRLLLLLLLK");
  TEST_REAL_SIMILAR(EntrapmentEvaluation::estimateTERatio(db, "ENTRAP_", "_DECOY", false, "Trypsin", 0, 7, 40), 3.0)
  // a suffix followed by something other than "_<n>" is not a decoy: its 2 peptides count as targets
  db.back().identifier = "P1_DECOY_X";
  TEST_REAL_SIMILAR(EntrapmentEvaluation::estimateTERatio(db, "ENTRAP_", "_DECOY", false, "Trypsin", 0, 7, 40), 5.0)
  db.pop_back();
  db[2].identifier = "DECOY_P1";
  // no entrapment: 1
  db.erase(db.begin() + 1);
  TEST_REAL_SIMILAR(EntrapmentEvaluation::estimateTERatio(db, "ENTRAP_", "DECOY_", true, "Trypsin", 0, 7, 40), 1.0)
  // I = L: the entrapment peptide CCCCLCCK equals the target peptide CCCCICCK, EEEEEEEK is entrapment-only
  db[0].sequence = "AAAAAAAKCCCCICCKDDDDDDDK";
  db.emplace_back("ENTRAP_P2", "", "CCCCLCCKEEEEEEEK");
  TEST_REAL_SIMILAR(EntrapmentEvaluation::estimateTERatio(db, "ENTRAP_", "DECOY_", true, "Trypsin", 0, 7, 40), 3.0)
}
END_SECTION

START_SECTION((static std::unordered_set<std::string> targetPeptides(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& entrapment_prefix, const std::string& decoy_string, bool decoy_is_prefix, const std::string& enzyme = "Trypsin", Size missed_cleavages = 2, Size min_length = 7, Size max_length = 40)))
{
  std::vector<FASTAFile::FASTAEntry> db;
  db.emplace_back("P1", "", "AAAAAAAKCCCCICCK");
  db.emplace_back("ENTRAP_P1", "", "EEEEEEEK");
  db.emplace_back("DECOY_P1", "", "GGGGGGGK");
  const auto peptides = EntrapmentEvaluation::targetPeptides(db, "ENTRAP_", "DECOY_", true, "Trypsin", 0, 7, 40);
  TEST_EQUAL(peptides.size(), 2)
  TEST_EQUAL(peptides.contains("AAAAAAAK"), true)
  TEST_EQUAL(peptides.contains("CCCCLCCK"), true) // I replaced by L
  TEST_EQUAL(peptides.contains("EEEEEEEK"), false)
  TEST_EQUAL(peptides.contains("GGGGGGGK"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
