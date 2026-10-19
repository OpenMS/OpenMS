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

START_SECTION((static std::map<std::string, Result> evaluate(const PeptideIdentificationList& pep_ids, double q_threshold, const std::string& entrapment_prefix, double te_ratio, const std::string& class_meta_value = "")))
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
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
