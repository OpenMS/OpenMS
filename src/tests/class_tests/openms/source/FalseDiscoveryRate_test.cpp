// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FalseDiscoveryRate, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FalseDiscoveryRate* ptr = nullptr;
FalseDiscoveryRate* nullPointer = nullptr;
START_SECTION(FalseDiscoveryRate())
{
  ptr = new FalseDiscoveryRate();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FalseDiscoveryRate())
{
  delete ptr;
}
END_SECTION

START_SECTION((void apply(PeptideIdentificationList &fwd_ids, PeptideIdentificationList &rev_ids)))
{
  ptr = new FalseDiscoveryRate();
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids;
  PeptideIdentificationList fwd_pep_ids, rev_pep_ids;
  std::string document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);
  ptr->apply(fwd_pep_ids, rev_pep_ids);
  TOLERANCE_ABSOLUTE(0.0001)
  for (PeptideIdentificationList::const_iterator it = fwd_pep_ids.begin(); it != fwd_pep_ids.end(); ++it)
  {
    if (!it->getHits().empty())
    {
      PeptideHit hit(*it->getHits().begin());
      double fdr(hit.getScore());
      double orig_score((double)hit.getMetaValue("XTandem_score"));
      
      if (orig_score >= 39.4)
      {
        TEST_REAL_SIMILAR(fdr, 0)
      }
      if (orig_score <= 37.9 + 0.0001 && orig_score >= 37.9 - 0.0001)
      {
        TEST_REAL_SIMILAR(fdr, 0.08)
      }
    }
  }
}
END_SECTION

START_SECTION((void apply(std::vector<ProteinIdentification> &fwd_ids, std::vector<ProteinIdentification> &rev_ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids;
  PeptideIdentificationList fwd_pep_ids, rev_pep_ids;
  std::string document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids_withProtScores.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids_withProtScores.idXML"), rev_prot_ids, rev_pep_ids, document_id);
  ptr->apply(fwd_prot_ids, rev_prot_ids);
  TOLERANCE_ABSOLUTE(0.001)
	
  for (vector<ProteinIdentification>::const_iterator prot_it = fwd_prot_ids.begin(); prot_it != fwd_prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));
        
        // it gets here, but neither of the conditions below are ever satisfied
        if (orig_score < -1.8)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
        if (orig_score == 0.0)
        {
          TEST_REAL_SIMILAR(fdr, 0.897384)
        }
      }
    }
  }
}
END_SECTION

START_SECTION((void apply(PeptideIdentificationList &id)))
{
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_OMSSA.idXML"), prot_ids, pep_ids);

  ptr->apply(pep_ids);
  TOLERANCE_ABSOLUTE(0.001)
	
  for (Size z = 1; z <= 4; ++z)
  {
    for (PeptideIdentificationList::const_iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
  	{
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
      	double fdr(pit->getScore());
      	double orig_score((double)pit->getMetaValue("OMSSA_score"));

      	if (orig_score <= 10e-4)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
      	if (orig_score >= 1000 && pit->getCharge() != 1)
        {
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }

    // target hit
    PeptideIdentification pep_id = pep_ids[0];
    PeptideHit pit = pep_id.getHits()[0];
    double fdr(pit.getScore());
    TEST_REAL_SIMILAR(fdr, 0.0730478589420655);

    // target+decoy hit considered as target
    pep_id = pep_ids[5];
    pit = pep_id.getHits()[0];
    fdr = pit.getScore();
    TEST_REAL_SIMILAR(fdr, 0.409926470588235);

    // decoy hit removed
    pep_id = pep_ids[9];
    TEST_EQUAL(pep_id.getHits().size(), 0)
  }
}
END_SECTION

START_SECTION((void apply(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  PeptideIdentificationList fwd_pep_ids, rev_pep_ids, pep_ids;
  std::string document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

  for (vector<ProteinIdentification>::const_iterator it = fwd_prot_ids.begin(); it != fwd_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "target");
    }
  }
  
  for (vector<ProteinIdentification>::const_iterator it = rev_prot_ids.begin(); it != rev_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "decoy");
    }
  }

  ptr->apply(prot_ids);

  TOLERANCE_ABSOLUTE(0.001)

  for (vector<ProteinIdentification>::const_iterator prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));

        if (orig_score < -1.8)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
        if (orig_score == -1.7)
        {
          TEST_REAL_SIMILAR(fdr, 0.0617284)
        }
        if (orig_score > -1.2)
        {
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }
  }
}
END_SECTION

START_SECTION((void applyPicked(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_picked_in.idXML"), prot_ids, pep_ids);

  ptr->applyPickedProteinFDR(prot_ids[0],"decoy_");

  TOLERANCE_ABSOLUTE(0.001)
  const auto& hits = prot_ids[0].getHits();
  TEST_REAL_SIMILAR(hits[0].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[1].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[2].getScore(),0.25)
  TEST_REAL_SIMILAR(hits[3].getScore(),0.4)
  TEST_REAL_SIMILAR(hits[4].getScore(),0.4)
  TEST_REAL_SIMILAR(hits[5].getScore(),0.5)
}
END_SECTION

START_SECTION((void applyPickedProteinFDR(ProteinIdentification& id, std::string decoy_string, bool prefix, bool groups_too)))
{
  // Hand-computed oracles for picked protein-GROUP q-values (defaults: q-values, conservative FDR
  // (D+1)/(T+1) over the voting groups at least as good as a score, then a cumulative minimum).
  // A group votes target if a member's target won its target/decoy pair, else decoy if a member's
  // decoy won it; groups whose members all lost do not vote. Decoys carry the prefix DECOY_.
  auto picked_run = [](bool higher_better, const vector<pair<string, double>>& hits,
                       const vector<pair<double, vector<string>>>& groups)
  {
    ProteinIdentification run;
    run.setScoreType("score");
    run.setHigherScoreBetter(higher_better);
    for (const auto& [acc, score] : hits)
    {
      ProteinHit hit(score, 1, acc, "");
      hit.setMetaValue("target_decoy", acc.rfind("DECOY_", 0) == 0 ? "decoy" : "target");
      run.getHits().push_back(hit);
    }
    for (const auto& [score, accs] : groups)
    {
      ProteinIdentification::ProteinGroup grp;
      grp.probability = score;
      grp.accessions = accs;
      run.getIndistinguishableProteins().push_back(grp);
    }
    return run;
  };
  auto group_qvalues = [](const ProteinIdentification& run)
  {
    vector<double> qs;
    for (const auto& grp : run.getIndistinguishableProteins()) qs.push_back(grp.probability);
    return qs;
  };
  FalseDiscoveryRate fdr;
  TOLERANCE_ABSOLUTE(1e-6)

  // (1) lower score = better. Pairs: A wins (0.01 < 0.02), DECOY_B wins (0.04 < 0.05), C has no decoy.
  // Votes: {A} 0.01 T, {C} 0.02 T, {DECOY_B} 0.04 D; {DECOY_A} and {B} lost and do not vote.
  // FDR: 0.01 -> 1/2, 0.02 -> 1/3, 0.04 -> 2/3; q: 0.01 -> 1/3, 0.02 -> 1/3, 0.04 -> 2/3.
  // {B} (0.05) is worse than every calibrated score: q = 2/3 (all calibrated groups accepted).
  // Before the fix, its lookup dereferenced the end of the score map (ASan stack-buffer-overflow).
  {
    ProteinIdentification run = picked_run(false,
      {{"A", 0.01}, {"DECOY_A", 0.02}, {"B", 0.05}, {"DECOY_B", 0.04}, {"C", 0.02}},
      {{0.01, {"A"}}, {0.02, {"DECOY_A"}}, {0.05, {"B"}}, {0.04, {"DECOY_B"}}, {0.02, {"C"}}});
    fdr.applyPickedProteinFDR(run, "DECOY_", true);
    vector<double> qs = group_qvalues(run);
    TEST_REAL_SIMILAR(qs[0], 1.0 / 3.0)
    TEST_REAL_SIMILAR(qs[1], 1.0 / 3.0)
    TEST_REAL_SIMILAR(qs[2], 2.0 / 3.0)
    TEST_REAL_SIMILAR(qs[3], 2.0 / 3.0)
    TEST_REAL_SIMILAR(qs[4], 1.0 / 3.0)
  }

  // (2) higher score = better. DECOY_A (0.95, in no group) beats A (0.9). Only {C} votes: 0.5 T -> q 1/2.
  // {A} (0.9) is better than every calibrated score: q = 1/2, the best calibrated q-value.
  // Before the fix, its lookup dereferenced the end of the score map.
  {
    ProteinIdentification run = picked_run(true,
      {{"A", 0.9}, {"DECOY_A", 0.95}, {"C", 0.5}},
      {{0.9, {"A"}}, {0.5, {"C"}}});
    fdr.applyPickedProteinFDR(run, "DECOY_", true);
    vector<double> qs = group_qvalues(run);
    TEST_REAL_SIMILAR(qs[0], 0.5)
    TEST_REAL_SIMILAR(qs[1], 0.5)
  }

  // (3) No group votes (the only group lost its pair): nothing to calibrate group q-values on.
  {
    ProteinIdentification run = picked_run(true, {{"A", 0.9}, {"DECOY_A", 0.95}}, {{0.9, {"A"}}});
    TEST_EXCEPTION(Exception::MissingInformation, fdr.applyPickedProteinFDR(run, "DECOY_", true))
  }

  // (4) Score direction. Lower score = better. Pairs: P1, P4 have no decoy; P2 wins (0.02 < 0.025);
  // DECOY_P3 wins (0.03 < 0.06). Votes: 0.01 T, 0.02 T, 0.03 D, 0.05 T.
  // FDR: 0.01 -> 1/2, 0.02 -> 1/3, 0.03 -> 2/3, 0.05 -> 1/2; q: 1/3, 1/3, 1/2, 1/2.
  // {DECOY_P2} (0.025) lies between 0.02 and 0.03: its closest equal-or-better calibrated score is
  // 0.02, so q = 1/3 (as for a higher-better score). Before the fix it got the WORSE neighbour's 1/2.
  // The same fixture with negated scores and higher = better must give the same q-values.
  {
    for (bool higher_better : {false, true})
    {
      const double s = higher_better ? -1.0 : 1.0;
      ProteinIdentification run = picked_run(higher_better,
        {{"P1", s * 0.01}, {"P2", s * 0.02}, {"DECOY_P2", s * 0.025}, {"DECOY_P3", s * 0.03}, {"P3", s * 0.06}, {"P4", s * 0.05}},
        {{s * 0.01, {"P1"}}, {s * 0.02, {"P2"}}, {s * 0.025, {"DECOY_P2"}}, {s * 0.03, {"DECOY_P3"}}, {s * 0.06, {"P3"}}, {s * 0.05, {"P4"}}});
      fdr.applyPickedProteinFDR(run, "DECOY_", true);
      vector<double> qs = group_qvalues(run);
      TEST_REAL_SIMILAR(qs[0], 1.0 / 3.0)
      TEST_REAL_SIMILAR(qs[1], 1.0 / 3.0)
      TEST_REAL_SIMILAR(qs[2], 1.0 / 3.0)
      TEST_REAL_SIMILAR(qs[3], 0.5)
      TEST_REAL_SIMILAR(qs[4], 0.5)
      TEST_REAL_SIMILAR(qs[5], 0.5)
    }
  }

  // (5) Mixed target/decoy group, one vote, independent of accession order. Higher score = better.
  // Pairs: Z, W have no decoy; Y ties with DECOY_Y (0.8) -> the target wins; DECOY_X (0.8) beats X (0.5).
  // {DECOY_X, Y} holds a winning decoy AND a winning target: it votes once, as target.
  // Votes: 0.9 T, 0.8 T, 0.7 T. FDR: 1/2, 1/3, 1/4; q: 1/4 for all. {DECOY_Y} (0.8) and {X} (0.5)
  // do not vote and get q = 1/4. Before the fix, the decoy-first order added a second (decoy) vote
  // at 0.8, giving q = 1/2 for every group, while the target-first order gave 1/4.
  {
    for (const vector<string>& mixed : {vector<string>{"DECOY_X", "Y"}, vector<string>{"Y", "DECOY_X"}})
    {
      ProteinIdentification run = picked_run(true,
        {{"Z", 0.9}, {"Y", 0.8}, {"DECOY_Y", 0.8}, {"X", 0.5}, {"DECOY_X", 0.8}, {"W", 0.7}},
        {{0.9, {"Z"}}, {0.8, mixed}, {0.8, {"DECOY_Y"}}, {0.5, {"X"}}, {0.7, {"W"}}});
      fdr.applyPickedProteinFDR(run, "DECOY_", true);
      for (double q : group_qvalues(run))
      {
        TEST_REAL_SIMILAR(q, 0.25)
      }
    }
  }

  // (6) Raw FDR (no_qvalues). Fixture (2): {A} (0.9) is better than every calibrated score, so no
  // voting group is accepted at its score: FDR (0+1)/(0+1) = 1 (a q-value takes the minimum, 1/2).
  // {C}: 1/2. Fixture (1), lower = better: conservative FDR 0.01 -> 1/2, 0.02 -> 1/3, 0.04 -> 2/3;
  // with conservative = false (D+1)/(D+T+1), 0.04 -> 2/4. {B} (0.05) gets the value of 0.04.
  {
    FalseDiscoveryRate raw;
    Param p = raw.getParameters();
    p.setValue("no_qvalues", "true");
    raw.setParameters(p);
    for (bool higher_better : {true, false})
    {
      const double s = higher_better ? 1.0 : -1.0;
      ProteinIdentification run = picked_run(higher_better,
        {{"A", s * 0.9}, {"DECOY_A", s * 0.95}, {"C", s * 0.5}},
        {{s * 0.9, {"A"}}, {s * 0.5, {"C"}}});
      raw.applyPickedProteinFDR(run, "DECOY_", true);
      vector<double> fdrs = group_qvalues(run);
      TEST_REAL_SIMILAR(fdrs[0], 1.0)
      TEST_REAL_SIMILAR(fdrs[1], 0.5)
    }
    for (bool conservative : {true, false})
    {
      p.setValue("conservative", conservative ? "true" : "false");
      raw.setParameters(p);
      ProteinIdentification run = picked_run(false,
        {{"A", 0.01}, {"DECOY_A", 0.02}, {"B", 0.05}, {"DECOY_B", 0.04}, {"C", 0.02}},
        {{0.01, {"A"}}, {0.02, {"DECOY_A"}}, {0.05, {"B"}}, {0.04, {"DECOY_B"}}, {0.02, {"C"}}});
      raw.applyPickedProteinFDR(run, "DECOY_", true);
      vector<double> fdrs = group_qvalues(run);
      const double last = conservative ? 2.0 / 3.0 : 0.5;
      TEST_REAL_SIMILAR(fdrs[0], 0.5)
      TEST_REAL_SIMILAR(fdrs[1], 1.0 / 3.0)
      TEST_REAL_SIMILAR(fdrs[2], last)
      TEST_REAL_SIMILAR(fdrs[3], last)
      TEST_REAL_SIMILAR(fdrs[4], 1.0 / 3.0)
    }
  }
}
END_SECTION


START_SECTION((void apply(std::vector<ProteinIdentification>& ids)))
{
  vector<ProteinIdentification> fwd_prot_ids, rev_prot_ids, prot_ids;
  PeptideIdentificationList fwd_pep_ids, rev_pep_ids, pep_ids;
  std::string document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_fwd_ids.idXML"), fwd_prot_ids, fwd_pep_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("XTandem_rev_ids.idXML"), rev_prot_ids, rev_pep_ids, document_id);

  for (vector<ProteinIdentification>::const_iterator it = fwd_prot_ids.begin(); it != fwd_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "target");
    }
  }

  for (vector<ProteinIdentification>::const_iterator it = rev_prot_ids.begin(); it != rev_prot_ids.end(); ++it)
  {
    prot_ids.push_back(*it);
    for (vector<ProteinHit>::iterator hit_it = prot_ids.back().getHits().begin(); hit_it != prot_ids.back().getHits().end(); ++hit_it)
    {
      hit_it->setMetaValue("target_decoy", "decoy");
    }
  }

  ptr->apply(prot_ids);

  TOLERANCE_ABSOLUTE(0.001)

  for (vector<ProteinIdentification>::const_iterator prot_it = prot_ids.begin(); prot_it != prot_ids.end(); ++prot_it)
  {
    if (!prot_it->getHits().empty())
    {
      for (vector<ProteinHit>::const_iterator it = prot_it->getHits().begin(); it != prot_it->getHits().end(); ++it)
      {
        ProteinHit hit(*it);
        double fdr(hit.getScore());
        double orig_score((double)hit.getMetaValue("XTandem_score"));

        if (orig_score < -1.8)
        {
          TEST_REAL_SIMILAR(fdr, 0)
        }
        if (orig_score == -1.7)
        {
          TEST_REAL_SIMILAR(fdr, 0.0617284)
        }
        if (orig_score > -1.2)
        {
          TEST_EQUAL(fdr > 0.1, true)
        }
      }
    }
  }
}
END_SECTION

START_SECTION((void applyBasicPeptideLevel(PeptideIdentificationList & ids)))
{
  // Regression test for issue #9767 finding #3:
  // addToPeptideScoreMap_() must keep the BEST score per unmodified peptide
  // sequence, regardless of the order in which duplicate sequences are seen.
  // Before the fix the try_emplace "inserted" flag was misused, so the map kept
  // the FIRST-seen score and the peptide-level q-value assigned to a duplicated
  // sequence depended on input order. We assert order-independence, which needs
  // no hand-computed FDR values and fails before / passes after the fix.
  auto build_ids = [](bool dup_worse_first) -> PeptideIdentificationList
  {
    PeptideIdentificationList ids;
    auto add = [&ids](const std::string& seq, double score, const std::string& td)
    {
      PeptideHit h;
      h.setSequence(AASequence::fromString(seq));
      h.setScore(score);
      h.setCharge(2);
      h.setMetaValue("target_decoy", td);
      PeptideIdentification pid;
      pid.setScoreType("score");
      pid.setHigherScoreBetter(true);
      pid.setHits({h});
      ids.push_back(pid);
    };
    // targets spanning the score range
    add("AAAAAAAK", 0.90, "target");
    add("CCCCCCCK", 0.80, "target");
    add("DDDDDDDK", 0.70, "target");
    // decoys sitting between the duplicate's two candidate scores (0.10 and 0.95)
    add("FFFFFFFK", 0.50, "decoy");
    add("HHHHHHHK", 0.30, "decoy");
    // the SAME sequence twice (both targets), in the requested order
    if (dup_worse_first)
    {
      add("EEEEEEEK", 0.10, "target");
      add("EEEEEEEK", 0.95, "target");
    }
    else
    {
      add("EEEEEEEK", 0.95, "target");
      add("EEEEEEEK", 0.10, "target");
    }
    return ids;
  };

  // peptide q-value assigned to a surviving EEEEEEEK identification (-1 if absent)
  auto dup_qvalue = [](const PeptideIdentificationList& ids) -> double
  {
    for (const auto& id : ids)
    {
      if (id.getHits().empty()) continue;
      if (id.getHits()[0].getSequence().toUnmodifiedString() == "EEEEEEEK")
      {
        return id.getHits()[0].getScore();
      }
    }
    return -1.0;
  };

  FalseDiscoveryRate fdr;

  PeptideIdentificationList ids_worse_first = build_ids(true);
  PeptideIdentificationList ids_better_first = build_ids(false);

  fdr.applyBasicPeptideLevel(ids_worse_first);
  fdr.applyBasicPeptideLevel(ids_better_first);

  const double q_worse_first = dup_qvalue(ids_worse_first);
  const double q_better_first = dup_qvalue(ids_better_first);

  // both orders must resolve to a real (found) q-value ...
  TEST_NOT_EQUAL(q_worse_first, -1.0)
  TEST_NOT_EQUAL(q_better_first, -1.0)
  // ... and, because the BEST score (0.95) is the representative regardless of
  // input order, the two runs must agree (before the fix they differed).
  TEST_REAL_SIMILAR(q_worse_first, q_better_first)

  // Label consistency: when the better duplicate is a TARGET but an earlier duplicate
  // of the same sequence was a DECOY, the stored representative must become
  // (best score, target) and not keep the earlier decoy label. Otherwise the result is
  // again order-dependent (a high-scoring "decoy" would distort the FDR). We assert
  // order-independence for a mixed-label duplicate.
  auto build_mixed = [](bool decoy_first, double decoy_score, double target_score) -> PeptideIdentificationList
  {
    PeptideIdentificationList ids;
    auto add = [&ids](const std::string& seq, double score, const std::string& td)
    {
      PeptideHit h;
      h.setSequence(AASequence::fromString(seq));
      h.setScore(score);
      h.setCharge(2);
      h.setMetaValue("target_decoy", td);
      PeptideIdentification pid;
      pid.setScoreType("score");
      pid.setHigherScoreBetter(true);
      pid.setHits({h});
      ids.push_back(pid);
    };
    add("AAAAAAAK", 0.90, "target");
    add("CCCCCCCK", 0.80, "target");
    add("DDDDDDDK", 0.70, "target");
    add("FFFFFFFK", 0.50, "decoy");
    add("HHHHHHHK", 0.30, "decoy");
    // same sequence twice with DIFFERENT labels, in the requested order
    if (decoy_first)
    {
      add("EEEEEEEK", decoy_score, "decoy");
      add("EEEEEEEK", target_score, "target");
    }
    else
    {
      add("EEEEEEEK", target_score, "target");
      add("EEEEEEEK", decoy_score, "decoy");
    }
    return ids;
  };

  // (a) strictly-better target (0.95) over an earlier decoy (0.10)
  PeptideIdentificationList ids_decoy_first = build_mixed(true, 0.10, 0.95);
  PeptideIdentificationList ids_target_first = build_mixed(false, 0.10, 0.95);
  fdr.applyBasicPeptideLevel(ids_decoy_first);
  fdr.applyBasicPeptideLevel(ids_target_first);
  const double q_decoy_first = dup_qvalue(ids_decoy_first);   // surviving target EEEEEEEK
  const double q_target_first = dup_qvalue(ids_target_first);
  TEST_NOT_EQUAL(q_decoy_first, -1.0)
  TEST_NOT_EQUAL(q_target_first, -1.0)
  TEST_REAL_SIMILAR(q_decoy_first, q_target_first)

  // (b) EQUAL scores (tie): a tie between a target and a decoy of the same sequence
  // must deterministically resolve to the target (as getPickedProteinScores_ does),
  // so the surviving target's q-value is again independent of input order.
  PeptideIdentificationList ids_tie_decoy_first = build_mixed(true, 0.95, 0.95);
  PeptideIdentificationList ids_tie_target_first = build_mixed(false, 0.95, 0.95);
  fdr.applyBasicPeptideLevel(ids_tie_decoy_first);
  fdr.applyBasicPeptideLevel(ids_tie_target_first);
  const double q_tie_decoy_first = dup_qvalue(ids_tie_decoy_first);
  const double q_tie_target_first = dup_qvalue(ids_tie_target_first);
  TEST_NOT_EQUAL(q_tie_decoy_first, -1.0)
  TEST_NOT_EQUAL(q_tie_target_first, -1.0)
  TEST_REAL_SIMILAR(q_tie_decoy_first, q_tie_target_first)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



