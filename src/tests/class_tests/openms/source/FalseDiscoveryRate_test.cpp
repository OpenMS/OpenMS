// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/NativeIdentificationTest.h>

#include <limits>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
// A consensus map of the PSMs of FalseDiscoveryRate_OMSSA.idXML: two of three are assigned (to features of one PSM),
// the others unassigned.
ConsensusMap omssaMap(const vector<ProteinIdentification>& prot_ids, const PeptideIdentificationList& pep_ids)
{
  ConsensusMap cmap;
  MapPeptideIdentifications ids;
  ids.proteins = prot_ids;
  for (Size i = 0; i < pep_ids.size(); ++i)
  {
    if (i % 3 == 0)
    {
      ids.unassigned.push_back(pep_ids[i]);
      continue;
    }
    ConsensusFeature feature;
    feature.setUniqueId(i + 1);
    ids.features[{cmap.size()}].push_back(pep_ids[i]);
    cmap.push_back(feature);
  }
  IdentificationDataConverter::importMapIdentifications(cmap, ids);
  return cmap;
}
/// The peptide identifications of the features of @p cmap
PeptideIdentificationList assignedPSMs(const ConsensusMap& cmap)
{
  PeptideIdentificationList result;
  for (const auto& [path, peptides] : IdentificationDataConverter::exportMapIdentifications(cmap).features)
    result.insert(result.end(), peptides.begin(), peptides.end());
  return result;
}
/// The unassigned peptide identifications of @p cmap
PeptideIdentificationList unassignedPSMs(const ConsensusMap& cmap)
{
  return IdentificationDataConverter::exportMapIdentifications(cmap).unassigned;
}
/// A consensus map with a feature for each identification of @p psms, of the protein run "run"
ConsensusMap psmMap(const PeptideIdentificationList& psms)
{
  ConsensusMap cmap;
  MapPeptideIdentifications ids;
  ProteinIdentification run;
  run.setIdentifier("run");
  ids.proteins = {run};
  for (const auto& psm : psms)
  {
    ConsensusFeature feature;
    feature.setUniqueId(cmap.size() + 1);
    ids.features[{cmap.size()}].push_back(psm);
    cmap.push_back(feature);
  }
  IdentificationDataConverter::importMapIdentifications(cmap, ids);
  return cmap;
}
/// The peptide identifications of the features of @p cmap, then the unassigned ones
PeptideIdentificationList allPSMs(const ConsensusMap& cmap)
{
  auto result = assignedPSMs(cmap);
  const auto unassigned = unassignedPSMs(cmap);
  result.insert(result.end(), unassigned.begin(), unassigned.end());
  return result;
}
} // namespace

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

  // (7) A group without a score (NaN) that takes no part in the competition gets the worst
  // calibrated q-value, in both score directions. Fixture (1) with {B}'s group score NaN: 2/3.
  {
    for (bool higher_better : {false, true})
    {
      const double s = higher_better ? -1.0 : 1.0;
      ProteinIdentification run = picked_run(higher_better,
        {{"A", s * 0.01}, {"DECOY_A", s * 0.02}, {"B", s * 0.05}, {"DECOY_B", s * 0.04}, {"C", s * 0.02}},
        {{s * 0.01, {"A"}}, {std::numeric_limits<double>::quiet_NaN(), {"B"}}, {s * 0.02, {"C"}}, {s * 0.04, {"DECOY_B"}}});
      fdr.applyPickedProteinFDR(run, "DECOY_", true);
      TEST_REAL_SIMILAR(group_qvalues(run)[1], 2.0 / 3.0)
    }
  }

  // (8) Prefix mode: a target accession that ends with the decoy string is still a target.
  // {XDECOY_} 0.9 T, {C} 0.5 T: q = 1/3 for both. Before the fix the decoy string was also removed
  // as a suffix, and the lookup of "X" threw std::out_of_range.
  {
    ProteinIdentification run = picked_run(true, {{"XDECOY_", 0.9}, {"C", 0.5}}, {{0.9, {"XDECOY_"}}, {0.5, {"C"}}});
    fdr.applyPickedProteinFDR(run, "DECOY_", true);
    vector<double> qs = group_qvalues(run);
    TEST_REAL_SIMILAR(qs[0], 1.0 / 3.0)
    TEST_REAL_SIMILAR(qs[1], 1.0 / 3.0)
  }

  // (9) Two equal-scoring hits of the same decoy (no target) stay a decoy: {DECOY_A} 0.9 D, {C} 0.5 T.
  // FDR 0.9 -> (1+1)/(0+1) = 2, 0.5 -> (1+1)/(1+1) = 1; q = 1 for both. Before the fix the tie made
  // the pair's winner a target, {DECOY_A} did not vote and both groups got 1/2.
  {
    ProteinIdentification run = picked_run(true, {{"DECOY_A", 0.9}, {"DECOY_A", 0.9}, {"C", 0.5}},
                                           {{0.9, {"DECOY_A"}}, {0.5, {"C"}}});
    fdr.applyPickedProteinFDR(run, "DECOY_", true);
    vector<double> qs = group_qvalues(run);
    TEST_REAL_SIMILAR(qs[0], 1.0)
    TEST_REAL_SIMILAR(qs[1], 1.0)
    // The tie also decides the protein hits' FDR (same picked scores, decoys removed by default):
    // C gets q = 1. Before the fix A counted as a target winner at 0.9 and C got 1/3.
    ABORT_IF(run.getHits().size() != 1)
    TEST_STRING_EQUAL(run.getHits()[0].getAccession(), "C")
    TEST_REAL_SIMILAR(run.getHits()[0].getScore(), 1.0)
  }
}
END_SECTION

START_SECTION((void applyBasic(ProteinIdentification & id, bool groups_too = true)))
{
  // Protein-group q-values for lower-is-better scores. Scores within 1e-12 count as one score; its
  // key is the best one. Groups: {T1} 0.01 T, {T2} 0.01 + 5e-13 T, {DECOY_D} 0.02 D, {T3} 0.03 T.
  // FDR (D+1)/(T+1): 0.01 -> 1/3, 0.02 -> 2/3, 0.03 -> 2/4; q: 1/3, 1/2, 1/2. {T2} gets the q-value of
  // its own score 0.01, 1/3; before the fix it got the one of the next worse score 0.02 (1/2).
  ProteinIdentification run;
  run.setScoreType("score");
  run.setHigherScoreBetter(false);
  for (const auto& [acc, score] : vector<pair<string, double>>{{"T1", 0.01}, {"T2", 0.01 + 5e-13}, {"DECOY_D", 0.02}, {"T3", 0.03}})
  {
    ProteinHit hit(score, 1, acc, "");
    hit.setMetaValue("target_decoy", acc.rfind("DECOY_", 0) == 0 ? "decoy" : "target");
    run.getHits().push_back(hit);
    ProteinIdentification::ProteinGroup grp;
    grp.probability = score;
    grp.accessions = {acc};
    run.getIndistinguishableProteins().push_back(grp);
  }
  FalseDiscoveryRate fdr;
  fdr.applyBasic(run, true);
  TOLERANCE_ABSOLUTE(1e-6)
  const auto& grps = run.getIndistinguishableProteins();
  TEST_REAL_SIMILAR(grps[0].probability, 1.0 / 3.0)
  TEST_REAL_SIMILAR(grps[1].probability, 1.0 / 3.0)
  TEST_REAL_SIMILAR(grps[2].probability, 0.5)
  TEST_REAL_SIMILAR(grps[3].probability, 0.5)
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

START_SECTION([EXTRA] protein FDR without indistinguishable groups skips the group pass)
{
  // Inference without group annotation (e.g. for '-protein_quantification strictly_unique_peptides')
  // leaves no groups. The group-level FDR then has nothing to score: it must be skipped, not warn about
  // missing target-decoy annotation, while the protein-level FDR still runs.
  auto make_run = []()
  {
    ProteinIdentification run;
    run.setScoreType("score");
    run.setHigherScoreBetter(true);
    const std::vector<std::pair<std::string, double>> proteins = {
      {"P1", 0.9}, {"P2", 0.8}, {"decoy_P3", 0.7}, {"P3", 0.6}, {"decoy_P1", 0.5}};
    for (const auto& [accession, score] : proteins)
    {
      ProteinHit hit;
      hit.setAccession(accession);
      hit.setScore(score);
      hit.setMetaValue("target_decoy", accession.rfind("decoy_", 0) == 0 ? "decoy" : "target");
      run.getHits().push_back(hit);
    }
    return run;
  };

  FalseDiscoveryRate fdr;
  for (bool picked : {false, true})
  {
    ProteinIdentification run = make_run();
    // Start from a clean log: lines written with '\n' wait in the buffer until the next sync (which
    // remove() does into the streams still attached), and a line repeating a cached one is not passed
    // on at all. An earlier section emits this very warning (a group FDR finding no competing group),
    // so without flushing and clearing the cache the capture would see its line, or miss ours.
    OPENMS_LOG_WARN.flush();
    OPENMS_LOG_WARN.rdbuf()->clearCache();
    std::ostringstream captured;
    OPENMS_LOG_WARN.insert(captured);
    try
    {
      if (picked) { fdr.applyPickedProteinFDR(run, "decoy_", true); }
      else { fdr.applyBasic(run); }
    }
    catch (...)
    {
      OPENMS_LOG_WARN.remove(captured);
      throw;
    }
    OPENMS_LOG_WARN.remove(captured);

    TEST_TRUE(captured.str().find("No scores extracted") == std::string::npos)
    TEST_TRUE(run.getIndistinguishableProteins().empty())
    TEST_EQUAL(run.getScoreType(), "q-value")
  }
}
END_SECTION

START_SECTION([EXTRA] FDRs and q-values on a ConsensusMap are lower-is-better whatever the search score)
{
  // The q-values are lower-is-better even after a higher-is-better search score (a hyperscore), so
  // filtering by them keeps the best hits.
  auto make_psms = []()
  {
    PeptideIdentificationList ids;
    const std::vector<std::tuple<std::string, double, std::string>> psms = {
      {"AAAAAAAK", 90.0, "target"}, {"CCCCCCCK", 80.0, "target"}, {"DDDDDDDK", 70.0, "target"},
      {"FFFFFFFK", 50.0, "decoy"}, {"EEEEEEEK", 40.0, "target"}, {"HHHHHHHK", 30.0, "decoy"}};
    for (const auto& [sequence, score, target_decoy] : psms)
    {
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(sequence));
      hit.setScore(score);
      hit.setCharge(2);
      hit.setMetaValue("target_decoy", target_decoy);
      PeptideIdentification id;
      id.setIdentifier("run");
      id.setScoreType("hyperscore");
      id.setHigherScoreBetter(true);
      id.setHits({hit});
      ids.push_back(id);
    }
    return ids;
  };

  // q-value of the hit with this sequence (-1 if it is gone)
  auto qvalue = [](const PeptideIdentificationList& ids, const std::string& sequence)
  {
    for (const auto& id : ids)
    {
      if (!id.getHits().empty() && id.getHits()[0].getSequence().toUnmodifiedString() == sequence)
      {
        return id.getHits()[0].getScore();
      }
    }
    return -1.0;
  };

  for (bool peptide_level : {false, true})
  {
    ConsensusMap cmap = psmMap(make_psms());
    FalseDiscoveryRate fdr;
    if (peptide_level) { fdr.applyBasicPeptideLevel(cmap, true); }
    else { fdr.applyBasic(cmap, true); }

    TEST_EQUAL(cmap.getIdentificationData().getPrimaryScoreDefinition()->higher_better, false)
    PeptideIdentificationList ids = assignedPSMs(cmap);
    for (const auto& id : ids)
    {
      if (id.getHits().empty()) { continue; } // a decoy, removed
      TEST_EQUAL(id.isHigherScoreBetter(), false)
    }

    // the best target has a lower q-value than the target ranked below a decoy ...
    const double q_best = qvalue(ids, "AAAAAAAK");
    const double q_worse = qvalue(ids, "EEEEEEEK");
    TEST_EQUAL(q_best < q_worse, true)

    // ... so filtering at the best target's q-value keeps it and removes the other one
    IDFilter::filterHitsByScore(ids, q_best);
    TEST_REAL_SIMILAR(qvalue(ids, "AAAAAAAK"), q_best)
    TEST_REAL_SIMILAR(qvalue(ids, "EEEEEEEK"), -1.0)
  }

  // the same for the PeptideIdentificationList overload of the peptide-level FDR
  PeptideIdentificationList ids = make_psms();
  FalseDiscoveryRate fdr;
  fdr.applyBasicPeptideLevel(ids);
  for (const auto& id : ids)
  {
    if (id.getHits().empty()) { continue; }
    TEST_EQUAL(id.isHigherScoreBetter(), false)
  }
}
END_SECTION

START_SECTION([EXTRA] peptide-level FDR scores the best hit of an identification and the hits tied with it)
{
  // Only the best hit of each identification, and a hit tied with it whose peptide has a q-value,
  // get a peptide q-value. Any other hit would keep its search score (0.1) under the q-value score
  // type and pass a q-value filter at 0.25. The best hit is kept even when it is not stored first.
  auto make_psms = []()
  {
    PeptideIdentificationList ids;
    const std::vector<std::tuple<std::string, double, std::string>> psms = {
      {"AAAAAAAK", 90.0, "target"}, {"CCCCCCCK", 80.0, "target"}, {"FFFFFFFK", 50.0, "decoy"}};
    for (const auto& [sequence, score, target_decoy] : psms)
    {
      std::vector<PeptideHit> hits(1);
      hits[0].setSequence(AASequence::fromString(sequence));
      hits[0].setScore(score);
      hits[0].setCharge(2);
      hits[0].setMetaValue("target_decoy", target_decoy);
      if (sequence == "AAAAAAAK")
      {
        hits.push_back(hits[0]); // a much worse candidate for the same spectrum
        hits[1].setSequence(AASequence::fromString("KKKKKKKK"));
        hits[1].setScore(0.1);
        hits.push_back(hits[0]); // tied, and the best hit of another identification: kept
        hits[2].setSequence(AASequence::fromString("CCCCCCCK"));
        hits.push_back(hits[0]); // tied, but its peptide has no q-value: removed
        hits[3].setSequence(AASequence::fromString("EEEEEEEK"));
      }
      if (sequence == "CCCCCCCK") // the same, but stored before the best hit
      {
        hits.insert(hits.begin(), hits[0]);
        hits[0].setSequence(AASequence::fromString("KKKKKKKK"));
        hits[0].setScore(0.1);
      }
      PeptideIdentification id;
      id.setIdentifier("run");
      id.setScoreType("hyperscore");
      id.setHigherScoreBetter(true);
      id.setHits(hits);
      ids.push_back(id);
    }
    return ids;
  };

  auto check = [](const PeptideIdentificationList& ids)
  {
    for (const auto& id : ids)
    {
      if (id.getHits().empty()) { continue; } // a decoy, removed
      std::vector<std::string> sequences;
      for (const auto& hit : id.getHits())
      {
        sequences.push_back(hit.getSequence().toUnmodifiedString());
        TEST_EQUAL(hit.getScore() <= 1.0, true) // a q-value, not the search score
      }
      if (sequences[0] == "AAAAAAAK") { TEST_EQUAL(ListUtils::concatenate(sequences, ","), "AAAAAAAK,CCCCCCCK") }
      else { TEST_EQUAL(sequences.size(), 1) }
      TEST_NOT_EQUAL(sequences[0], "KKKKKKKK")
    }
  };

  FalseDiscoveryRate fdr;
  ConsensusMap cmap = psmMap(make_psms());
  fdr.applyBasicPeptideLevel(cmap, true);
  check(assignedPSMs(cmap));

  PeptideIdentificationList ids = make_psms();
  fdr.applyBasicPeptideLevel(ids);
  check(ids);
}
END_SECTION

delete ptr;

START_SECTION((void applyBasic(ConsensusMap & cmap, bool include_unassigned_peptides)))
{
  using Internal::ClassTest::NativeIdentificationDetail::peptideDifference;
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_OMSSA.idXML"), prot_ids, pep_ids);
  FalseDiscoveryRate fdr;

  // With the unassigned PSMs, the PSMs of the map get what they get together.
  auto cmap = omssaMap(prot_ids, pep_ids);
  const auto assigned_count = assignedPSMs(cmap).size();
  auto expected = allPSMs(cmap);
  fdr.applyBasic(prot_ids, expected);
  fdr.applyBasic(cmap, true);
  TEST_EQUAL(cmap.getIdentificationData().getPrimaryScoreDefinition()->name, "q-value")
  TEST_EQUAL(unassignedPSMs(cmap).size(), pep_ids.size() - assigned_count) // the unassigned ones stay
  auto actual = allPSMs(cmap);
  TEST_EQUAL(peptideDifference(expected, actual), "")
  TEST_EQUAL(actual[0].getScoreType(), "q-value")
  TEST_TRUE(actual[0].getHits()[0].metaValueExists("OMSSA_score"))

  // Without, the FDR is that of the assigned PSMs; the unassigned ones get the q-value of their score, keep their decoys.
  cmap = omssaMap(prot_ids, pep_ids);
  expected = assignedPSMs(cmap);
  TEST_EQUAL(expected.size(), assigned_count)
  fdr.applyBasic(prot_ids, expected);
  fdr.applyBasic(cmap, false);
  TEST_EQUAL(peptideDifference(expected, assignedPSMs(cmap)), "")
  std::map<double, double> q_values; // OMSSA score -> q-value of the assigned PSMs
  for (const auto& id : expected)
    for (const auto& hit : id.getHits())
      q_values[hit.getMetaValue("OMSSA_score")] = hit.getScore();
  Size decoys = 0, shared_scores = 0;
  for (const auto& id : unassignedPSMs(cmap))
  {
    TEST_EQUAL(id.getScoreType(), "q-value")
    TEST_EQUAL(id.isHigherScoreBetter(), false)
    for (const auto& hit : id.getHits())
    {
      TEST_EQUAL(hit.getScore() >= 0.0 && hit.getScore() <= 1.0, true)
      if (hit.getMetaValue("target_decoy") == "decoy") ++decoys;
      const auto found = q_values.find(hit.getMetaValue("OMSSA_score"));
      if (found == q_values.end()) continue;
      ++shared_scores;
      TEST_REAL_SIMILAR(hit.getScore(), found->second)
    }
  }
  TEST_NOT_EQUAL(decoys, 0)
  TEST_NOT_EQUAL(shared_scores, 0)
}
END_SECTION

START_SECTION((double rocN(const ConsensusMap& ids, Size fp_cutoff, bool include_unassigned_peptides) const))
{
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_OMSSA.idXML"), prot_ids, pep_ids);
  FalseDiscoveryRate fdr;
  const auto cmap = omssaMap(prot_ids, pep_ids);
  const auto all = allPSMs(cmap);
  const auto assigned = assignedPSMs(cmap);
  TEST_REAL_SIMILAR(fdr.rocN(cmap, 50, true), fdr.rocN(all, 50))
  TEST_REAL_SIMILAR(fdr.rocN(cmap, 50, false), fdr.rocN(assigned, 50))
  TEST_REAL_SIMILAR(fdr.rocN(cmap, 50, prot_ids[1].getIdentifier(), true), fdr.rocN(all, 50, prot_ids[1].getIdentifier()))
}
END_SECTION

START_SECTION((void applyBasicPeptideLevel(ConsensusMap & map, bool include_unassigned)))
{
  vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FalseDiscoveryRate_OMSSA.idXML"), prot_ids, pep_ids);
  FalseDiscoveryRate fdr;
  // With the unassigned PSMs, (best) hits get the q-values of their peptides as with the PSMs together.
  auto cmap = omssaMap(prot_ids, pep_ids);
  auto expected = allPSMs(cmap);
  fdr.applyBasicPeptideLevel(expected);
  fdr.applyBasicPeptideLevel(cmap, true);
  const auto actual = allPSMs(cmap);
  ABORT_IF(actual.size() != expected.size())
  Size compared = 0;
  for (Size i = 0; i < actual.size(); ++i)
  {
    TEST_EQUAL(actual[i].getHits().size(), expected[i].getHits().size())
    if (expected[i].getHits().empty() || actual[i].getHits().empty()) continue;
    TEST_EQUAL(actual[i].getScoreType(), expected[i].getScoreType())
    TEST_REAL_SIMILAR(actual[i].getHits()[0].getScore(), expected[i].getHits()[0].getScore())
    TEST_EQUAL(actual[i].getHits()[0].getMetaValue("OMSSA"), expected[i].getHits()[0].getMetaValue("OMSSA"))
    ++compared;
  }
  TEST_NOT_EQUAL(compared, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



