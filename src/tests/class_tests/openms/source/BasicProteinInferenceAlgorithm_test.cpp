// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

namespace
{
  using ID = IdentificationData;

  /**
    A consensus map with identification data: run "search" (PEPs) with database sequences P1-P3. Feature 0 links
    PEPTIDEA (P1, 0.01), feature 1 PEPTIDEB (P1 and P2, 0.02), and PEPTIDEC (P3, 0.05) is unassigned.
  */
  ConsensusMap nativeMap()
  {
    ConsensusMap map;
    auto& run = map.getIdentificationData().addRun("search");
    ID::ScoreDefinition score;
    score.name = "Posterior Error Probability";
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    ID::Database database;
    database.path = "proteins.fasta";
    const auto db = run.addDatabase(database);
    run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P1", ID::TargetDecoy::TARGET}, {db, "P2", ID::TargetDecoy::TARGET},
                                                                {db, "P3", ID::TargetDecoy::TARGET}});
    const auto source = run.addSource({});
    const auto add = [&](const std::string& sequence, const std::vector<std::string>& proteins, double pep) {
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      for (const auto& protein : proteins) match.sequence_evidence.push_back({db, protein, std::nullopt, std::nullopt, 0, 0});
      return ID::MatchReference {run.getUuid(), run.addMatch(run.addIdentification(source, ID::Observation {}), match, {pep})};
    };
    ConsensusFeature f0, f1;
    f0.setUniqueId(1);
    f0.addIDMatch(add("PEPTIDEA", {"P1"}, 0.01));
    f1.setUniqueId(2);
    f1.addIDMatch(add("PEPTIDEB", {"P1", "P2"}, 0.02));
    add("PEPTIDEC", {"P3"}, 0.05);
    map.push_back(f0);
    map.push_back(f1);
    return map;
  }

  /// accession:score of the proteins of the inference result
  std::string proteinScores(const ConsensusMap& map)
  {
    std::vector<std::string> scores;
    for (const auto& hit : map.getIdentificationData().getInferenceResults().at(0).proteins.getHits())
      scores.push_back(hit.getAccession() + ":" + StringUtils::toStr(hit.getScore()));
    return ListUtils::concatenate(scores, ",");
  }
} // namespace

START_TEST(BasicProteinInferenceAlgorithm, "$Id$")

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "false");
      bpia.setParameters(p);
      bpia.run(peps, prots);
      TEST_EQUAL(prots[0].getHits()[0].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[1].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 0.9)

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID without shared peps)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("use_shared_peptides","false");
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "false");
      bpia.setParameters(p);
      bpia.run(peps, prots);
      TEST_EQUAL(prots[0].getHits()[0].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[1].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 0.9)

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID with grouping)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "true");
      bpia.setParameters(p);
      bpia.run(peps, prots);
      TEST_EQUAL(prots[0].getHits()[0].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[1].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 0.9)

      TEST_EQUAL(prots[0].getIndistinguishableProteins().size(), 4);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[0].probability, 0.9);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[1].probability, 0.8);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[2].probability, 0.6);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[3].probability, 0.6);

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID with grouping plus resolution)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "true");
      p.setValue("greedy_group_resolution", "true");
      bpia.setParameters(p);
      bpia.run(peps, prots);

      TEST_EQUAL(prots[0].getHits().size(), 4)
      TEST_EQUAL(prots[0].getHits()[0].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[1].getScore(), 0.6)
      TEST_EQUAL(prots[0].getHits()[2].getScore(), 0.8)
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 0.9)

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 1)

      TEST_EQUAL(prots[0].getIndistinguishableProteins().size(), 3);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[0].probability, 0.9);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[1].probability, 0.8);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[2].probability, 0.6);
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID with grouping and user score)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "true");
      p.setValue("score_type", "RAW");  // should use the XTandem score meta value      
      bpia.setParameters(p);

      TEST_EQUAL(peps[0].getScoreType(), "Posterior Error Probability"); // check if main score is PEP
      bpia.run(peps, prots);
      TEST_EQUAL(peps[0].getScoreType(), "Posterior Error Probability"); // check if main score has been reset again to PEP
      
      TEST_EQUAL(prots[0].getHits()[0].getScore(), 2.5)
      TEST_EQUAL(prots[0].getHits()[1].getScore(), 2.5)
      TEST_EQUAL(prots[0].getHits()[2].getScore(), -std::numeric_limits<double>::infinity())
      TEST_EQUAL(prots[0].getHits()[3].getScore(), 5.0)
      TEST_EQUAL(prots[0].getHits()[4].getScore(), 2.5)
      TEST_EQUAL(prots[0].getHits()[5].getScore(), 10.0)

      TEST_EQUAL(prots[0].getIndistinguishableProteins().size(), 4);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[0].probability, 10.0);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[1].probability, 5.0);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[2].probability, 2.5);
      TEST_EQUAL(prots[0].getIndistinguishableProteins()[3].probability, 2.5);

      TEST_EQUAL(prots[0].getHits()[0].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[1].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[2].getMetaValue("nr_found_peptides"), 0)
      TEST_EQUAL(prots[0].getHits()[3].getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits()[4].getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits()[5].getMetaValue("nr_found_peptides"), 1)
    }
    END_SECTION

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID with grouping plus resolution and user set score type)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BasicProteinInferenceAlgorithm bpia;
      Param p = bpia.getParameters();
      p.setValue("min_peptides_per_protein", 0);
      p.setValue("annotate_indistinguishable_groups", "true");
      p.setValue("greedy_group_resolution", "true");
      p.setValue("score_type", "RAW");  // should use the XTandem score meta value
      bpia.setParameters(p);

      TEST_EQUAL(peps[0].getScoreType(), "Posterior Error Probability"); // check if main score is PEP
      bpia.run(peps, prots);
      TEST_EQUAL(peps[0].getScoreType(), "Posterior Error Probability"); // check if main score has been reset again to PEP

      TEST_EQUAL(prots[0].getHits().size(), 4)
      TEST_EQUAL(prots[0].getHits().at(0).getScore(), 2.5)
      TEST_EQUAL(prots[0].getHits().at(1).getScore(), 2.5)
      TEST_EQUAL(prots[0].getHits().at(2).getScore(), 5.0) 
      TEST_EQUAL(prots[0].getHits().at(3).getScore(), 10.0)

      TEST_EQUAL(prots[0].getHits().at(0).getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits().at(1).getMetaValue("nr_found_peptides"), 1)
      TEST_EQUAL(prots[0].getHits().at(2).getMetaValue("nr_found_peptides"), 2)
      TEST_EQUAL(prots[0].getHits().at(3).getMetaValue("nr_found_peptides"), 1)    

      TEST_EQUAL(prots[0].getIndistinguishableProteins().size(), 3);
      TEST_EQUAL(prots[0].getIndistinguishableProteins().at(0).probability, 10);
      TEST_EQUAL(prots[0].getIndistinguishableProteins().at(1).probability, 5.0);
      TEST_EQUAL(prots[0].getIndistinguishableProteins().at(2).probability, 2.5);      
    }
    END_SECTION

    START_SECTION(static void annotateIndistinguishableGroups(ProteinIdentification& proteins, const PeptideIdentificationList& peptides, Size use_top_psms = 1, bool add_singletons = true))
    {
      // P1 and P2 share every top hit; P3 has its own. The second-ranked hit of the fourth
      // spectrum and the hit from another run would each set P1 resp. P2 apart, so they must
      // only count when asked for (all hits) resp. never.
      auto make_proteins = []()
      {
        ProteinIdentification proteins;
        proteins.setIdentifier("run");
        proteins.insertHit(ProteinHit(0.9, 1, "P1", ""));
        proteins.insertHit(ProteinHit(0.8, 2, "P2", ""));
        proteins.insertHit(ProteinHit(0.7, 3, "P3", ""));
        return proteins;
      };
      auto make_hit = [](const std::string& sequence, UInt rank, const std::vector<std::string>& accessions)
      {
        PeptideHit hit(1.0 / rank, rank, 2, AASequence::fromString(sequence));
        for (const auto& accession : accessions) hit.addPeptideEvidence(PeptideEvidence(accession));
        return hit;
      };
      auto make_id = [](const std::string& run, const std::vector<PeptideHit>& hits)
      {
        PeptideIdentification id;
        id.setIdentifier(run);
        id.setHits(hits);
        return id;
      };
      PeptideIdentificationList peptides;
      peptides.push_back(make_id("run", {make_hit("PEPTIDEA", 1, {"P1", "P2"})}));
      peptides.push_back(make_id("run", {make_hit("PEPTIDEK", 1, {"P1", "P2"})}));
      peptides.push_back(make_id("run", {make_hit("PEPTIDEC", 1, {"P3"})}));
      peptides.push_back(make_id("run", {make_hit("PEPTIDED", 1, {"P3"}), make_hit("PEPTIDEE", 2, {"P1"})}));
      peptides.push_back(make_id("other run", {make_hit("PEPTIDEF", 1, {"P2"})}));

      // the groups in a canonical order: accessions sorted within, groups sorted
      auto groups = [](const ProteinIdentification& proteins)
      {
        std::vector<std::string> result;
        for (const auto& group : proteins.getIndistinguishableProteins())
        {
          std::vector<std::string> accessions = group.accessions;
          std::sort(accessions.begin(), accessions.end());
          result.push_back(ListUtils::concatenate(accessions, ","));
        }
        std::sort(result.begin(), result.end());
        return ListUtils::concatenate(result, "|");
      };

      ProteinIdentification proteins = make_proteins();
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides);
      TEST_EQUAL(groups(proteins), "P1,P2|P3")
      // nothing but the groups changes
      TEST_EQUAL(proteins.getHits().size(), 3)
      TEST_REAL_SIMILAR(proteins.getHits()[0].getScore(), 0.9)
      TEST_EQUAL(peptides[3].getHits().size(), 2)

      proteins = make_proteins();
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides, 1, false);
      TEST_EQUAL(groups(proteins), "P1,P2")

      proteins = make_proteins();
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides, 0, true);
      TEST_EQUAL(groups(proteins), "P1|P2|P3")

      proteins = make_proteins();
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides, 0, false);
      TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 0)

      // the groups replace those the run had, also from an earlier call, instead of adding to them
      ProteinIdentification::ProteinGroup stale;
      stale.accessions = {"P3", "P1"};
      proteins = make_proteins();
      proteins.getIndistinguishableProteins().push_back(stale);
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides);
      BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, peptides);
      TEST_EQUAL(groups(proteins), "P1,P2|P3")

      // on failure the run keeps its groups
      proteins = make_proteins();
      proteins.getIndistinguishableProteins().push_back(stale);
      TEST_EXCEPTION(Exception::MissingInformation,
                     BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(proteins, PeptideIdentificationList()))
      TEST_EQUAL(groups(proteins), "P1,P3")
    }
    END_SECTION

    START_SECTION((void run(ConsensusMap& cmap, bool include_unassigned) const))
    {
      // the proteins of the run are scored; the result is an inference result for the run
      auto map = nativeMap();
      BasicProteinInferenceAlgorithm bpia;
      bpia.run(map, true);
      const auto& data = map.getIdentificationData();
      TEST_EQUAL(data.getInferenceResults().size(), 1)
      ABORT_IF(data.getInferenceResults().size() != 1)
      const auto& result = data.getInferenceResults()[0];
      TEST_EQUAL(result.inputs.size(), 1)
      TEST_EQUAL(result.inputs[0].run_uuid, data.getRuns()[0].getUuid())
      TEST_EQUAL(result.proteins.getIdentifier(), "search")
      TEST_EQUAL(result.proteins.getInferenceEngine(), "TOPPProteinInference")
      TEST_EQUAL(result.proteins.getScoreType(), "Posterior Probability")
      TEST_EQUAL(proteinScores(map), "P1:0.99,P2:0.98,P3:0.95")
      TEST_EQUAL(result.proteins.getIndistinguishableProteins().size(), 3)

      // inference again replaces the result
      bpia.run(map, true);
      TEST_EQUAL(data.getInferenceResults().size(), 1)
      TEST_EQUAL(proteinScores(map), "P1:0.99,P2:0.98,P3:0.95")

      // the result is the legacy protein run
      IdentificationDataConverter::exportConsensusIDs(map);
      TEST_EQUAL(map.getProteinIdentifications().size(), 1)
      ABORT_IF(map.getProteinIdentifications().size() != 1)
      TEST_EQUAL(map.getProteinIdentifications()[0].getInferenceEngine(), "TOPPProteinInference")
      TEST_EQUAL(map.getProteinIdentifications()[0].getHits().size(), 3)

      // without the unassigned identification, P3 has no peptide (min_peptides_per_protein 1): it is removed, and
      // so is the unassigned match that referred to it
      map = nativeMap();
      bpia.run(map, false);
      TEST_EQUAL(proteinScores(map), "P1:0.99,P2:0.98")
      TEST_EQUAL(map.getIdentificationData().getRuns()[0].getNumberOfMatches(), 2)
    }
    END_SECTION

    START_SECTION([EXTRA] void run(ConsensusMap& cmap, bool include_unassigned) const with greedy group resolution)
    {
      // the shared peptide goes to the better protein: P2 loses its reference and is removed
      auto map = nativeMap();
      BasicProteinInferenceAlgorithm bpia;
      Param params = bpia.getParameters();
      params.setValue("greedy_group_resolution", "true");
      bpia.setParameters(params);
      bpia.run(map, true);
      TEST_EQUAL(proteinScores(map), "P1:0.99,P3:0.95")
      std::vector<std::string> references;
      for (const auto& query : map.getIdentificationData().getRuns()[0].getSources()[0].identifications)
      {
        for (const auto& match : query.getMatches())
        {
          const auto accessions = match.extractProteinAccessionsSet();
          references.push_back(match.representation + ":" + ListUtils::concatenate(std::vector<std::string>(accessions.begin(), accessions.end()), "+"));
        }
      }
      TEST_EQUAL(ListUtils::concatenate(references, ","), "PEPTIDEA:P1,PEPTIDEB:P1,PEPTIDEC:P3")
    }
    END_SECTION

    START_SECTION([EXTRA] void run(ConsensusMap& cmap, ProteinIdentification& prot_id, bool include_unassigned) const)
    {
      BasicProteinInferenceAlgorithm bpia;
      // with identification data: the run gets the proteins of the result
      auto map = nativeMap();
      ProteinIdentification proteins;
      bpia.run(map, proteins, true);
      TEST_EQUAL(proteins.getHits().size(), 3)
      TEST_EQUAL(proteins.getScoreType(), "Posterior Probability")

      // with peptide identifications: the only protein run of the map gets the result
      map = nativeMap();
      IdentificationDataConverter::exportConsensusIDs(map);
      ProteinIdentification other;
      TEST_EXCEPTION(Exception::InvalidParameter, bpia.run(map, other, true))
      auto& run = map.getProteinIdentifications()[0];
      bpia.run(map, run, true);
      TEST_EQUAL(&map.getProteinIdentifications()[0] == &run, true)
      TEST_EQUAL(run.getScoreType(), "Posterior Probability")
      TEST_EQUAL(run.getInferenceEngine(), "TOPPProteinInference")
      TEST_EQUAL(map.getIdentificationData().empty(), true)
    }
    END_SECTION

    START_SECTION([EXTRA] void run(ConsensusMap& cmap, bool include_unassigned) const needs one protein run)
    {
      // two runs in two protein runs; merging pools them
      auto map = nativeMap();
      auto second = nativeMap();
      map.getIdentificationData().merge(second.getIdentificationData());
      BasicProteinInferenceAlgorithm bpia;
      TEST_EXCEPTION(Exception::InvalidParameter, bpia.run(map, true))
    }
    END_SECTION

END_TEST
