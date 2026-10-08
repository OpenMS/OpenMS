// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;

namespace
{
  using ID = IdentificationData;

  /**
    A consensus map with identification data: run "search" (posterior error probabilities) with database sequences
    P1-P4. Features link PEPTIDEA (P1, 0.01), PEPTIDEB (P1 and P2, 0.02) and PEPTIDEC (P3, 0.05); PEPTIDED (P4, 0.1)
    is unassigned.
  */
  ConsensusMap nativeMap(const std::string& score_name = "Posterior Error Probability")
  {
    ConsensusMap map;
    auto& run = map.getIdentificationData().addRun("search");
    ID::ScoreDefinition score;
    score.name = score_name;
    score.higher_better = false;
    run.setPrimaryScore(run.addScore(score));
    ID::Database database;
    database.path = "proteins.fasta";
    const auto db = run.addDatabase(database);
    run.setDatabaseSequences(std::vector<ID::DatabaseSequence> {{db, "P1", ID::TargetDecoy::TARGET}, {db, "P2", ID::TargetDecoy::TARGET},
                                                                {db, "P3", ID::TargetDecoy::TARGET}, {db, "P4", ID::TargetDecoy::TARGET}});
    const auto source = run.addSource({});
    const auto add = [&](const std::string& sequence, const std::vector<std::string>& proteins, double pep) {
      ID::MatchData match;
      match.representation = sequence;
      match.charge = 2;
      match.target_decoy = ID::TargetDecoy::TARGET;
      for (const auto& protein : proteins) match.sequence_evidence.push_back({db, protein, std::nullopt, std::nullopt, 0, 0});
      return ID::MatchReference {run.getUuid(), run.addMatch(run.addIdentification(source, ID::Observation {}), match, {pep})};
    };
    ConsensusFeature f0, f1, f2;
    f0.setUniqueId(1);
    f0.addIDMatch(add("PEPTIDEA", {"P1"}, 0.01));
    f1.setUniqueId(2);
    f1.addIDMatch(add("PEPTIDEB", {"P1", "P2"}, 0.02));
    f2.setUniqueId(3);
    f2.addIDMatch(add("PEPTIDEC", {"P3"}, 0.05));
    add("PEPTIDED", {"P4"}, 0.1);
    map.push_back(f0);
    map.push_back(f1);
    map.push_back(f2);
    return map;
  }

  /// Model parameters as in the TOPP tests, so there is no grid search
  void setModel(BayesianProteinInferenceAlgorithm& bpia, const std::string& cutoff = "")
  {
    Param p = bpia.getParameters();
    p.setValue("model_parameters:prot_prior", 0.7);
    p.setValue("model_parameters:pep_emission", 0.1);
    p.setValue("model_parameters:pep_spurious_emission", 0.001);
    if (!cutoff.empty()) p.setValue("psm_probability_cutoff", std::stod(cutoff));
    bpia.setParameters(p);
  }

  /// The inference of the peptide identification path on the identifications of @p map that features link
  void reference(const ConsensusMap& map, bool greedy, std::vector<ProteinIdentification>& proteins, PeptideIdentificationList& peptides)
  {
    IdentificationDataConverter::exportIDs(map.getIdentificationData(), proteins, peptides);
    PeptideIdentificationList assigned;
    for (const auto& peptide : peptides)
    {
      if (peptide.getHits()[0].getSequence().toString() != "PEPTIDED") assigned.push_back(peptide);
    }
    peptides = assigned;
    BayesianProteinInferenceAlgorithm bpia;
    setModel(bpia);
    bpia.inferPosteriorProbabilities(proteins, peptides, greedy);
  }

  /// accession -> score of the proteins of the inference result
  std::map<std::string, double> proteinScores(const ConsensusMap& map)
  {
    std::map<std::string, double> scores;
    for (const auto& hit : map.getIdentificationData().getInferenceResults().at(0).proteins.getHits()) scores[hit.getAccession()] = hit.getScore();
    return scores;
  }

  /// sequence -> primary score of the matches
  std::map<std::string, double> matchScores(const ConsensusMap& map)
  {
    std::map<std::string, double> scores;
    for (const auto& run : map.getIdentificationData().getRuns())
    {
      const auto score = run.bindScore(*run.getPrimaryScore());
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches()) scores[match.representation] = *score(match);
    }
    return scores;
  }

  /// sequence -> referenced proteins of the matches
  std::map<std::string, std::set<std::string>> matchProteins(const ConsensusMap& map)
  {
    std::map<std::string, std::set<std::string>> proteins;
    for (const auto& run : map.getIdentificationData().getRuns())
      for (const auto& source : run.getSources())
        for (const auto& query : source.identifications)
          for (const auto& match : query.getMatches()) proteins[match.representation] = match.extractProteinAccessionsSet();
    return proteins;
  }
} // namespace

START_TEST(BayesianProteinInferenceAlgorithm, "$Id$")

    START_SECTION(BayesianProteinInferenceAlgorithm on Protein Peptide ID)
    {
      vector<ProteinIdentification> prots;
      PeptideIdentificationList peps;
      IdXMLFile idf;
      idf.load(OPENMS_GET_TEST_DATA_PATH("newMergerTest_out.idXML"),prots,peps);
      BayesianProteinInferenceAlgorithm bpia;
      bpia.inferPosteriorProbabilities(prots,peps,false);
    }
    END_SECTION

    TOLERANCE_ABSOLUTE(0.002)
    TOLERANCE_RELATIVE(1.002)
    START_SECTION(BayesianProteinInferenceAlgorithm test)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("update_PSM_probabilities", "false");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 9)
          TEST_EQUAL(peps[0].getHits()[0].getScore(), 0.6)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.624641)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.648346)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.9);
          p.setValue("model_parameters:prot_prior", 0.3);
          p.setValue("model_parameters:pep_spurious_emission", 0.1);
          p.setValue("model_parameters:pep_prior", 0.3);
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 9)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.827132)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.755653)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.580705)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 filter)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.9);
          p.setValue("model_parameters:prot_prior", 0.3);
          p.setValue("model_parameters:pep_spurious_emission", 0.1);
          p.setValue("model_parameters:pep_prior", 0.3);
          p.setValue("psm_probability_cutoff",0.61);
          //TODO setParams needs to update the filter function or we need to make a member.
          //p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 8)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.77821544)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.787325)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.609742)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 regularize)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.9);
          p.setValue("model_parameters:prot_prior", 0.3);
          p.setValue("model_parameters:pep_spurious_emission", 0.1);
          p.setValue("model_parameters:pep_prior", 0.3);
          //p.setValue("loopy_belief_propagation:p_norm_inference", -1.)
          p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 9)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.779291)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.684165)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.458033)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 regularize max-product)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.9);
          p.setValue("model_parameters:prot_prior", 0.3);
          p.setValue("model_parameters:pep_spurious_emission", 0.1);
          p.setValue("model_parameters:pep_prior", 0.3);
          p.setValue("loopy_belief_propagation:p_norm_inference", -1.);
          p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 9)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.83848989)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(),   0.784666)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(),  0.548296)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 max-product)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.9);
          p.setValue("model_parameters:prot_prior", 0.3);
          p.setValue("model_parameters:pep_spurious_emission", 0.1);
          p.setValue("model_parameters:pep_prior", 0.3);
          p.setValue("loopy_belief_propagation:p_norm_inference", -1.);
          //p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 9)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.9117111)
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.879245)
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.708133)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 super-easy)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_2_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.7);
          p.setValue("model_parameters:prot_prior", 0.5);
          p.setValue("model_parameters:pep_spurious_emission", 0.0);
          p.setValue("model_parameters:pep_prior", 0.5);
          p.setValue("loopy_belief_propagation:dampening_lambda", 0.0);
          p.setValue("loopy_belief_propagation:p_norm_inference", 1.);
          //p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 3)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.843211)
          TEST_REAL_SIMILAR(peps[1].getHits()[0].getScore(), 0.944383)
          TEST_REAL_SIMILAR(peps[2].getHits()[0].getScore(), 0.701081)
          std::cout << prots[0].getHits()[0].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.883060)
          std::cout << prots[0].getHits()[1].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.519786)
          std::cout << prots[0].getHits()[2].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[2].getScore(), 0.775994)
        }
    END_SECTION

    START_SECTION(BayesianProteinInferenceAlgorithm test2 mini-loop)
        {
          vector<ProteinIdentification> prots;
          PeptideIdentificationList peps;
          IdXMLFile idf;
          idf.load(OPENMS_GET_TEST_DATA_PATH("BayesianProteinInference_3_test.idXML"),prots,peps);
          BayesianProteinInferenceAlgorithm bpia;
          Param p = bpia.getParameters();
          p.setValue("model_parameters:pep_emission", 0.7);
          p.setValue("model_parameters:prot_prior", 0.5);
          p.setValue("model_parameters:pep_spurious_emission", 0.0);
          p.setValue("model_parameters:pep_prior", 0.5);
          p.setValue("loopy_belief_propagation:dampening_lambda", 0.0);
          p.setValue("loopy_belief_propagation:p_norm_inference", 1.);
          //p.setValue("model_parameters:regularize","true");
          bpia.setParameters(p);
          bpia.inferPosteriorProbabilities(prots,peps,false);
          TEST_EQUAL(peps.size(), 3)
          TEST_REAL_SIMILAR(peps[0].getHits()[0].getScore(), 0.934571)
          TEST_REAL_SIMILAR(peps[1].getHits()[0].getScore(), 0.944383)
          TEST_REAL_SIMILAR(peps[2].getHits()[0].getScore(), 0.701081)
          std::cout << prots[0].getHits()[0].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[0].getScore(), 0.675421)
          std::cout << prots[0].getHits()[1].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[1].getScore(), 0.675421)
          std::cout << prots[0].getHits()[2].getAccession() << std::endl;
          TEST_REAL_SIMILAR(prots[0].getHits()[2].getScore(), 0.775994)
        }
    END_SECTION

    START_SECTION((void inferPosteriorProbabilities(ConsensusMap& cmap, bool greedy_group_resolution, std::optional<const ExperimentalDesign> exp_des)))
    {
      // inference on identification data agrees with inference on the (exported) peptide identifications
      auto map = nativeMap();
      BayesianProteinInferenceAlgorithm bpia;
      setModel(bpia);
      bpia.inferPosteriorProbabilities(map, false);
      std::vector<ProteinIdentification> ref_proteins;
      PeptideIdentificationList ref_peptides;
      reference(nativeMap(), false, ref_proteins, ref_peptides);

      const auto& data = map.getIdentificationData();
      TEST_EQUAL(data.getInferenceResults().size(), 1)
      ABORT_IF(data.getInferenceResults().size() != 1)
      const auto& result = data.getInferenceResults()[0];
      TEST_EQUAL(result.inputs.size(), 1)
      TEST_EQUAL(result.inputs[0].run_uuid, data.getRuns()[0].getUuid())
      TEST_EQUAL(result.proteins.getInferenceEngine(), "Epifany")
      TEST_EQUAL(result.proteins.getScoreType(), "Posterior Probability")
      // the protein of the unassigned PSM is not inferred: it is kept with score 0
      const auto scores = proteinScores(map);
      TEST_EQUAL(scores.size(), 4)
      TEST_EQUAL(ref_proteins[0].getHits().size(), 3)
      for (const auto& hit : ref_proteins[0].getHits())
      {
        TEST_REAL_SIMILAR(scores.at(hit.getAccession()), hit.getScore())
      }
      TEST_REAL_SIMILAR(scores.at("P4"), 0.0)
      TEST_EQUAL(scores.at("P1") > scores.at("P2"), true)
      TEST_EQUAL(result.proteins.getIndistinguishableProteins().size(), 4)

      // the PSMs have posterior probabilities, updated by inference; the unassigned one is converted only
      TEST_EQUAL(data.getPrimaryScoreDefinition()->name, "Posterior Probability")
      TEST_EQUAL(data.getPrimaryScoreDefinition()->higher_better, true)
      const auto psms = matchScores(map);
      TEST_EQUAL(psms.size(), 4)
      TEST_EQUAL(ref_peptides.size(), 3)
      for (const auto& peptide : ref_peptides)
      {
        TEST_REAL_SIMILAR(psms.at(peptide.getHits()[0].getSequence().toString()), peptide.getHits()[0].getScore())
      }
      TEST_REAL_SIMILAR(psms.at("PEPTIDED"), 0.9)
      TEST_NOT_EQUAL(psms.at("PEPTIDEA"), 0.99)

      // the result is the legacy protein run
      IdentificationDataConverter::exportConsensusIDs(map);
      TEST_EQUAL(map.getProteinIdentifications().size(), 1)
      ABORT_IF(map.getProteinIdentifications().size() != 1)
      TEST_EQUAL(map.getProteinIdentifications()[0].getInferenceEngine(), "Epifany")
      TEST_EQUAL(map.getProteinIdentifications()[0].getHits().size(), 4)
      TEST_EQUAL(map[0].getPeptideIdentifications()[0].getScoreType(), "Posterior Probability")

      // a map with peptide identifications is converted for inference: same result
      auto legacy = nativeMap();
      IdentificationDataConverter::exportConsensusIDs(legacy);
      bpia.inferPosteriorProbabilities(legacy, false);
      TEST_EQUAL(legacy.getIdentificationData().empty(), true)
      TEST_EQUAL(legacy.getProteinIdentifications().size(), 1)
      ABORT_IF(legacy.getProteinIdentifications().size() != 1)
      for (const auto& hit : legacy.getProteinIdentifications()[0].getHits())
      {
        TEST_REAL_SIMILAR(scores.at(hit.getAccession()), hit.getScore())
      }
      TEST_REAL_SIMILAR(legacy[0].getPeptideIdentifications()[0].getHits()[0].getScore(), psms.at("PEPTIDEA"))
    }
    END_SECTION

    START_SECTION([EXTRA] inferPosteriorProbabilities(ConsensusMap&) filters PSMs by probability)
    {
      // PEPTIDEC (probability 0.95) is below the cutoff: its match is erased, and P3 is left without evidence
      auto map = nativeMap();
      BayesianProteinInferenceAlgorithm bpia;
      setModel(bpia, "0.96");
      bpia.inferPosteriorProbabilities(map, false);
      const auto psms = matchScores(map);
      TEST_EQUAL(psms.size(), 2)
      TEST_EQUAL(psms.count("PEPTIDEC"), 0)
      TEST_EQUAL(psms.count("PEPTIDED"), 0)
      const auto scores = proteinScores(map);
      TEST_EQUAL(scores.size(), 2)
      TEST_EQUAL(scores.count("P3"), 0)
      // the feature without its match links no identification
      TEST_EQUAL(map[2].getIDMatches().empty(), true)
      TEST_EQUAL(map[2].getIDQueries().empty(), true)
    }
    END_SECTION

    START_SECTION([EXTRA] inferPosteriorProbabilities(ConsensusMap&) with greedy group resolution)
    {
      // the shared peptide keeps the protein that the peptide identification path keeps
      auto map = nativeMap();
      BayesianProteinInferenceAlgorithm bpia;
      setModel(bpia);
      bpia.inferPosteriorProbabilities(map, true);
      std::vector<ProteinIdentification> ref_proteins;
      PeptideIdentificationList ref_peptides;
      reference(nativeMap(), true, ref_proteins, ref_peptides);
      const auto proteins = matchProteins(map);
      for (const auto& peptide : ref_peptides)
      {
        const auto& hit = peptide.getHits()[0];
        const auto expected = hit.extractProteinAccessionsSet();
        const auto& actual = proteins.at(hit.getSequence().toString());
        TEST_EQUAL(ListUtils::concatenate(std::vector<std::string>(actual.begin(), actual.end()), "+"),
                   ListUtils::concatenate(std::vector<std::string>(expected.begin(), expected.end()), "+"))
      }
      TEST_EQUAL(proteins.at("PEPTIDEB").size(), 1)
    }
    END_SECTION

    START_SECTION([EXTRA] inferPosteriorProbabilities(ConsensusMap&) needs one protein run and probabilities)
    {
      BayesianProteinInferenceAlgorithm bpia;
      setModel(bpia);
      // no identifications
      ConsensusMap empty;
      TEST_EXCEPTION(Exception::MissingInformation, bpia.inferPosteriorProbabilities(empty, false))
      // two runs in two protein runs (merging pools them)
      auto map = nativeMap();
      auto second = nativeMap();
      map.getIdentificationData().merge(second.getIdentificationData());
      TEST_EXCEPTION(Exception::MissingInformation, bpia.inferPosteriorProbabilities(map, false))
      // two legacy protein runs
      auto legacy = nativeMap();
      IdentificationDataConverter::exportConsensusIDs(legacy);
      legacy.getProteinIdentifications().push_back(legacy.getProteinIdentifications()[0]);
      legacy.getProteinIdentifications()[1].setIdentifier("other");
      TEST_EXCEPTION(Exception::MissingInformation, bpia.inferPosteriorProbabilities(legacy, false))
      // no (error) probabilities
      auto scored = nativeMap("XTandem");
      TEST_EXCEPTION(Exception::MissingInformation, bpia.inferPosteriorProbabilities(scored, false))
    }
    END_SECTION

END_TEST
