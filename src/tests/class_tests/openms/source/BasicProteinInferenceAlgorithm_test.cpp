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

#include <algorithm>

using namespace OpenMS;
using namespace std;

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

END_TEST
