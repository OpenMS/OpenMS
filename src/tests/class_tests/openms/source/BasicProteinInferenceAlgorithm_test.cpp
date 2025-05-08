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

using namespace OpenMS;
using namespace std;

START_TEST(BasicProteinInferenceAlgorithm, "$Id$")

    START_SECTION(BasicProteinInferenceAlgorithm on Protein Peptide ID)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
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
      vector<PeptideIdentification> peps;
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
      vector<PeptideIdentification> peps;
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
      vector<PeptideIdentification> peps;
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
      vector<PeptideIdentification> peps;
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
      vector<PeptideIdentification> peps;
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

END_TEST
