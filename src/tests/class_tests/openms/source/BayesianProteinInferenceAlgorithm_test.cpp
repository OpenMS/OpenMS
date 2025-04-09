// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;

START_TEST(BayesianProteinInferenceAlgorithm, "$Id$")

    START_SECTION(BayesianProteinInferenceAlgorithm on Protein Peptide ID)
    {
      vector<ProteinIdentification> prots;
      vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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
          vector<PeptideIdentification> peps;
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

END_TEST
