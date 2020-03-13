// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
          bpia.inferPosteriorProbabilities(prots,peps);
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
