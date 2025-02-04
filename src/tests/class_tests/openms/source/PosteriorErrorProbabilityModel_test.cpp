// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <vector>
#include <iostream>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

START_TEST(PosteriorErrorProbabilityModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PosteriorErrorProbabilityModel* ptr = nullptr;
PosteriorErrorProbabilityModel* nullPointer = nullptr;
START_SECTION(PosteriorErrorProbabilityModel())
{
	ptr = new PosteriorErrorProbabilityModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~PosteriorErrorProbabilityModel()))
{
  delete ptr;
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void fit( std::vector<double>& search_engine_scores)))
NOT_TESTABLE
//tested below
END_SECTION

START_SECTION((void fit( std::vector<double>& search_engine_scores, std::vector<double>& probabilities)))
	ptr = new PosteriorErrorProbabilityModel();
{
	
	// ------- This code was used for the test file: ------------
	// Use actual Gaussian data to see if fitting works
	//random_device device_random_;
 	//default_random_engine generator_(device_random_());

 	// Gaussian mean and SD, mixture of 2.
 	//normal_distribution<> distribution_1_(1.5, 0.5);
 	//normal_distribution<> distribution_2_(3.5, 1.0);
	// ----------------------------------------------------------

 	vector<double> rand_score_vector;

 	CsvFile gauss_mix (OPENMS_GET_TEST_DATA_PATH("GaussMix_2_1D.csv"), ';');
 	StringList gauss_mix_strings;
 	gauss_mix.getRow(0, gauss_mix_strings);

 	// Load mixture of 2 Gaussians (1D) from provided csv
 	for (StringList::const_iterator it = gauss_mix_strings.begin(); it != gauss_mix_strings.end(); ++it)
 	{
 		if(!it->empty())
 		{
 			rand_score_vector.push_back(it->toDouble());
 		}
 	}

 	TEST_EQUAL(rand_score_vector.size(),2000)

	// Class expects sorted scores
	sort(rand_score_vector.begin(), rand_score_vector.end());
	
	vector<double> probabilities;
	Param param;
	param.setValue("number_of_bins", 10);
	param.setValue("incorrectly_assigned","Gauss");
	ptr->setParameters(param);
	ptr->fit(rand_score_vector, probabilities, "none");
	
	Size i(0),j(1);
	TOLERANCE_ABSOLUTE(0.5)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().x0 , 3.5)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().sigma, 1.0)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().x0, 1.5)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().sigma, 0.5)
	TEST_REAL_SIMILAR(ptr->getNegativePrior(), 0.5)
	TOLERANCE_ABSOLUTE(0.001)
	while(i < rand_score_vector.size() && j < rand_score_vector.size())
	{
		cout<<"i: "<<rand_score_vector[i] << ", j: "<<rand_score_vector[j]<<endl;
		cout<<"pi:"<<probabilities[i] <<", j: "<<probabilities[j]<<endl;
		if(rand_score_vector[i] <= rand_score_vector[j])
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
			TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[i]), probabilities[i])
			TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[j]), probabilities[j])
		}
		else
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
			TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[i]), probabilities[i])
			TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[j]), probabilities[j])
		}
	++i;
	++j;
	}
}

{
	vector<double> score_vector;
	score_vector.push_back(-0.39);
	score_vector.push_back(0.06);
	score_vector.push_back(0.12);
	score_vector.push_back(0.48);
	score_vector.push_back(0.94);
	score_vector.push_back(1.01);
	score_vector.push_back(1.67);
	score_vector.push_back(1.68);
	score_vector.push_back(1.76);
	score_vector.push_back(1.80);
	score_vector.push_back(2.44);
	score_vector.push_back(3.25);
	score_vector.push_back(3.72);
	score_vector.push_back(4.12);
	score_vector.push_back(4.28);
	score_vector.push_back(4.60);
	score_vector.push_back(4.92);
	score_vector.push_back(5.28);
	score_vector.push_back(5.53);
	score_vector.push_back(6.22);
	
	
	vector<double> probabilities;
	Param param;
	param.setValue("number_of_bins", 10);
	param.setValue("incorrectly_assigned","Gumbel");

	ptr->setParameters(param);
	ptr->fit(score_vector, probabilities, "none");
	
	Size i(0),j(1);
	TOLERANCE_ABSOLUTE(0.5)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().x0 , 4.62)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().sigma, 0.87)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().x0, 1.06)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().sigma, 0.77)
	TEST_REAL_SIMILAR(ptr->getNegativePrior(), 0.546)
	TOLERANCE_ABSOLUTE(0.001)
	while(i < score_vector.size() && j < score_vector.size())
	{
		cout<<"i: "<<score_vector[i] << ", j: "<<score_vector[j]<<endl;
		cout<<"pi:"<<probabilities[i] <<", j: "<<probabilities[j]<<endl;
		if(score_vector[i] <= score_vector[j])
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
			TEST_REAL_SIMILAR(ptr->computeProbability(score_vector[i]), probabilities[i])
			TEST_REAL_SIMILAR(ptr->computeProbability(score_vector[j]), probabilities[j])
		}
		else
		{
			TEST_EQUAL(probabilities[i] >= probabilities[j],true)
			TEST_REAL_SIMILAR(ptr->computeProbability(score_vector[i]), probabilities[i])
			TEST_REAL_SIMILAR(ptr->computeProbability(score_vector[j]), probabilities[j])
		}
	++i;
	++j;
	}
}

END_SECTION

START_SECTION((void fillLogDensities(std::vector<double>& x_scores,std::vector<double>& incorrect_density,std::vector<double>& correct_density)))
NOT_TESTABLE
//tested in fit
END_SECTION
START_SECTION((double computeLogLikelihood(std::vector<double>& incorrect_density, std::vector<double>& correct_density)))
NOT_TESTABLE
//tested in fit
END_SECTION
START_SECTION((double getGauss(double x,const GaussFitter::GaussFitResult& params)))
NOT_TESTABLE
//tested in fit
END_SECTION
START_SECTION((double getGumbel(double x,const GaussFitter::GaussFitResult& params)))
NOT_TESTABLE
//tested in fit
END_SECTION

START_SECTION((GaussFitter::GaussFitResult getCorrectlyAssignedFitResult() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((GaussFitter::GaussFitResult getIncorrectlyAssignedFitResult() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((double getNegativePrior() const))
//tested in fit
NOT_TESTABLE
END_SECTION

START_SECTION((double getSmallestScore() const))
TEST_REAL_SIMILAR(ptr->getSmallestScore(), -0.39)
END_SECTION

START_SECTION((const String getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const))
String gumbel = ptr->getGumbelGnuplotFormula(ptr->getIncorrectlyAssignedFitResult());
//approx. f(x) = (1/0.907832") * exp(( 1.48185 - x)/0.907832) * exp(-exp(( 1.48185 - x)/0.907832))"
	cout<<gumbel<<endl;
	TEST_EQUAL(gumbel.hasSubstring("(1/0.90"), true)
	TEST_EQUAL(gumbel.hasSubstring("exp(( 1.47"), true)
	TEST_EQUAL(gumbel.hasSubstring(") * exp(-exp(("), true)
END_SECTION
				
START_SECTION((const String getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const))
String gauss = ptr->getGaussGnuplotFormula(ptr->getCorrectlyAssignedFitResult());
//g(x)=0.444131 * exp(-(x - 5.05539) ** 2 / 2 / (0.898253) ** 2)
	TEST_EQUAL(gauss.hasSubstring(" * exp(-(x - "), true)
	TEST_EQUAL(gauss.hasSubstring(") ** 2 / 2 / ("), true)
	TEST_EQUAL(gauss.hasSubstring(") ** 2)"), true)
END_SECTION

    START_SECTION(fitWithGumbel)
        {
          // ------- This code was used for the test file: ------------
          // Use actual Gaussian data to see if fitting works
          //random_device device_random_;
          //default_random_engine generator_(device_random_());

          // Gaussian mean and SD, mixture of 2.
          //normal_distribution<> distribution_1_(1.5, 0.5);
          //normal_distribution<> distribution_2_(3.5, 1.0);
          // ----------------------------------------------------------

          vector<double> rand_score_vector;

          CsvFile gauss_mix (OPENMS_GET_TEST_DATA_PATH("GumbelGaussMix_2_1D.csv"));
          StringList gauss_mix_strings;
          gauss_mix.getRow(0, gauss_mix_strings);

          // Load mixture of Gumbel and Gaussian (1D) from provided csv
          for (StringList::const_iterator it = gauss_mix_strings.begin(); it != gauss_mix_strings.end(); ++it)
          {
            if(!it->empty())
            {
              rand_score_vector.push_back(it->toDouble());
            }
          }

          TEST_EQUAL(rand_score_vector.size(),2000)

          // Class expects sorted scores
          sort(rand_score_vector.begin(), rand_score_vector.end());

          vector<double> probabilities;
          Param param;
          param.setValue("number_of_bins", 10);
          param.setValue("incorrectly_assigned","Gumbel");
          ptr->setParameters(param);
          ptr->fitGumbelGauss(rand_score_vector, "none");
          double smallest = ptr->getSmallestScore();

          TOLERANCE_ABSOLUTE(0.5)
          TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().x0 , 8.-smallest)
          TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().sigma, 3.5)
          TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedGumbelFitResult().a, 2.-smallest)
          TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedGumbelFitResult().b, .6)
          TEST_REAL_SIMILAR(ptr->getNegativePrior(), 0.6)
          TOLERANCE_ABSOLUTE(0.001)
          /*while(i < rand_score_vector.size() && j < rand_score_vector.size())
          {
            cout<<"i: "<<rand_score_vector[i] << ", j: "<<rand_score_vector[j]<<endl;
            cout<<"pi:"<<probabilities[i] <<", j: "<<probabilities[j]<<endl;
            if(rand_score_vector[i] <= rand_score_vector[j])
            {
              TEST_EQUAL(probabilities[i] >= probabilities[j],true)
              TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[i]), probabilities[i])
              TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[j]), probabilities[j])
            }
            else
            {
              TEST_EQUAL(probabilities[i] >= probabilities[j],true)
              TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[i]), probabilities[i])
              TEST_REAL_SIMILAR(ptr->computeProbability(rand_score_vector[j]), probabilities[j])
            }
            ++i;
            ++j;
          }*/
        }
    END_SECTION

START_SECTION((const String getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect, const GaussFitter::GaussFitResult& correct) const))
NOT_TESTABLE
delete ptr;
END_SECTION

START_SECTION((double computeProbability(double score)))
NOT_TESTABLE
//tested in fit
END_SECTION
START_SECTION((TextFile* InitPlots(std::vector<double> & x_scores)))
NOT_TESTABLE
//tested in fit
END_SECTION
START_SECTION((void	plotTargetDecoyEstimation(std::vector<double> &target,std::vector<double> & decoy)))
NOT_TESTABLE
//not yet tested
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



