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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <algorithm>
#include <cmath>
#include <random>
#include <sstream>
#include <vector>
#include <iostream>
///////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

namespace
{
  // portable pseudo-random numbers (the std::*_distribution classes are implementation-specific)
  struct PortableRNG
  {
    std::mt19937 gen;
    explicit PortableRNG(unsigned seed) : gen(seed) {}
    double uniform() { return (double(gen()) + 0.5) / 4294967296.0; }
    double gumbel(double a, double b) { return a - b * std::log(-std::log(uniform())); }
    double normal(double mu, double sigma)
    {
      const double u1 = uniform(), u2 = uniform();
      return mu + sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * Constants::PI * u2);
    }
  };

  // simulated target-decoy competition: incorrect targets and decoys ~ Gumbel(2, 0.6), correct targets ~ N(5, 1)
  struct TDData
  {
    std::vector<double> targets, decoys;
  };

  TDData simulateTDC(unsigned seed, Size n_incorrect, Size n_correct, Size n_decoys)
  {
    PortableRNG rng(seed);
    TDData data;
    for (Size i = 0; i < n_incorrect; ++i) { data.targets.push_back(rng.gumbel(2.0, 0.6)); }
    for (Size i = 0; i < n_correct; ++i) { data.targets.push_back(rng.normal(5.0, 1.0)); }
    for (Size i = 0; i < n_decoys; ++i) { data.decoys.push_back(rng.gumbel(2.0, 0.6)); }
    return data;
  }

  double gumbelPDF(double x, double a, double b)
  {
    const double z = (x - a) / b;
    return std::exp(-z - std::exp(-z)) / b;
  }

  double gaussPDF(double x, double mu, double sigma)
  {
    const double z = (x - mu) / sigma;
    return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * Constants::PI));
  }

  // true PEP of a target with score s in simulateTDC()
  double truePEP(double s, double n_incorrect, double n_correct)
  {
    const double f0 = n_incorrect * gumbelPDF(s, 2.0, 0.6), f1 = n_correct * gaussPDF(s, 5.0, 1.0);
    return f0 / (f0 + f1);
  }

  PeptideIdentification makeIdentification(const std::string& run, const std::string& score_type, bool higher_better, const std::string& spectrum_reference)
  {
    PeptideIdentification pep;
    pep.setIdentifier(run);
    pep.setScoreType(score_type);
    pep.setHigherScoreBetter(higher_better);
    pep.setSpectrumReference(spectrum_reference);
    return pep;
  }

  PeptideHit makeHit(double score, Int charge, const std::string& target_decoy)
  {
    PeptideHit hit;
    hit.setScore(score);
    hit.setCharge(charge);
    hit.setMetaValue("target_decoy", target_decoy); // "" is used to test robustness against invalid annotations
    return hit;
  }
}

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
 			rand_score_vector.push_back(StringUtils::toDouble(*it));
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
	// The incorrect component is a Gumbel distribution fitted by maximum likelihood (x0/sigma = location/scale).
	// (Up to OpenMS 3.5, a Gaussian was fitted and its mean/sd were plugged into a Gumbel density for the probabilities.)
	// Reference: an independent EM implementation converges to the same (non-degenerate) maximum:
	// location 1.17, scale 0.97, mean 5.16, sd 0.85, prior 0.61 (scores shifted by |smallest| + 0.001).
	TOLERANCE_ABSOLUTE(0.05)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().x0 , 5.16)
	TEST_REAL_SIMILAR(ptr->getCorrectlyAssignedFitResult().sigma, 0.85)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().x0, 1.17)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedFitResult().sigma, 0.97)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedGumbelFitResult().a, 1.17)
	TEST_REAL_SIMILAR(ptr->getIncorrectlyAssignedGumbelFitResult().b, 0.97)
	TEST_REAL_SIMILAR(ptr->getNegativePrior(), 0.61)
	TEST_EQUAL(ptr->getModelType() == PosteriorErrorProbabilityModel::ModelType::MIXTURE, true)
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

START_SECTION((const std::string getGumbelGnuplotFormula(const GaussFitter::GaussFitResult& params) const))
std::string gumbel = ptr->getGumbelGnuplotFormula(ptr->getIncorrectlyAssignedFitResult());
//approx. f(x) = (1/0.965686) * exp(( 1.17282- x)/0.965686) * exp(-exp((1.17282 - x)/0.965686))
	cout<<gumbel<<endl;
	std::ostringstream scale, location;
	scale << "(1/" << ptr->getIncorrectlyAssignedFitResult().sigma << ")";
	location << "exp(( " << ptr->getIncorrectlyAssignedFitResult().x0 << "- x)";
	TEST_EQUAL(StringUtils::hasSubstring(gumbel, scale.str()), true)
	TEST_EQUAL(StringUtils::hasSubstring(gumbel, location.str()), true)
	TEST_EQUAL(StringUtils::hasSubstring(gumbel, "(1/0.96"), true)
	TEST_EQUAL(StringUtils::hasSubstring(gumbel, "exp(( 1.17"), true)
	TEST_EQUAL(StringUtils::hasSubstring(gumbel, ") * exp(-exp(("), true)
END_SECTION
				
START_SECTION((const std::string getGaussGnuplotFormula(const GaussFitter::GaussFitResult& params) const))
std::string gauss = ptr->getGaussGnuplotFormula(ptr->getCorrectlyAssignedFitResult());
//g(x)=0.444131 * exp(-(x - 5.05539) ** 2 / 2 / (0.898253) ** 2)
	TEST_EQUAL(StringUtils::hasSubstring(gauss, " * exp(-(x - "), true)
	TEST_EQUAL(StringUtils::hasSubstring(gauss, ") ** 2 / 2 / ("), true)
	TEST_EQUAL(StringUtils::hasSubstring(gauss, ") ** 2)"), true)
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
              rand_score_vector.push_back(StringUtils::toDouble(*it));
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

START_SECTION((const std::string getBothGnuplotFormula(const GaussFitter::GaussFitResult& incorrect, const GaussFitter::GaussFitResult& correct) const))
NOT_TESTABLE
delete ptr;
END_SECTION

START_SECTION((double computeProbability(double score)))
{
  // The probabilities must be computed from the very same densities that were fitted
  // (up to OpenMS 3.5, a Gaussian fit was plugged into a Gumbel density).
  vector<double> scores;
  CsvFile gumbel_gauss_mix(OPENMS_GET_TEST_DATA_PATH("GumbelGaussMix_2_1D.csv"));
  StringList strings;
  gumbel_gauss_mix.getRow(0, strings);
  for (const auto& s : strings)
  {
    if (!s.empty()) scores.push_back(StringUtils::toDouble(s));
  }
  for (const std::string incorrect : {"Gumbel", "Gauss"})
  {
    PosteriorErrorProbabilityModel model;
    Param p = model.getParameters();
    p.setValue("incorrectly_assigned", incorrect);
    model.setParameters(p);
    vector<double> fit_scores(scores);
    TEST_EQUAL(model.fit(fit_scores, "none"), true)
    const double shift = std::fabs(model.getSmallestScore()) + 0.001;
    const double prior = model.getNegativePrior();
    const auto inc = model.getIncorrectlyAssignedFitResult();
    const auto cor = model.getCorrectlyAssignedFitResult();
    TEST_EQUAL(inc.x0 < cor.x0, true)
    // scores between the mode of the incorrect and the mean of the correct component (no tail clamping)
    for (double s = inc.x0 - shift + 0.5; s < cor.x0 - shift; s += 1.0)
    {
      const double x = s + shift;
      const double f0 = (incorrect == "Gumbel") ? gumbelPDF(x, inc.x0, inc.sigma) : gaussPDF(x, inc.x0, inc.sigma);
      const double f1 = gaussPDF(x, cor.x0, cor.sigma);
      TOLERANCE_ABSOLUTE(1e-9)
      TEST_REAL_SIMILAR(model.computeProbability(s), prior * f0 / (prior * f0 + (1 - prior) * f1))
      vector<double> inc_density, cor_density;
      model.fillDensities({x}, inc_density, cor_density);
      TEST_REAL_SIMILAR(inc_density[0], f0)
      TEST_REAL_SIMILAR(cor_density[0], f1)
    }
    if (incorrect == "Gumbel")
    {
      // the generating parameters (Gumbel(2, 0.6), N(8, 3.5), prior 0.6) are recovered
      TOLERANCE_ABSOLUTE(0.15)
      TEST_REAL_SIMILAR(inc.x0 - shift, 2.0)
      TEST_REAL_SIMILAR(inc.sigma, 0.6)
      TEST_REAL_SIMILAR(model.getIncorrectlyAssignedGumbelFitResult().a, inc.x0)
      TEST_REAL_SIMILAR(model.getNegativePrior(), 0.6)
      TOLERANCE_ABSOLUTE(0.3)
      TEST_REAL_SIMILAR(cor.x0 - shift, 8.0)
      TEST_REAL_SIMILAR(cor.sigma, 3.5)
    }
    // probabilities do not increase with the score, also in the tails
    double previous = 1.0;
    bool monotone = true;
    for (double s = -20.0; s < 60.0; s += 0.05)
    {
      const double pep = model.computeProbability(s);
      monotone = monotone && pep <= previous + 1e-12 && pep >= 0.0 && pep <= 1.0;
      previous = pep;
    }
    TEST_EQUAL(monotone, true)
  }

  // non-finite scores are ignored for fitting
  {
    PosteriorErrorProbabilityModel model;
    vector<double> fit_scores(scores);
    fit_scores.push_back(std::numeric_limits<double>::infinity());
    fit_scores.push_back(-std::numeric_limits<double>::infinity());
    fit_scores.push_back(std::numeric_limits<double>::quiet_NaN());
    PosteriorErrorProbabilityModel reference;
    vector<double> reference_scores(scores);
    TEST_EQUAL(model.fit(fit_scores, "none"), true)
    TEST_EQUAL(reference.fit(reference_scores, "none"), true)
    TOLERANCE_ABSOLUTE(1e-12)
    TEST_REAL_SIMILAR(model.getNegativePrior(), reference.getNegativePrior())
    TEST_REAL_SIMILAR(model.computeProbability(5.0), reference.computeProbability(5.0))
    TEST_REAL_SIMILAR(model.computeProbability(std::numeric_limits<double>::infinity()), 0.0)
    TEST_REAL_SIMILAR(model.computeProbability(-std::numeric_limits<double>::infinity()), 1.0)
  }
}
END_SECTION

START_SECTION((bool fitTargetDecoy(const std::vector<double>& target_scores, const std::vector<double>& decoy_scores, const std::string& outlier_handling)))
{
  const Size n0 = 4000, n1 = 2000;
  // target-decoy competition: one decoy per incorrect target
  TDData data = simulateTDC(42, n0, n1, n0);

  // non-parametric estimate (default with enough decoys)
  PosteriorErrorProbabilityModel model;
  TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::NONE, true)
  TEST_EQUAL(model.fitTargetDecoy(data.targets, data.decoys, "ignore_iqr_outliers"), true)
  TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC, true)
  TEST_REAL_SIMILAR(model.getTargetDecoyIncorrectFraction(), double(n0) / double(n0 + n1))
  TEST_EQUAL(model.hasCalibrationWarning(), false)

  double previous = 1.0;
  bool monotone = true;
  for (double s = -5.0; s < 15.0; s += 0.01)
  {
    const double pep = model.computeProbability(s);
    monotone = monotone && pep <= previous + 1e-12 && pep >= 0.0 && pep <= 1.0;
    previous = pep;
  }
  TEST_EQUAL(monotone, true)

  // calibration: pointwise (the step function of the isotonic estimate is coarsest where the true PEP is steepest)
  TOLERANCE_ABSOLUTE(0.1)
  for (double s : {2.5, 3.0, 3.5, 4.0, 4.5, 5.0})
  {
    TEST_REAL_SIMILAR(model.computeProbability(s), truePEP(s, n0, n1))
  }
  // ... and in aggregate: the PEPs sum up to the number of incorrect targets
  double sum_pep = 0.0;
  for (double t : data.targets) { sum_pep += model.computeProbability(t); }
  TOLERANCE_ABSOLUTE(0.05)
  TEST_REAL_SIMILAR(sum_pep / n0, 1.0)
  // targets above all decoys get a small, but positive PEP: at most the Jeffreys estimate for the targets above the best decoy
  {
    const double best_decoy = *std::max_element(data.decoys.begin(), data.decoys.end());
    const double n_above = double(std::count_if(data.targets.begin(), data.targets.end(), [best_decoy](double t) { return t > best_decoy; }));
    const double p = 0.5 / (n_above + 1.0);
    TEST_EQUAL(model.computeProbability(1e6) > 0.0, true)
    TEST_EQUAL(model.computeProbability(1e6) <= p / (1.0 - p) + 1e-12, true)
  }

  // the estimate only depends on the ordering of the scores: invariant to monotone transformations
  {
    TDData transformed = data;
    for (double& t : transformed.targets) { t = std::exp(t / 2.0); }
    for (double& d : transformed.decoys) { d = std::exp(d / 2.0); }
    PosteriorErrorProbabilityModel model2;
    TEST_EQUAL(model2.fitTargetDecoy(transformed.targets, transformed.decoys, "none"), true)
    auto curve1 = model.getIsotonicPEPCurve(), curve2 = model2.getIsotonicPEPCurve();
    TEST_EQUAL(curve1.size(), curve2.size())
    double max_difference = 0.0;
    for (Size i = 0; i < std::min(curve1.size(), curve2.size()); ++i)
    {
      max_difference = std::max(max_difference, std::fabs(curve1[i].second - curve2[i].second));
    }
    TOLERANCE_ABSOLUTE(1e-12)
    TEST_REAL_SIMILAR(max_difference, 0.0)
  }

  // discrete scores: tied scores are pooled, the result does not depend on the order of the input
  {
    TDData discrete = data;
    for (double& t : discrete.targets) { t = std::round(t * 2.0) / 2.0; }
    for (double& d : discrete.decoys) { d = std::round(d * 2.0) / 2.0; }
    PosteriorErrorProbabilityModel model3, model4;
    TEST_EQUAL(model3.fitTargetDecoy(discrete.targets, discrete.decoys, "none"), true)
    std::reverse(discrete.targets.begin(), discrete.targets.end());
    std::reverse(discrete.decoys.begin(), discrete.decoys.end());
    TEST_EQUAL(model4.fitTargetDecoy(discrete.targets, discrete.decoys, "none"), true)
    auto curve3 = model3.getIsotonicPEPCurve(), curve4 = model4.getIsotonicPEPCurve();
    TEST_EQUAL(curve3.size(), curve4.size())
    bool identical = curve3.size() == curve4.size();
    for (Size i = 0; identical && i < curve3.size(); ++i)
    {
      identical = curve3[i] == curve4[i];
    }
    TEST_EQUAL(identical, true)
    TOLERANCE_ABSOLUTE(0.1)
    TEST_REAL_SIMILAR(model3.computeProbability(3.5), truePEP(3.5, n0, n1))
  }

  // several decoys per target: the number of incorrect targets is #decoys / decoy_ratio
  {
    PosteriorErrorProbabilityModel model5;
    Param p = model5.getParameters();
    p.setValue("decoy_ratio", 2.0);
    model5.setParameters(p);
    TEST_EQUAL(model5.fitTargetDecoy(data.targets, data.decoys, "none"), true)
    auto curve1 = model.getIsotonicPEPCurve(), curve5 = model5.getIsotonicPEPCurve();
    TEST_EQUAL(curve1.size(), curve5.size())
    TOLERANCE_ABSOLUTE(1e-12)
    for (Size i = 0; i < std::min(curve1.size(), curve5.size()); ++i)
    {
      if (curve1[i].second < 1.0)
      {
        TEST_REAL_SIMILAR(curve5[i].second, curve1[i].second / 2.0)
      }
    }
    TEST_REAL_SIMILAR(model5.getTargetDecoyIncorrectFraction(), double(n0) / double(n0 + n1) / 2.0)
  }

  // semi-supervised mixture model
  {
    PosteriorErrorProbabilityModel model6;
    Param p = model6.getParameters();
    p.setValue("decoy_handling", "semi_supervised");
    model6.setParameters(p);
    TEST_EQUAL(model6.fitTargetDecoy(data.targets, data.decoys, "none"), true)
    TEST_EQUAL(model6.getModelType() == PosteriorErrorProbabilityModel::ModelType::SEMI_SUPERVISED_MIXTURE, true)
    TOLERANCE_ABSOLUTE(0.03)
    TEST_REAL_SIMILAR(model6.getNegativePrior(), double(n0) / double(n0 + n1))
    TEST_REAL_SIMILAR(model6.getIncorrectlyAssignedGumbelFitResult().b, 0.6)
    TOLERANCE_ABSOLUTE(0.1)
    TEST_REAL_SIMILAR(model6.getCorrectlyAssignedFitResult().sigma, 1.0)
    TOLERANCE_ABSOLUTE(0.05)
    for (double s : {2.5, 3.0, 3.5, 4.0, 4.5, 5.0})
    {
      TEST_REAL_SIMILAR(model6.computeProbability(s), truePEP(s, n0, n1))
    }
    TEST_EQUAL(model6.hasCalibrationWarning(), false)
  }

  // 'auto': the estimator depends on the number of decoys
  {
    const std::vector<std::pair<Size, PosteriorErrorProbabilityModel::ModelType>> expected =
    {
      {5, PosteriorErrorProbabilityModel::ModelType::MIXTURE},
      {20, PosteriorErrorProbabilityModel::ModelType::SEMI_SUPERVISED_MIXTURE},
      {60, PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC}
    };
    for (const auto& [n_decoys, type] : expected)
    {
      TDData small = simulateTDC(7, 400, 200, n_decoys);
      PosteriorErrorProbabilityModel model7;
      TEST_EQUAL(model7.fitTargetDecoy(small.targets, small.decoys, "none"), true)
      TEST_EQUAL(model7.getModelType() == type, true)
    }
    // 'ignore': unsupervised model of targets and decoys
    PosteriorErrorProbabilityModel model8;
    Param p = model8.getParameters();
    p.setValue("decoy_handling", "ignore");
    model8.setParameters(p);
    TEST_EQUAL(model8.fitTargetDecoy(data.targets, data.decoys, "none"), true)
    TEST_EQUAL(model8.getModelType() == PosteriorErrorProbabilityModel::ModelType::MIXTURE, true)
    // without targets, nothing can be fitted
    PosteriorErrorProbabilityModel model9;
    TEST_EQUAL(model9.fitTargetDecoy({}, data.decoys, "none"), false)
    TEST_EQUAL(model9.getModelType() == PosteriorErrorProbabilityModel::ModelType::NONE, true)
  }
}
END_SECTION

START_SECTION((static std::map<std::string, std::vector<std::vector<double>>> extractAndTransformScores(const std::vector<ProteinIdentification> & protein_ids, const PeptideIdentificationList & peptide_ids, const bool split_charge, const bool top_hits_only)))
{
  std::vector<ProteinIdentification> proteins(4);
  proteins[0].setIdentifier("xtandem");
  proteins[0].setSearchEngine("X! Tandem"); // pepXML engine name
  proteins[1].setIdentifier("qvalues");
  proteins[1].setSearchEngine("SomeRescorer"); // main score is a q-value: unusable
  proteins[2].setIdentifier("generic");
  proteins[2].setSearchEngine("SomeEngine"); // unknown engine with an e-value
  proteins[3].setIdentifier("raw");
  proteins[3].setSearchEngine("OtherEngine"); // unknown engine with a raw score

  PeptideIdentificationList peptides;
  // X! Tandem: main score hyperscore, e-value in meta value 'expect'
  const std::vector<std::pair<double, std::string>> xtandem = {{1e-6, "target"}, {1e-3, "decoy"}, {1e-2, "target+decoy"}, {0.5, ""}, {2.0, "decoy"}, {0.0, "target"}};
  for (Size i = 0; i < xtandem.size(); ++i)
  {
    PeptideIdentification pep = makeIdentification("xtandem", "XTandem", true, "scan=" + StringUtils::toStr(i));
    PeptideHit hit = makeHit(10.0 - double(i), 2, xtandem[i].second);
    hit.setMetaValue("expect", xtandem[i].first);
    PeptideHit second = makeHit(-1.0, 2, "target"); // rank 2
    second.setMetaValue("expect", 50.0);
    pep.setHits({hit, second});
    peptides.push_back(pep);
  }
  // separate decoy search of spectrum "scan=0": competition with the target top hit above (decoy scores worse)
  {
    PeptideIdentification pep = makeIdentification("xtandem", "XTandem", true, "scan=0");
    PeptideHit hit = makeHit(1.0, 2, "decoy");
    hit.setMetaValue("expect", 1e-1);
    pep.setHits({hit});
    peptides.push_back(pep);
  }
  for (Size i = 0; i < 4; ++i)
  {
    PeptideIdentification pep = makeIdentification("qvalues", "q-value", false, "");
    pep.setHits({makeHit(0.01 * double(i), 2, "target")});
    peptides.push_back(pep);
  }
  for (Size i = 0; i < 4; ++i)
  {
    PeptideIdentification pep = makeIdentification("generic", "E-value", false, "");
    pep.setHits({makeHit(std::pow(10.0, -double(i)), 3, i == 3 ? "decoy" : "target")});
    peptides.push_back(pep);
  }
  for (Size i = 0; i < 4; ++i)
  {
    PeptideIdentification pep = makeIdentification("raw", "my_distance", false, "");
    pep.setHits({makeHit(double(i), 1, "target")});
    peptides.push_back(pep);
  }

  auto all_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(proteins, peptides, false, false);
  TEST_EQUAL(all_scores.size(), 3)
  TEST_EQUAL(all_scores.count("XTANDEM"), 1)
  TEST_EQUAL(all_scores.count("SOMEENGINE"), 1)
  TEST_EQUAL(all_scores.count("OTHERENGINE"), 1)
  TEST_EQUAL(all_scores.count("SOMERESCORER"), 0)

  const auto& xt = all_scores["XTANDEM"];
  // all hits of all identifications (7 identifications, 6 with two hits)
  TEST_EQUAL(xt[0].size(), 13)
  // top hits: targets (incl. target+decoy; unknown label excluded), decoys; spectrum scan=0 is identified twice, the target wins
  TEST_EQUAL(xt[1].size(), 3)
  TEST_EQUAL(xt[2].size(), 2)
  TOLERANCE_ABSOLUTE(1e-9)
  TEST_REAL_SIMILAR(*std::max_element(xt[1].begin(), xt[1].end()), -std::log10(std::numeric_limits<double>::min())) // e-value 0
  TEST_REAL_SIMILAR(*std::min_element(xt[1].begin(), xt[1].end()), 2.0) // 1e-2
  TEST_REAL_SIMILAR(*std::max_element(xt[2].begin(), xt[2].end()), 3.0) // 1e-3

  // generic e-value: -log10
  const auto& generic = all_scores["SOMEENGINE"];
  TEST_EQUAL(generic[0].size(), 4)
  TEST_REAL_SIMILAR(*std::max_element(generic[0].begin(), generic[0].end()), 3.0)
  // generic lower-is-better raw score: negated
  const auto& raw = all_scores["OTHERENGINE"];
  TEST_REAL_SIMILAR(*std::max_element(raw[0].begin(), raw[0].end()), 0.0)
  TEST_REAL_SIMILAR(*std::min_element(raw[0].begin(), raw[0].end()), -3.0)

  // split by charge
  auto by_charge = PosteriorErrorProbabilityModel::extractAndTransformScores(proteins, peptides, true, true);
  TEST_EQUAL(by_charge.count("XTANDEM,2"), 1)
  TEST_EQUAL(by_charge.count("SOMEENGINE,3"), 1)
  TEST_EQUAL(by_charge.count("OTHERENGINE,1"), 1)
  TEST_EQUAL(by_charge["XTANDEM,2"][0].size(), 7) // top hits only

  // deprecated overload: no target/decoy scores if not requested
  auto deprecated = PosteriorErrorProbabilityModel::extractAndTransformScores(proteins, peptides, false, false, false, 0.05);
  TEST_EQUAL(deprecated["XTANDEM"][1].empty() && deprecated["XTANDEM"][2].empty(), true)

  // search engine keys
  TEST_EQUAL(PosteriorErrorProbabilityModel::getSearchEngineKey(proteins[0]), "XTANDEM")
  ProteinIdentification msgf;
  msgf.setSearchEngine("MS-GF+");
  TEST_EQUAL(PosteriorErrorProbabilityModel::getSearchEngineKey(msgf), "MSGFPLUS")
  msgf.setSearchEngine("MSGFPlus");
  TEST_EQUAL(PosteriorErrorProbabilityModel::getSearchEngineKey(msgf), "MSGFPLUS")

  // update scores: only identifications of the fitted search engine are updated
  PosteriorErrorProbabilityModel model;
  TDData data = simulateTDC(3, 400, 200, 400);
  TEST_EQUAL(model.fit(data.targets, "none"), true)
  bool unable_to_fit = true, not_well_fit = true;
  std::set<Size> updated;
  PosteriorErrorProbabilityModel::updateScores(model, "XTANDEM", -1, false, false, proteins, peptides, unable_to_fit, not_well_fit, updated);
  TEST_EQUAL(unable_to_fit, false)
  TEST_EQUAL(updated.size(), 7)
  bool only_xtandem = true, valid = true;
  for (Size index : updated)
  {
    only_xtandem = only_xtandem && peptides[index].getIdentifier() == "xtandem";
    for (const PeptideHit& hit : peptides[index].getHits())
    {
      valid = valid && hit.getScore() >= 0.0 && hit.getScore() <= 1.0 && hit.metaValueExists("XTandem_score");
    }
  }
  TEST_EQUAL(only_xtandem, true)
  TEST_EQUAL(valid, true)
}
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



