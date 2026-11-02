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
  // targets above all decoys get small, but positive PEPs, which sum up to the expected number of incorrect targets above the best decoy (1)
  {
    const double best_decoy = *std::max_element(data.decoys.begin(), data.decoys.end());
    double sum_above = 0.0;
    for (double t : data.targets)
    {
      if (t > best_decoy) { sum_above += model.computeProbability(t); }
    }
    TEST_EQUAL(model.computeProbability(1e6) > 0.0, true)
    TOLERANCE_ABSOLUTE(1e-9)
    TEST_REAL_SIMILAR(sum_above, 1.0)
  }
  // the "+1" of the targets above the best decoy is pooled with the following block if necessary (no inflation, no ties):
  // best first: 9 targets, decoy, 19 targets, decoy, then targets and decoys alternating
  {
    std::vector<double> targets, decoys;
    double s = 200.0;
    for (int i = 0; i < 9; ++i) { targets.push_back(s--); }
    decoys.push_back(s--);
    for (int i = 0; i < 19; ++i) { targets.push_back(s--); }
    decoys.push_back(s--);
    for (int i = 0; i < 40; ++i) { targets.push_back(s--); decoys.push_back(s--); }
    PosteriorErrorProbabilityModel pooled;
    Param p = pooled.getParameters();
    p.setValue("decoy_handling", "isotonic");
    pooled.setParameters(p);
    TEST_EQUAL(pooled.fitTargetDecoy(targets, decoys, "none"), true)
    double sum_top = 0.0, previous = 0.0;
    bool strictly_decreasing = true;
    for (int i = 27; i >= 0; --i) // the 28 targets above the second decoy, from worst to best
    {
      const double pep = pooled.computeProbability(targets[i]);
      strictly_decreasing = strictly_decreasing && (i == 27 || pep < previous);
      previous = pep;
      sum_top += pep;
    }
    // the 9 targets above the best decoy are pooled with the next 19 targets and the decoy (PEP at the mean score 186: 2 / 28)
    TOLERANCE_ABSOLUTE(1e-9)
    TEST_REAL_SIMILAR(pooled.computeProbability(186.0), 2.0 / 28.0)
    // target-decoy estimate: (#decoys + 1) = 2 incorrect targets (up to the small shift of the linear ramp within the block)
    TOLERANCE_ABSOLUTE(0.02)
    TEST_REAL_SIMILAR(sum_top, 2.0)
    TEST_EQUAL(strictly_decreasing, true)
  }

  // PSMs are ordered by score (no ties of different scores), except for PEPs of 1
  {
    std::vector<double> sorted_targets(data.targets);
    std::sort(sorted_targets.begin(), sorted_targets.end());
    bool strictly_decreasing = true;
    for (Size i = 1; i < sorted_targets.size(); ++i)
    {
      const double previous = model.computeProbability(sorted_targets[i - 1]), current = model.computeProbability(sorted_targets[i]);
      if (sorted_targets[i] > sorted_targets[i - 1] && previous < 1.0)
      {
        strictly_decreasing = strictly_decreasing && current < previous;
      }
    }
    TEST_EQUAL(strictly_decreasing, true)
  }

  // the estimate only depends on the ordering of the scores: invariant to monotone transformations
  {
    TDData transformed = data;
    for (double& t : transformed.targets) { t = std::exp(t / 2.0); }
    for (double& d : transformed.decoys) { d = std::exp(d / 2.0); }
    PosteriorErrorProbabilityModel model2;
    TEST_EQUAL(model2.fitTargetDecoy(transformed.targets, transformed.decoys, "none"), true)
    // same blocks (the linear ordering of the PSMs within a block depends on the score scale, but not the blocks themselves) ...
    TEST_EQUAL(model.getIsotonicPEPCurve().size(), model2.getIsotonicPEPCurve().size())
    // ... hence the same total PEP (each block keeps its mean PEP) and the same ordering of the PSMs by PEP
    std::vector<double> all(data.targets), all_transformed(transformed.targets);
    all.insert(all.end(), data.decoys.begin(), data.decoys.end());
    all_transformed.insert(all_transformed.end(), transformed.decoys.begin(), transformed.decoys.end());
    double sum1 = 0.0, sum2 = 0.0;
    std::vector<std::pair<double, Size>> order1, order2;
    for (Size i = 0; i < all.size(); ++i)
    {
      const double pep1 = model.computeProbability(all[i]), pep2 = model2.computeProbability(all_transformed[i]);
      sum1 += pep1;
      sum2 += pep2;
      order1.emplace_back(pep1, i);
      order2.emplace_back(pep2, i);
    }
    TOLERANCE_ABSOLUTE(1e-6)
    TEST_REAL_SIMILAR(sum1, sum2)
    std::sort(order1.begin(), order1.end());
    std::sort(order2.begin(), order2.end());
    bool same_order = true;
    for (Size i = 0; i < order1.size(); ++i)
    {
      same_order = same_order && order1[i].second == order2[i].second;
    }
    TEST_EQUAL(same_order, true)
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

  // two decoys per target (each incorrect spectrum is won by a decoy with probability 2/3): the number of incorrect targets is #decoys / decoy_ratio
  {
    TDData two_decoys = simulateTDC(43, n0, n1, 2 * n0);
    PosteriorErrorProbabilityModel model5;
    Param p = model5.getParameters();
    p.setValue("decoy_ratio", 2.0);
    model5.setParameters(p);
    TEST_EQUAL(model5.fitTargetDecoy(two_decoys.targets, two_decoys.decoys, "none"), true)
    TOLERANCE_ABSOLUTE(1e-9)
    TEST_REAL_SIMILAR(model5.getTargetDecoyIncorrectFraction(), double(n0) / double(n0 + n1))
    double sum_pep5 = 0.0, sum_above5 = 0.0;
    const double best_decoy = *std::max_element(two_decoys.decoys.begin(), two_decoys.decoys.end());
    for (double t : two_decoys.targets)
    {
      sum_pep5 += model5.computeProbability(t);
      if (t > best_decoy) { sum_above5 += model5.computeProbability(t); }
    }
    TEST_REAL_SIMILAR(sum_above5, 0.5) // 1 / decoy_ratio
    TOLERANCE_ABSOLUTE(0.05)
    TEST_REAL_SIMILAR(sum_pep5 / n0, 1.0)
    TEST_EQUAL(model5.hasCalibrationWarning(), false)
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
      {0, PosteriorErrorProbabilityModel::ModelType::MIXTURE},
      {5, PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC},
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

  // calibration check: with a wrong decoy_ratio, there are 4 times more decoys than the (correct) PEPs of the semi-supervised model imply
  {
    PosteriorErrorProbabilityModel model10;
    Param p = model10.getParameters();
    p.setValue("decoy_handling", "semi_supervised");
    p.setValue("decoy_ratio", 0.25);
    model10.setParameters(p);
    TEST_EQUAL(model10.fitTargetDecoy(data.targets, data.decoys, "none"), true)
    TEST_EQUAL(model10.hasCalibrationWarning(), true)
  }
}
END_SECTION

START_SECTION(([EXTRA] degenerate input))
{
  auto check_monotone = [](const PosteriorErrorProbabilityModel& model)
  {
    double previous = 1.0;
    bool ok = true;
    for (double s = -20.0; s < 20.0; s += 0.01)
    {
      const double pep = model.computeProbability(s);
      ok = ok && pep <= previous + 1e-12 && pep >= 0.0 && pep <= 1.0;
      previous = pep;
    }
    return ok;
  };
  const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();

  // NaN values are sorted to the end
  {
    PosteriorErrorProbabilityModel model;
    std::vector<double> scores = {3.0, nan, 1.0, 2.0, nan, 0.5};
    model.fit(scores, "none");
    TEST_EQUAL(std::is_sorted(scores.begin(), scores.begin() + 4), true)
    TEST_EQUAL(std::isnan(scores[4]) && std::isnan(scores[5]), true)
  }

  // unsupervised: too few scores, no spread, non-finite scores only
  for (std::vector<double> scores : std::vector<std::vector<double>>{{}, {1.0}, {1.0, 2.0}, std::vector<double>(1000, 4.2), {nan, inf, -inf}})
  {
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fit(scores, "none"), false)
    TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::NONE, true)
  }

  // target-decoy: all scores identical -> no information, PEP 1
  {
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy(std::vector<double>(500, 3.0), std::vector<double>(500, 3.0), "none"), true)
    TEST_REAL_SIMILAR(model.computeProbability(3.0), 1.0)
    TEST_REAL_SIMILAR(model.computeProbability(100.0), 1.0)
  }

  TDData null_data = simulateTDC(11, 2000, 0, 2000);
  // decoys score better than targets (e.g. wrong score orientation): PEP 1 for all targets
  {
    std::vector<double> better_decoys;
    for (double t : null_data.targets) { better_decoys.push_back(t + 10.0); }
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy(null_data.decoys, better_decoys, "none"), true)
    double max_pep = 0.0, min_pep = 1.0;
    for (double t : null_data.decoys)
    {
      min_pep = std::min(min_pep, model.computeProbability(t));
      max_pep = std::max(max_pep, model.computeProbability(t));
    }
    TEST_REAL_SIMILAR(min_pep, 1.0)
    TEST_EQUAL(check_monotone(model), true)
  }

  // no correct targets at all: the default (non-parametric) estimate assigns no target a low PEP
  {
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy(null_data.targets, null_data.decoys, "none"), true)
    TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC, true)
    Size low_pep = std::count_if(null_data.targets.begin(), null_data.targets.end(), [&model](double t) { return model.computeProbability(t) < 0.05; });
    TEST_EQUAL(low_pep, 0)
    TEST_EQUAL(check_monotone(model), true)
  }

  // perfect separation: small (but positive) PEPs for targets above all decoys, PEP 1 below
  {
    TDData separated = simulateTDC(12, 0, 500, 500);
    for (double& t : separated.targets) { t += 10.0; }
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy(separated.targets, separated.decoys, "none"), true)
    const double best_decoy = *std::max_element(separated.decoys.begin(), separated.decoys.end());
    TEST_EQUAL(model.computeProbability(best_decoy + 5.0) > 0.0, true)
    TEST_EQUAL(model.computeProbability(best_decoy + 5.0) < 0.01, true)
    TEST_REAL_SIMILAR(model.computeProbability(best_decoy), 1.0)
    TEST_EQUAL(check_monotone(model), true)
  }

  // a single target: no evidence that it is correct (above all decoys, it is incorrect with probability 1/2 under target-decoy competition)
  {
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy({100.0}, null_data.decoys, "none"), true)
    TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC, true)
    TEST_REAL_SIMILAR(model.computeProbability(100.0), 1.0)
    TEST_EQUAL(model.fitTargetDecoy({5.0}, {3.0}, "none"), true) // one target, one decoy
    TEST_REAL_SIMILAR(model.computeProbability(5.0), 1.0)
    TEST_REAL_SIMILAR(model.computeProbability(3.0), 1.0)
    Param p = model.getParameters();
    p.setValue("decoy_ratio", 4.0); // 4 decoys per target: an incorrect PSM scores above the best of the decoys with probability 1/5
    model.setParameters(p);
    TEST_EQUAL(model.fitTargetDecoy({100.0}, null_data.decoys, "none"), true)
    TEST_REAL_SIMILAR(model.computeProbability(100.0), 0.25)
  }

  // a single decoy (among many correct targets): still informative, about one incorrect target is expected above it
  {
    TDData one_decoy = simulateTDC(13, 0, 1000, 1);
    one_decoy.decoys[0] = 3.5;
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy(one_decoy.targets, one_decoy.decoys, "none"), true)
    TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC, true)
    double sum_pep = 0.0, sum_above = 0.0;
    for (double t : one_decoy.targets)
    {
      sum_pep += model.computeProbability(t);
      if (t > 3.5) { sum_above += model.computeProbability(t); }
    }
    TOLERANCE_ABSOLUTE(1e-9)
    TEST_REAL_SIMILAR(sum_above, 1.0)
    TEST_EQUAL(sum_pep < 3.0, true) // an unsupervised mixture model would assign much larger PEPs to (many of) these correct targets
    TEST_EQUAL(check_monotone(model), true)
    // forced semi-supervised mixture with too few decoys: non-parametric estimate
    Param p = model.getParameters();
    p.setValue("decoy_handling", "semi_supervised");
    model.setParameters(p);
    TEST_EQUAL(model.fitTargetDecoy(one_decoy.targets, one_decoy.decoys, "none"), true)
    TEST_EQUAL(model.getModelType() == PosteriorErrorProbabilityModel::ModelType::TARGET_DECOY_ISOTONIC, true)
  }

  // non-finite scores are ignored, no (finite) target is not enough
  {
    PosteriorErrorProbabilityModel model;
    TEST_EQUAL(model.fitTargetDecoy({nan, inf}, null_data.decoys, "none"), false)
    std::vector<double> targets(null_data.targets);
    targets.insert(targets.end(), {nan, inf, -inf});
    TEST_EQUAL(model.fitTargetDecoy(targets, null_data.decoys, "none"), true)
    TEST_EQUAL(check_monotone(model), true)
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
  std::set<std::pair<Size, Size>> updated;
  PosteriorErrorProbabilityModel::updateScores(model, "XTANDEM", -1, false, false, proteins, peptides, unable_to_fit, not_well_fit, updated);
  TEST_EQUAL(unable_to_fit, false)
  TEST_EQUAL(updated.size(), 13) // all hits of the X! Tandem run
  bool only_xtandem = true, valid = true;
  for (const auto& [pep_index, hit_index] : updated)
  {
    only_xtandem = only_xtandem && peptides[pep_index].getIdentifier() == "xtandem";
    const PeptideHit& hit = peptides[pep_index].getHits()[hit_index];
    valid = valid && hit.getScore() >= 0.0 && hit.getScore() <= 1.0 && hit.metaValueExists("XTandem_score");
  }
  TEST_EQUAL(only_xtandem, true)
  TEST_EQUAL(valid, true)

  // separate target and decoy searches of the same spectra file, stored as separate runs: top hits compete per spectrum
  {
    std::vector<ProteinIdentification> runs(2);
    runs[0].setIdentifier("target_search");
    runs[1].setIdentifier("decoy_search");
    for (ProteinIdentification& run : runs)
    {
      run.setSearchEngine("Comet");
      run.setPrimaryMSRunPath({"/data/sample.mzML"});
    }
    PeptideIdentificationList separate;
    for (Size i = 0; i < 20; ++i)
    {
      PeptideIdentification target = makeIdentification("target_search", "expect", false, "scan=" + StringUtils::toStr(i));
      target.setHits({makeHit(std::pow(10.0, -double(i % 5)), 2, "target")});
      PeptideIdentification decoy = makeIdentification("decoy_search", "expect", false, "scan=" + StringUtils::toStr(i));
      decoy.setHits({makeHit(0.5, 2, "decoy")}); // e-value 0.5: the target wins unless its e-value is 1 (i % 5 == 0)
      separate.push_back(target);
      separate.push_back(decoy);
    }
    auto competed = PosteriorErrorProbabilityModel::extractAndTransformScores(runs, separate, false, true);
    TEST_EQUAL(competed["COMET"][1].size(), 16)
    TEST_EQUAL(competed["COMET"][2].size(), 4)
    TEST_EQUAL(competed["COMET"][0].size(), 40) // the unsupervised scores contain all hits
  }

  // supported engine: one score type per engine; either the dedicated score (hits without it are skipped) or, if no hit has it, the main score
  {
    std::vector<ProteinIdentification> runs(2);
    runs[0].setIdentifier("with_evalues");
    runs[1].setIdentifier("without_evalues");
    runs[0].setSearchEngine("XTandem");
    runs[1].setSearchEngine("XTandem");
    PeptideIdentificationList ids;
    for (Size i = 0; i < 4; ++i)
    {
      PeptideIdentification pep = makeIdentification("with_evalues", "XTandem", true, "");
      PeptideHit hit = makeHit(40.0 + double(i), 2, "target");
      if (i % 2 == 0) { hit.setMetaValue("E-Value", 1e-5); }
      pep.setHits({hit});
      ids.push_back(pep);
      PeptideIdentification pep2 = makeIdentification("without_evalues", "XTandem", true, "");
      pep2.setHits({makeHit(40.0 + double(i), 3, "target")});
      ids.push_back(pep2);
    }
    // both runs belong to the same engine (and model): only the hits with e-values are used, no mix of e-values and hyperscores
    auto scores = PosteriorErrorProbabilityModel::extractAndTransformScores(runs, ids, false, false);
    TEST_EQUAL(scores["XTANDEM"][0].size(), 2)
    TEST_REAL_SIMILAR(*std::max_element(scores["XTANDEM"][0].begin(), scores["XTANDEM"][0].end()), 5.0)
    // without any e-values: hyperscore (main score) for all hits
    std::vector<ProteinIdentification> run_without_evalues(1, runs[1]);
    auto main_scores = PosteriorErrorProbabilityModel::extractAndTransformScores(run_without_evalues, ids, false, false);
    TEST_EQUAL(main_scores["XTANDEM"][0].size(), 4)
    TEST_REAL_SIMILAR(*std::max_element(main_scores["XTANDEM"][0].begin(), main_scores["XTANDEM"][0].end()), 43.0)
  }

  // MS-GF+: a single score type for all hits (the most preferred one available), even if hits have different aliases
  {
    std::vector<ProteinIdentification> runs(1);
    runs[0].setIdentifier("msgf");
    runs[0].setSearchEngine("MSGFPlus");
    PeptideIdentificationList ids;
    for (Size i = 0; i < 4; ++i)
    {
      PeptideIdentification pep = makeIdentification("msgf", "SpecEValue", false, "scan=" + StringUtils::toStr(i));
      PeptideHit hit = makeHit(1e-12, 2, "target");
      if (i < 3)
      {
        hit.setMetaValue("MS:1002053", 1e-5); // EValue (preferred); the last hit only has a SpecEValue
      }
      pep.setHits({hit});
      ids.push_back(pep);
    }
    auto scores = PosteriorErrorProbabilityModel::extractAndTransformScores(runs, ids, false, false);
    TEST_EQUAL(scores["MSGFPLUS"][0].size(), 3)
    TEST_REAL_SIMILAR(*std::max_element(scores["MSGFPLUS"][0].begin(), scores["MSGFPLUS"][0].end()), 5.0)
  }

  // competition per spectrum also between hits of different charges (with split_charge), but not if the spectra file is ambiguous
  {
    std::vector<ProteinIdentification> runs(2);
    runs[0].setIdentifier("target_search");
    runs[1].setIdentifier("decoy_search");
    for (ProteinIdentification& run : runs)
    {
      run.setSearchEngine("Comet");
      run.setPrimaryMSRunPath({"/data/sample.mzML"});
    }
    PeptideIdentificationList separate;
    for (Size i = 0; i < 20; ++i)
    {
      PeptideIdentification target = makeIdentification("target_search", "expect", false, "scan=" + StringUtils::toStr(i));
      target.setHits({makeHit(1e-3, 2, "target")});
      PeptideIdentification decoy = makeIdentification("decoy_search", "expect", false, "scan=" + StringUtils::toStr(i));
      decoy.setHits({makeHit(i < 5 ? 1e-4 : 1e-2, 3, "decoy")}); // other charge state; the decoy wins for 5 spectra
      separate.push_back(target);
      separate.push_back(decoy);
    }
    auto by_charge = PosteriorErrorProbabilityModel::extractAndTransformScores(runs, separate, true, true);
    TEST_EQUAL(by_charge["COMET,2"][1].size(), 15)
    TEST_EQUAL(by_charge["COMET,2"][2].size(), 0)
    TEST_EQUAL(by_charge["COMET,3"][1].size(), 0)
    TEST_EQUAL(by_charge["COMET,3"][2].size(), 5)
    TEST_EQUAL(by_charge["COMET,3"][0].size(), 20) // all top hits for the unsupervised model

    // one run with two spectra files, identifications without id_merge_index: the file of a spectrum is unknown, no competition
    std::vector<ProteinIdentification> merged(1);
    merged[0].setIdentifier("merged");
    merged[0].setSearchEngine("Comet");
    merged[0].setPrimaryMSRunPath({"/data/a.mzML", "/data/b.mzML"});
    PeptideIdentificationList ids;
    for (Size i = 0; i < 10; ++i)
    {
      PeptideIdentification first = makeIdentification("merged", "expect", false, "scan=" + StringUtils::toStr(i));
      first.setHits({makeHit(1e-3, 2, "target")});
      PeptideIdentification second = makeIdentification("merged", "expect", false, "scan=" + StringUtils::toStr(i)); // same native ID, other file
      second.setHits({makeHit(1e-2, 2, "decoy")});
      ids.push_back(first);
      ids.push_back(second);
    }
    auto not_competed = PosteriorErrorProbabilityModel::extractAndTransformScores(merged, ids, false, true);
    TEST_EQUAL(not_competed["COMET"][1].size(), 10)
    TEST_EQUAL(not_competed["COMET"][2].size(), 10)
  }

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



