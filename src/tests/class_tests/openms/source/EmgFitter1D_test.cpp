// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>

///////////////////////////

using namespace OpenMS;
using namespace std;
using boost::normal_distribution;

START_TEST(EmgFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EmgFitter1D* ptr = nullptr;
EmgFitter1D* nullPointer = nullptr;
START_SECTION(EmgFitter1D())
{
  ptr = new EmgFitter1D();
  TEST_EQUAL(ptr->getName(), "EmgFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((EmgFitter1D(const  EmgFitter1D &source)))
	EmgFitter1D emgf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	emgf1.setParameters(param);

	EmgFitter1D emgf2(emgf1);
  EmgFitter1D emgf3;
	emgf3.setParameters(param);
  emgf1 = EmgFitter1D();
	TEST_EQUAL(emgf3.getParameters(), emgf2.getParameters())
END_SECTION

START_SECTION((virtual ~EmgFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual EmgFitter1D& operator=(const  EmgFitter1D &source)))
	EmgFitter1D emgf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
	emgf1.setParameters(param);

  EmgFitter1D emgf2;
  emgf2 = emgf1;

  EmgFitter1D emgf3;
	emgf3.setParameters(param);

  emgf1 = EmgFitter1D();
	TEST_EQUAL(emgf3.getParameters(), emgf2.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	//create data via a model
	EmgModel em;
  em.setInterpolationStep(0.2);
  Param tmp;
  tmp.setValue("bounding_box:min", 678.9);
  tmp.setValue("bounding_box:max", 789.0);
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("emg:height",100000.0);
  tmp.setValue("emg:width",5.0);
  tmp.setValue("emg:symmetry",5.0);
  tmp.setValue("emg:retention",725.0);
  em.setParameters(tmp);
  EmgModel::SamplesType samples;
  em.getSamples(samples);
  //fit the data
  EmgFitter1D ef = EmgFitter1D();
	std::unique_ptr<InterpolationModel> em_fitted;
	EmgFitter1D::QualityType correlation = ef.fit1d(samples, em_fitted);

	//check the fitted model on the exact data
	TEST_REAL_SIMILAR(correlation, 1);
	TEST_REAL_SIMILAR((double)em_fitted->getParameters().getValue("emg:height"), 100000.0);
	TEST_REAL_SIMILAR((double)em_fitted->getParameters().getValue("emg:width"), 5.0);
	TEST_REAL_SIMILAR((double)em_fitted->getParameters().getValue("emg:symmetry"), 5.0);
	TEST_REAL_SIMILAR((double)em_fitted->getParameters().getValue("emg:retention"), 725.0);


	//shake the samples a little with varying variance (difficult test for fitter)
	EmgModel::SamplesType unexact_samples;
	boost::mt19937 rng;//random number generator
	for(unsigned i=0; i<samples.size(); ++i)
  {

	  boost::normal_distribution<double> distInt (samples.at(i).getIntensity(), samples.at(i).getIntensity()/100); //use sample intensity as mean
	  Peak1D p (samples.at(i).getPosition()[0], distInt(rng));
	  std::cout << "point: (" <<  samples.at(i).getPosition()[0] << ", " << samples.at(i).getIntensity() << ") -> (" <<  p.getPosition()[0] << ", " << p.getIntensity() << ")" << std::endl;
	  unexact_samples.push_back(p);
  }
  //fit the data
  EmgFitter1D ef1 = EmgFitter1D();
  std::unique_ptr<InterpolationModel> em_fitted1;
  EmgFitter1D::QualityType correlation1 = ef1.fit1d(unexact_samples, em_fitted1);
  TOLERANCE_RELATIVE(1.01)
	TEST_REAL_SIMILAR(correlation1, 1);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:height"), 100000.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:width"), 5.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:symmetry"), 5.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:retention"), 725.0);

END_SECTION

// Test that demonstrates the RT-outside-range bug when bounds are disabled
// This reproduces issue #6239 where FeatureFinderMRM reported RT outside range
START_SECTION([EXTRA] fit1d_rt_without_bounds_can_exceed_range)
{
  // Create a peak near the upper edge of the data range
  // The data is truncated after RT=200, simulating a chromatogram that ends abruptly
  // Without the boundary constraint, the optimizer can shift RT beyond 200

  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 100.0);
  tmp.setValue("bounding_box:max", 300.0);  // Model goes to 300
  tmp.setValue("statistics:mean", 190.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 195.0);  // Peak near upper edge
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  // Truncate the data at RT=200 to simulate chromatogram ending abruptly
  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() <= 200.0)
    {
      truncated_samples.push_back(p);
    }
  }

  // Test WITHOUT boundary constraints (old behavior) - RT can go outside range
  EmgFitter1D ef_no_bounds;
  Param p_no_bounds;
  p_no_bounds.setValue("fit_bounds:enabled", "false");
  ef_no_bounds.setParameters(p_no_bounds);

  std::unique_ptr<InterpolationModel> em_fitted_no_bounds;
  EmgFitter1D::QualityType correlation_no_bounds = ef_no_bounds.fit1d(truncated_samples, em_fitted_no_bounds);

  double fitted_rt_no_bounds = (double)em_fitted_no_bounds->getParameters().getValue("emg:retention");

  std::cout << "WITHOUT bounds: fitted RT = " << fitted_rt_no_bounds
            << " (data range: 100-200, true peak: 195)" << std::endl;

  // Without bounds, the optimizer may report RT outside the data range
  // We don't assert a specific value here since the unbounded behavior is data-dependent
  // but we document that it CAN exceed the range
  TEST_EQUAL(correlation_no_bounds > 0.5, true)  // Should still produce some fit
}
END_SECTION

// Test that the fitter constrains RT to be within the data range when bounds are enabled
START_SECTION([EXTRA] fit1d_rt_with_bounds_stays_in_range)
{
  // Same setup as above but with bounds enabled (default)
  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 100.0);
  tmp.setValue("bounding_box:max", 300.0);
  tmp.setValue("statistics:mean", 190.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 195.0);
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() <= 200.0)
    {
      truncated_samples.push_back(p);
    }
  }

  // Test WITH boundary constraints (new default behavior) - RT stays in range
  EmgFitter1D ef_with_bounds;
  // fit_bounds:enabled defaults to true, but let's be explicit
  Param p_with_bounds;
  p_with_bounds.setValue("fit_bounds:enabled", "true");
  ef_with_bounds.setParameters(p_with_bounds);

  std::unique_ptr<InterpolationModel> em_fitted_with_bounds;
  EmgFitter1D::QualityType correlation_with_bounds = ef_with_bounds.fit1d(truncated_samples, em_fitted_with_bounds);

  double fitted_rt_with_bounds = (double)em_fitted_with_bounds->getParameters().getValue("emg:retention");

  std::cout << "WITH bounds: fitted RT = " << fitted_rt_with_bounds
            << " (data range: 100-200, true peak: 195)" << std::endl;

  // With bounds, the fitted RT MUST be within the data range [100, 200]
  TEST_EQUAL(fitted_rt_with_bounds >= 100.0, true)
  TEST_EQUAL(fitted_rt_with_bounds <= 200.0, true)

  // The fit should still be reasonable
  TEST_EQUAL(correlation_with_bounds > 0.9, true)
}
END_SECTION

// Test boundary constraint at the lower edge
START_SECTION([EXTRA] fit1d_rt_lower_bound)
{
  // Create a peak near the lower edge of the data range
  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 0.0);
  tmp.setValue("bounding_box:max", 200.0);
  tmp.setValue("statistics:mean", 15.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 10.0);  // Peak near lower edge
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  // Truncate the data at RT=5 at the lower end
  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() >= 5.0)
    {
      truncated_samples.push_back(p);
    }
  }

  // With bounds enabled (default)
  EmgFitter1D ef;
  std::unique_ptr<InterpolationModel> em_fitted;
  EmgFitter1D::QualityType correlation = ef.fit1d(truncated_samples, em_fitted);

  double fitted_rt = (double)em_fitted->getParameters().getValue("emg:retention");

  // The fitted RT must be within the truncated data range [5, 200]
  TEST_EQUAL(fitted_rt >= 5.0, true)
  TEST_EQUAL(fitted_rt <= 200.0, true)

  TEST_EQUAL(correlation > 0.9, true)

  std::cout << "Lower bound test: fitted RT = " << fitted_rt
            << " (data range: 5-200, true peak: 10)" << std::endl;
}
END_SECTION

// Test that default behavior (no explicit parameter set) uses bounds
START_SECTION([EXTRA] fit1d_default_uses_bounds)
{
  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 100.0);
  tmp.setValue("bounding_box:max", 300.0);
  tmp.setValue("statistics:mean", 190.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 195.0);
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() <= 200.0)
    {
      truncated_samples.push_back(p);
    }
  }

  // Default constructor - should use bounds by default
  EmgFitter1D ef_default;
  std::unique_ptr<InterpolationModel> em_fitted;
  ef_default.fit1d(truncated_samples, em_fitted);

  double fitted_rt = (double)em_fitted->getParameters().getValue("emg:retention");

  // Default should constrain RT to data range
  TEST_EQUAL(fitted_rt >= 100.0, true)
  TEST_EQUAL(fitted_rt <= 200.0, true)

  std::cout << "Default behavior test: fitted RT = " << fitted_rt
            << " (expected within 100-200)" << std::endl;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



