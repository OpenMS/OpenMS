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

// Test that the fitter constrains RT to be within the data range (fixes issue #6239)
// This reproduces the scenario where FeatureFinderMRM reported RT outside range
// The key: the TRUE peak center is BEYOND the truncation point, so without
// constraints the optimizer would naturally try to fit RT outside the data range.
START_SECTION([EXTRA] fit1d_rt_upper_bound)
{
  // Create a peak where the TRUE center is BEYOND the data range
  // This simulates a chromatogram that was cut off mid-peak
  // True RT=210, but data only goes to 200 - we only see the rising edge
  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 100.0);
  tmp.setValue("bounding_box:max", 300.0);  // Model extends to 300
  tmp.setValue("statistics:mean", 200.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 210.0);  // True peak is BEYOND truncation at 200
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  // Truncate the data at RT=200 - we only see the rising edge of the peak
  // Without constraints, the optimizer would fit RT=210 (true value, but outside data)
  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() <= 200.0)
    {
      truncated_samples.push_back(p);
    }
  }

  // Fit - boundary constraints are automatically applied based on data range
  EmgFitter1D ef;
  std::unique_ptr<InterpolationModel> em_fitted;
  ef.fit1d(truncated_samples, em_fitted);

  double fitted_rt = (double)em_fitted->getParameters().getValue("emg:retention");

  std::cout << "Upper bound test: fitted RT = " << fitted_rt
            << " (data range: 100-200, true peak: 210)" << std::endl;

  // The fitted RT MUST be within the data range [100, 200]
  // Without the fix, this would fail because optimizer would find RT=210
  TEST_EQUAL(fitted_rt >= 100.0, true)
  TEST_EQUAL(fitted_rt <= 200.0, true)
}
END_SECTION

// Test boundary constraint at the lower edge
// Similar to upper bound: true peak is BELOW the data range start
START_SECTION([EXTRA] fit1d_rt_lower_bound)
{
  // Create a peak where the TRUE center is BELOW the data range
  // True RT=2, but data only starts at 10 - we only see the falling edge
  EmgModel em;
  em.setInterpolationStep(1.0);
  Param tmp;
  tmp.setValue("bounding_box:min", 0.0);
  tmp.setValue("bounding_box:max", 200.0);
  tmp.setValue("statistics:mean", 10.0);
  tmp.setValue("statistics:variance", 2.0);
  tmp.setValue("emg:height", 50000.0);
  tmp.setValue("emg:width", 5.0);
  tmp.setValue("emg:symmetry", 5.0);
  tmp.setValue("emg:retention", 2.0);  // True peak is BELOW data start at 10
  em.setParameters(tmp);

  EmgModel::SamplesType full_samples;
  em.getSamples(full_samples);

  // Truncate the data at RT=10 - we only see the falling edge of the peak
  EmgModel::SamplesType truncated_samples;
  for (const auto& p : full_samples)
  {
    if (p.getPos() >= 10.0)
    {
      truncated_samples.push_back(p);
    }
  }

  EmgFitter1D ef;
  std::unique_ptr<InterpolationModel> em_fitted;
  ef.fit1d(truncated_samples, em_fitted);

  double fitted_rt = (double)em_fitted->getParameters().getValue("emg:retention");

  std::cout << "Lower bound test: fitted RT = " << fitted_rt
            << " (data range: 10-200, true peak: 2)" << std::endl;

  // The fitted RT must be within the truncated data range [10, 200]
  // Without the fix, this might fail because optimizer would find RT<10
  TEST_EQUAL(fitted_rt >= 10.0, true)
  TEST_EQUAL(fitted_rt <= 200.0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



