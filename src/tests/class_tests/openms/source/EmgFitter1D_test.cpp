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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



