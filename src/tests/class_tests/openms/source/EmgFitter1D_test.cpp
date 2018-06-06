// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>
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
	InterpolationModel* em_fitted = nullptr;
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
  InterpolationModel* em_fitted1 = nullptr;
  EmgFitter1D::QualityType correlation1 = ef1.fit1d(unexact_samples, em_fitted1);
  TOLERANCE_RELATIVE(1.01)
	TEST_REAL_SIMILAR(correlation1, 1);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:height"), 100000.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:width"), 5.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:symmetry"), 5.0);
  TEST_REAL_SIMILAR((double)em_fitted1->getParameters().getValue("emg:retention"), 725.0);

END_SECTION

START_SECTION((Fitter1D* create()))
{
  Fitter1D* ptr = EmgFitter1D::create();
  TEST_EQUAL(ptr->getName(), "EmgFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((const String getProductName()))
{
  TEST_EQUAL(EmgFitter1D::getProductName(),"EmgFitter1D")
  TEST_EQUAL(EmgFitter1D().getName(),"EmgFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



