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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/EGHFitter1D.h>
#include <OpenMS/SIMULATION/EGHModel.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <cmath>//toberemoved
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EGHFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EGHFitter1D* ptr = nullptr;
EGHFitter1D* nullPointer = nullptr;
START_SECTION(EGHFitter1D())
{
	ptr = new EGHFitter1D();
	TEST_EQUAL(ptr->getName(), "EGHFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((EGHFitter1D(const EGHFitter1D &source)))
{
  EGHFitter1D eghf1;

  Param param;
  param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
  eghf1.setParameters(param);

  EGHFitter1D eghf2(eghf1);
  EGHFitter1D eghf3;
  eghf3.setParameters(param);
  eghf1 = EGHFitter1D();
  TEST_EQUAL(eghf3.getParameters(), eghf2.getParameters())
}
END_SECTION

START_SECTION((virtual ~EGHFitter1D()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual EGHFitter1D& operator=(const EGHFitter1D &source)))
{
  EGHFitter1D eghf1;

  Param param;
  param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "max_iteration", 500 );
  param.setValue( "deltaAbsError", 0.0001 );
  param.setValue( "deltaRelError", 0.0001 );
  eghf1.setParameters(param);

  EGHFitter1D eghf2;
  eghf2 = eghf1;
  EGHFitter1D eghf3;
  eghf3.setParameters(param);
  eghf1 = EGHFitter1D();
  TEST_EQUAL(eghf3.getParameters(), eghf2.getParameters())
}
END_SECTION

START_SECTION((QualityType fit1d(const RawDataArrayType &range, InterpolationModel *&model)))
{
  EGHModel base_model;

  Param tmp;
  tmp.setValue( "statistics:variance", 1.0 );
  tmp.setValue( "statistics:mean", 1000.0 );

  tmp.setValue( "egh:height", 1000.0f );
  tmp.setValue( "egh:retention", 1000.0f);

  tmp.setValue("egh:guess_parameter", "false"); // disable guessing of parameters from A/B
  tmp.setValue("egh:tau", 72.0);
  tmp.setValue("egh:sigma_square", 3606.0);

  //p.setValue("egh:A", 50.0);
  //p.setValue("egh:B", 90.0);

  base_model.setParameters(tmp);

  // Raw data point type
  typedef Peak1D PeakType;
  // Raw data container type using for the temporary storage of the input data
  typedef std::vector<PeakType> RawDataArrayType;

  RawDataArrayType data_to_fit;

  for (double x = 800.0; x < 1200.0; x += 0.1)
  {
    PeakType p;
    p.setPos(x);
    p.setIntensity(base_model.getIntensity(x));
    data_to_fit.push_back(p);
  }

  // make some noise
  boost::random::mt19937 rnd_gen_ (0.0);
  boost::uniform_real<float> udist (-0.1, 0.1);
  for (Size i = 0; i < data_to_fit.size(); ++i)
  {
    float distort = std::exp(udist(rnd_gen_));
    data_to_fit[i].setIntensity(data_to_fit[i].getIntensity()
        * distort);
  }

  double egh_quality;
  Param egh_param;
  EGHFitter1D egh_fitter;

  // Set parameter for fitter
  egh_fitter.setParameters( egh_param );

  InterpolationModel* fitted_egh_model = nullptr;

  // Construct model for rt
  egh_quality = egh_fitter.fit1d(data_to_fit, fitted_egh_model);

  TOLERANCE_ABSOLUTE(5.0)
  TEST_REAL_SIMILAR(egh_quality, 0.996313)
  TEST_REAL_SIMILAR(fitted_egh_model->getParameters().getValue("egh:tau"), 72.0)
  TEST_REAL_SIMILAR(fitted_egh_model->getParameters().getValue("egh:sigma_square"), 3606.0)

}
END_SECTION

START_SECTION((static Fitter1D* create()))
{
  Fitter1D* ptr = EGHFitter1D::create();
  TEST_EQUAL(ptr->getName(), "EGHFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(EGHFitter1D::getProductName(),"EGHFitter1D")
  TEST_EQUAL(EGHFitter1D().getName(),"EGHFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



