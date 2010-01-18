// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/EGHFitter1D.h>
#include <OpenMS/SIMULATION/EGHModel.h>
///////////////////////////

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace OpenMS;
using namespace std;

START_TEST(EGHFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EGHFitter1D* ptr = 0;
START_SECTION(EGHFitter1D())
{
	ptr = new EGHFitter1D();
	TEST_EQUAL(ptr->getName(), "EGHFitter1D")
	TEST_NOT_EQUAL(ptr, 0)
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

  for (DoubleReal x = 800.0; x < 1200.0; x += 0.1)
  {
    PeakType p;
    p.setPos(x);
    p.setIntensity(base_model.getIntensity(x));
    data_to_fit.push_back(p);
  }

  // make some noise
  gsl_rng_default_seed = 0.0;
  gsl_rng* rnd_gen_ = gsl_rng_alloc(gsl_rng_mt19937);
  DoubleReal distortion = 0.1;

  for (Size i = 0; i < data_to_fit.size(); ++i)
  {
    DoubleReal distort = exp(gsl_ran_flat(rnd_gen_, -distortion,
        +distortion));
    data_to_fit[i].setIntensity(data_to_fit[i].getIntensity()
        * distort);
  }

  DoubleReal egh_quality;
  Param egh_param;
  EGHFitter1D egh_fitter;

  // Set parameter for fitter
  egh_fitter.setParameters( egh_param );

  InterpolationModel* fitted_egh_model = 0;

  // Construct model for rt
  egh_quality = egh_fitter.fit1d(data_to_fit, fitted_egh_model);

  TOLERANCE_ABSOLUTE(5.0)
  TEST_REAL_SIMILAR(fitted_egh_model->getParameters().getValue("egh:tau"), 72.0)
  TEST_REAL_SIMILAR(fitted_egh_model->getParameters().getValue("egh:sigma_square"), 3606.0)

}
END_SECTION

START_SECTION((static Fitter1D* create()))
{
  Fitter1D* ptr = EGHFitter1D::create();
  TEST_EQUAL(ptr->getName(), "EGHFitter1D")
  TEST_NOT_EQUAL(ptr, 0)
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



