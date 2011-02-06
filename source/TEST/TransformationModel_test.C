// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>

///////////////////////////

START_TEST(TransformationDescription, "$Id: TransformationModel_test.C 7583 2010-11-19 17:20:48Z hendrikweisser $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModel* ptr = 0;
START_SECTION((TransformationModel()))
{
	ptr = new TransformationModel();
	TEST_NOT_EQUAL(ptr, 0);
}
END_SECTION

START_SECTION((TransformationModel(const DataPoints&, const Param&)))
{
	ptr = new TransformationModel(TransformationModel::DataPoints(), Param());
	TEST_NOT_EQUAL(ptr, 0);
}
END_SECTION

START_SECTION((~TransformationModel()))
{
	delete ptr;
}
END_SECTION

TransformationModel::DataPoints data, empty;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 2.0));
data.push_back(make_pair(1.0, 4.0));

START_SECTION((DoubleReal evaluate(DoubleReal)))
{
	// null model (identity):
	TransformationModel tm;
	TEST_REAL_SIMILAR(tm.evaluate(-3.14159), -3.14159);
	TEST_REAL_SIMILAR(tm.evaluate(0.0), 0.0);
	TEST_REAL_SIMILAR(tm.evaluate(12345678.9), 12345678.9);

	// linear model:
	TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelLinear 
								 lm(empty, Param())); // need data

	TransformationModelLinear lm(data, Param());
	TEST_REAL_SIMILAR(lm.evaluate(-0.5), 0.0);
	TEST_REAL_SIMILAR(lm.evaluate(0.0), 1.0);
	TEST_REAL_SIMILAR(lm.evaluate(0.5), 2.0);
	TEST_REAL_SIMILAR(lm.evaluate(1.0), 3.0);
	TEST_REAL_SIMILAR(lm.evaluate(1.5), 4.0);
	
	// interpolation model:
	TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelInterpolated
								 im(empty, Param())); // need data
	
  data.push_back(make_pair(2.0, 2.0));
	TransformationModelInterpolated im(data, Param());
	// interpolation:
	TEST_REAL_SIMILAR(im.evaluate(0.0), 1.0);
	TEST_REAL_SIMILAR(im.evaluate(0.5), 2.0);
	TEST_REAL_SIMILAR(im.evaluate(1.0), 3.0);
	TEST_REAL_SIMILAR(im.evaluate(1.5), 2.5);
	TEST_REAL_SIMILAR(im.evaluate(2.0), 2.0);
	// extrapolation:
	TEST_REAL_SIMILAR(im.evaluate(-0.5), 0.75);
	TEST_REAL_SIMILAR(im.evaluate(2.5), 2.25);

	// B-spline model:
	TEST_EXCEPTION(Exception::IllegalArgument, 
								 TransformationModelBSpline bm(data, Param())); // need param.
	Param params;
	params.setValue("num_breakpoints", 4);
	TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelBSpline 
								 bm(empty, params)); // need data

	data.clear();
  data.push_back(make_pair(1.2, 5.2));
  data.push_back(make_pair(3.2, 7.3));
  data.push_back(make_pair(2.2, 6.25));
  data.push_back(make_pair(2.2, 3.1));
  data.push_back(make_pair(2.2, 7.25));
  data.push_back(make_pair(3.0, 8.5));
  data.push_back(make_pair(3.1, 4.7));
  data.push_back(make_pair(1.7, 6.0));
  data.push_back(make_pair(2.9, 4.7));
  data.push_back(make_pair(4.2, 5.0));
  data.push_back(make_pair(3.7, -2.4));

	TransformationModelBSpline bm(data, params);

#if 0
	// Since the numbers in this test were verified by manual (in fact, visual)
	// inspection, here is a recipe how this was done:
	//
	// To grep for output, "pairs:" and "spline:", you might use a command line
	// like this:
	// make TransformationDescription_test && ctest -V -R TransformationDescription_test &&  ctest -V -R TransformationDescription_test | grep pairs: > points.dat && ctest -V -R TransformationDescription_test | grep spline: > bla.dat
	// To have a look at the results using gnuplot:
	// gnuplot -  # (to start gnuplot)
	// plot './bla.dat' u 5:6, './points.dat' u 5:6

	for (UInt i = 0; i < data.size(); ++i)
	{
		STATUS("data: " << data[i].first << " " << data[i].second);
	}

	for (Int i = -10; i <= 60; i+=5)
	{
		DoubleReal value = i;
		value /= 10;
		DoubleReal image = bm.evaluate(value);
		STATUS("spline: " << value << " " << image);
	}
#endif

	DoubleReal sample_values[] = {-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6};
	DoubleReal sample_images[] = {-14.3519415123977, -9.91557518507088, -5.4792088577441, -1.04284253041731, 3.39352379690948, 6.4561466812738, 5.4858954730427, 6.14659387774751, 6.77299727168147, 0.646024122587505, -1.13062259235381, 18.3842099268184, 40.7826815802615, 63.1811532337045, 85.5796248871476};
	for (Size i = 0; i < sizeof(sample_values)/sizeof(*sample_values); ++i)
	{
		DoubleReal x = bm.evaluate(sample_values[i]);
		TEST_REAL_SIMILAR(x, sample_images[i]);
	}

}
END_SECTION

START_SECTION((void getParameters(Param&) const))
{
	TransformationModel tm;
	Param p_in, p_out;
	tm.getParameters(p_out);
	TEST_EQUAL(p_out.empty(), true);

	p_in.clear();
	p_in.setValue("symmetric_regression", "true");
	TransformationModelLinear lm(data, p_in);
	lm.getParameters(p_out);
	TEST_EQUAL(p_out, p_in);
	p_in.clear();
	p_in.setValue("slope", 12.3);
	p_in.setValue("intercept", -45.6);
	TransformationModelLinear lm2(empty, p_in);
	lm2.getParameters(p_out);
	TEST_EQUAL(p_out, p_in);

	p_in.clear();
	p_in.setValue("interpolation_type", "polynomial");
	TransformationModelInterpolated im(data, p_in);
	im.getParameters(p_out);
	TEST_EQUAL(p_out, p_in);

	p_in.clear();
	p_in.setValue("num_breakpoints", 4);
	TransformationModelBSpline bm(data, p_in);
	bm.getParameters(p_out);
	TEST_EQUAL(p_out, p_in);
}
END_SECTION

START_SECTION(([EXTRA] void getInterpolationTypes(StringList&) const))
{
	StringList result;
	TransformationModelInterpolated::getInterpolationTypes(result);
	TEST_EQUAL(result.concatenate(","), "linear,polynomial,cspline,akima");
}
END_SECTION

START_SECTION(([EXTRA] void getParameters(DoubleReal&, DoubleReal&)))
{
	Param param;
	param.setValue("slope", 12.3);
	param.setValue("intercept", -45.6);
	TransformationModelLinear lm(empty, param);
	DoubleReal slope, intercept;
	lm.getParameters(slope, intercept);
	TEST_REAL_SIMILAR(param.getValue("slope"), slope);
	TEST_REAL_SIMILAR(param.getValue("intercept"), intercept);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
