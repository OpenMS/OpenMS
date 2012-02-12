// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <vector>

///////////////////////////

START_TEST(TransformationDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


TransformationDescription* ptr = 0;
TransformationDescription* nullPointer = 0;
START_SECTION((TransformationDescription()))
	ptr = new TransformationDescription;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~TransformationDescription()))
	delete ptr;
END_SECTION

TransformationDescription::DataPoints data;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 3.0));

START_SECTION((TransformationDescription(const DataPoints& data)))
{
	ptr = new TransformationDescription(data);
  TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getDataPoints() == data, true);
	delete ptr;
}
END_SECTION

START_SECTION((const DataPoints& getDataPoints() const))
{
	TransformationDescription td;
	TEST_EQUAL(td.getDataPoints().empty(), true);
}
END_SECTION

START_SECTION((void setDataPoints(const DataPoints& data)))
{
	TransformationDescription td;
	td.fitModel("identity", Param());
	TEST_EQUAL(td.getModelType(), "identity");
	td.setDataPoints(data);
	// setting data points clears the model:
	TEST_EQUAL(td.getModelType(), "none");
	TEST_EQUAL(td.getDataPoints().size(), 2)
	TEST_EQUAL(td.getDataPoints() == data, true);

	TransformationDescription::DataPoints empty;
  empty.clear();
	td.setDataPoints(empty);
	TEST_EQUAL(td.getDataPoints().empty(), true)
}
END_SECTION

START_SECTION((DoubleReal apply(DoubleReal value) const))
{
	TransformationDescription td;
	TEST_EQUAL(td.apply(-0.5), -0.5);
	TEST_EQUAL(td.apply(1000), 1000);
	// tested further together with "fitModel"
}
END_SECTION

START_SECTION((const String& getModelType() const))
{
	TransformationDescription td;
	TEST_EQUAL(td.getModelType(), "none");
}
END_SECTION

START_SECTION((static void getModelTypes(StringList& result)))
{
	StringList result;
	TransformationDescription::getModelTypes(result);
	TEST_EQUAL(result.size(), 3);
	TEST_EQUAL(result[0], "linear");
	TEST_EQUAL(result[1], "b_spline");
	TEST_EQUAL(result[2], "interpolated");
}
END_SECTION


START_SECTION((void fitModel(const String& model_type, const Param& params=Param())))
{
	TransformationDescription td(data);
	Param params;
	td.fitModel("linear", params);
	TEST_EQUAL(td.getModelType(), "linear");
	TEST_REAL_SIMILAR(td.apply(0.0), 1.0);
	TEST_REAL_SIMILAR(td.apply(0.5), 2.0);
	TEST_REAL_SIMILAR(td.apply(1.0), 3.0);
	
  // special model type for reference files:
	td.fitModel("identity", Param());
	TEST_EQUAL(td.getModelType(), "identity");
	TEST_REAL_SIMILAR(td.apply(0.0), 0.0);
	TEST_REAL_SIMILAR(td.apply(0.5), 0.5);
	TEST_REAL_SIMILAR(td.apply(1.0), 1.0);
	// can't fit a different model to an "identity" transformation:
	td.fitModel("linear", params);
	TEST_EQUAL(td.getModelType(), "identity");
}
END_SECTION

START_SECTION((void getModelParameters(Param& params) const))
{
	TransformationDescription td;
	Param params;
	td.getModelParameters(params);
	TEST_EQUAL(params, Param());
	params.setValue("slope", 2.5);
	params.setValue("intercept", -100.0);
	const Param const_params = params;
	td.fitModel("linear", const_params);
	td.getModelParameters(params);
	TEST_EQUAL(params, const_params);
}
END_SECTION

START_SECTION((TransformationDescription(const TransformationDescription& rhs)))
{
	TransformationDescription td(data);
	td.fitModel("linear", Param());
	TransformationDescription td2 = td;
	TEST_EQUAL(td.getModelType(), td2.getModelType());
	TEST_EQUAL(td.getDataPoints() == td2.getDataPoints(), true);
	Param params, params2;
	td.getModelParameters(params);
	td2.getModelParameters(params2);
	TEST_EQUAL(params, params2);
}
END_SECTION

START_SECTION((TransformationDescription& operator=(const TransformationDescription& rhs)))
{
	TransformationDescription td(data);
	td.fitModel("linear", Param());
	TransformationDescription td2;
	td2 = td;
	TEST_EQUAL(td.getModelType(), td2.getModelType());
	TEST_EQUAL(td.getDataPoints() == td2.getDataPoints(), true);
	Param params, params2;
	td.getModelParameters(params);
	td2.getModelParameters(params2);
	TEST_EQUAL(params, params2);
}
END_SECTION

START_SECTION((void invert()))
{
	TransformationDescription td;
	DoubleReal value = 57.12;
	Param params;

	// test null transformation:
	td.fitModel("none", params);
	td.invert();
	TEST_EQUAL(td.getModelType(), "none");

	// test linear transformation:
	params.setValue("slope", 2.0);
	params.setValue("intercept", 47.12);
	td.fitModel("linear", params);

	td.invert();
	TEST_EQUAL(td.getModelType(), "linear");
	td.getModelParameters(params);
	TEST_REAL_SIMILAR(params.getValue("slope"), 0.5);
	TEST_REAL_SIMILAR(params.getValue("intercept"), -23.56);
	TEST_REAL_SIMILAR(td.apply(value), 5.0);

	// test inversion of data points:
	td.setDataPoints(data);
	td.invert();
	TEST_EQUAL(td.getDataPoints()[0].first, data[0].second);
	TEST_EQUAL(td.getDataPoints()[0].second, data[0].first);
	TEST_EQUAL(td.getDataPoints()[1].first, data[1].second);
	TEST_EQUAL(td.getDataPoints()[1].second, data[1].first);
	td.invert();
	TEST_EQUAL(td.getDataPoints() == data, true);

	// test interpolated-linear transformation:
	params.clear();
	params.setValue("interpolation_type", "linear");
	td.fitModel("interpolated", params);
	td.invert();
	TEST_EQUAL(td.getModelType(), "interpolated");
	// pairs have changed...
	TEST_EQUAL(td.getDataPoints() != data, true);
	td.invert();
	// ... now they're back to the original:
	TEST_EQUAL(td.getDataPoints() == data, true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
