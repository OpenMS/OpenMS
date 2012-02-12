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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/TransformationXMLFile.h>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationXMLFile* ptr = 0;
TransformationXMLFile* nullPointer = 0;

START_SECTION((TransformationXMLFile()))
{
	ptr = new TransformationXMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer);
}
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
{
	TransformationXMLFile f;
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML")),true);	
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML")),true);	
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_3.trafoXML")),false);	
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML")),true);	
}
END_SECTION


START_SECTION(void load(const String& filename, TransformationDescription& transformation))
{
	TransformationDescription trafo;
	TransformationXMLFile trafo_xml;
	Param params;

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"), trafo);
	TEST_STRING_EQUAL(trafo.getModelType(), "none");
	trafo.getModelParameters(params);
	TEST_EQUAL(params.empty(), true);

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"),trafo);
	TEST_STRING_EQUAL(trafo.getModelType(), "linear");
	trafo.getModelParameters(params);
	TEST_EQUAL(params.size(), 2);
	TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
	TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML"), trafo);
	TEST_STRING_EQUAL(trafo.getModelType(), "interpolated");
	trafo.getModelParameters(params);
	TEST_EQUAL(params.getValue("interpolation_type"), "linear");
	TEST_EQUAL(trafo.getDataPoints().size(), 3);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[0].first, 1.2);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[1].first, 2.2);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[2].first, 3.2);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[0].second, 5.2);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[1].second, 6.25);
	TEST_REAL_SIMILAR(trafo.getDataPoints()[2].second, 7.3);
}
END_SECTION

START_SECTION(void store(String filename, const TransformationDescription& transformation))
{
	TransformationDescription trafo, trafo2;
	TransformationXMLFile trafo_xml;

	String tmp_file_none;
	Param params;
	trafo.fitModel("none", params);
	NEW_TMP_FILE(tmp_file_none);
	trafo_xml.store(tmp_file_none, trafo);
	trafo_xml.load(tmp_file_none, trafo2);
	TEST_STRING_EQUAL(trafo2.getModelType(), "none");
	trafo2.getModelParameters(params);
	TEST_EQUAL(params.empty(), true);
	{
	  DoubleReal pre_image = 234255132.43212;
	  DoubleReal image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
	}
	String tmp_file_linear;
	NEW_TMP_FILE(tmp_file_linear);
	params.setValue("slope", 3.141592653589793238);
	params.setValue("intercept", 2.718281828459045235);
	trafo.fitModel("linear", params);
	trafo_xml.store(tmp_file_linear, trafo);
	trafo_xml.load(tmp_file_linear, trafo2);
	TEST_STRING_EQUAL(trafo.getModelType(), "linear");
	params.clear();
	trafo2.getModelParameters(params);
	TEST_EQUAL(params.size(), 2);
	TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
	TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);
  {
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

	String tmp_file_pairs;
	NEW_TMP_FILE(tmp_file_pairs);
	TransformationDescription::DataPoints pairs;
	pairs.push_back(make_pair(1.2f,5.2f));
	pairs.push_back(make_pair(2.2f,6.25f));
	pairs.push_back(make_pair(3.2f,7.3f));
	trafo.setDataPoints(pairs);
	params.clear();
	params.setValue("interpolation_type", "linear");
	trafo.fitModel("interpolated", params);
	trafo_xml.store(tmp_file_pairs, trafo);
	trafo_xml.load(tmp_file_pairs, trafo2);
	TEST_STRING_EQUAL(trafo2.getModelType(), "interpolated");
	trafo2.getModelParameters(params);
	TEST_EQUAL(params.size(), 1);
	TEST_STRING_EQUAL(params.getValue("interpolation_type"), "linear");
	TEST_EQUAL(trafo2.getDataPoints().size(), 3);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].first, 1.2);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[1].first, 2.2);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[2].first, 3.2);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].second, 5.2);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[1].second, 6.25);
	TEST_REAL_SIMILAR(trafo2.getDataPoints()[2].second, 7.3);
  {
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

  TEST_EXCEPTION(Exception::IllegalArgument, trafo.fitModel("mumble_pfrwoarpfz"));

	String tmp_file_bspline;
	NEW_TMP_FILE(tmp_file_bspline);
  pairs.clear();
  pairs.push_back(make_pair(1.2f,5.2f));
  pairs.push_back(make_pair(3.2f,7.3f));
  pairs.push_back(make_pair(2.2f,6.25f));
  pairs.push_back(make_pair(2.2f,3.1f));
  pairs.push_back(make_pair(2.2f,7.25f));
  pairs.push_back(make_pair(3.0f,8.5f));
  pairs.push_back(make_pair(3.1f,4.7f));
  pairs.push_back(make_pair(1.7f,6.0f));
  pairs.push_back(make_pair(2.9f,4.7f));
  pairs.push_back(make_pair(4.2f,5.0f));
  pairs.push_back(make_pair(3.7f,-2.4f));
  trafo.setDataPoints(pairs);
	params.clear();
	params.setValue("num_breakpoints", 4);
	params.setValue("break_positions", "uniform");
	trafo.fitModel("b_spline", params);
  trafo_xml.store(tmp_file_pairs,trafo);
  trafo_xml.load(tmp_file_pairs,trafo2);
  TEST_STRING_EQUAL(trafo2.getModelType(), "b_spline");
	params.clear();
	trafo2.getModelParameters(params);
  TEST_EQUAL(params.getValue("num_breakpoints"), 4);
  TEST_EQUAL(params.getValue("break_positions"), "uniform");
  TEST_EQUAL(params.size(), 2);
  TEST_EQUAL(trafo2.getDataPoints().size(), 11);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].first, 1.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].second, 5.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].first, 3.7);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].second, -2.4);
  for ( Int breaks = 0; breaks < 10; ++breaks)
  {
    if ( breaks == 1) continue;
		params.setValue("num_breakpoints", breaks);
		trafo.fitModel("b_spline", params);
    STATUS("breaks: " << breaks);
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
