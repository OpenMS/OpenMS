// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(FASTAFile, "$Id: TransformationXMLFile_test.C 6054 2009-09-29 10:03:45Z cbielow $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationXMLFile* ptr;
START_SECTION((TransformationXMLFile()))
{
	ptr = new TransformationXMLFile();
	TEST_NOT_EQUAL(ptr,0);
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

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"),trafo);
	TEST_STRING_EQUAL(trafo.getName(),"none");
	TEST_EQUAL(trafo.getParameters().empty(),true);

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"),trafo);
	TEST_STRING_EQUAL(trafo.getName(),"linear");
	TEST_EQUAL(trafo.getParameters().size(),2);
	TEST_REAL_SIMILAR(trafo.getParam("slope"),3.141592653589793238);
	TEST_REAL_SIMILAR(trafo.getParam("intercept"),2.718281828459045235);

	trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML"),trafo);
	TEST_STRING_EQUAL(trafo.getName(),"interpolated_linear");
	TEST_EQUAL(trafo.getParameters().size(),0);
	TEST_EQUAL(trafo.getPairs().size(),3);
	TEST_REAL_SIMILAR(trafo.getPairs()[0].first,1.2);
	TEST_REAL_SIMILAR(trafo.getPairs()[1].first,2.2);
	TEST_REAL_SIMILAR(trafo.getPairs()[2].first,3.2);
	TEST_REAL_SIMILAR(trafo.getPairs()[0].second,5.2);
	TEST_REAL_SIMILAR(trafo.getPairs()[1].second,6.25);
	TEST_REAL_SIMILAR(trafo.getPairs()[2].second,7.3);
}
END_SECTION

START_SECTION(void store(String filename, const TransformationDescription& transformation))
{
	TransformationDescription trafo,trafo2;
	TransformationXMLFile trafo_xml;

	String tmp_file_1;
	NEW_TMP_FILE(tmp_file_1);
	TEST_EXCEPTION(Exception::IllegalArgument,trafo_xml.store(tmp_file_1,trafo));

	String tmp_file_none;
	trafo.setName("none");
	NEW_TMP_FILE(tmp_file_none);
	trafo_xml.store(tmp_file_none,trafo);
	trafo_xml.load(tmp_file_none,trafo2);
	TEST_STRING_EQUAL(trafo2.getName(),"none");
	TEST_EQUAL(trafo2.getParameters().empty(),true);
	{
    // The actual transformation will be constructed when it is applied for the first time, so let us try this out.
	  DoubleReal pre_image = 234255132.43212;
	  DoubleReal image = pre_image;
	  trafo.apply(image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
	}
	
	String tmp_file_linear;
	NEW_TMP_FILE(tmp_file_linear);
	trafo.clear();
	trafo.setName("linear");
	trafo.setParam("slope",3.141592653589793238);
	trafo.setParam("intercept",2.718281828459045235);
	trafo_xml.store(tmp_file_linear,trafo);
	trafo_xml.load(tmp_file_linear,trafo2);
	TEST_STRING_EQUAL(trafo.getName(),"linear");
	TEST_EQUAL(trafo2.getParameters().size(),2);
	TEST_REAL_SIMILAR(trafo2.getParam("slope"),3.141592653589793238);
	TEST_REAL_SIMILAR(trafo2.getParam("intercept"),2.718281828459045235);
  {
    // The actual transformation will be constructed when it is applied for the first time, so let us try this out.
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = pre_image;
    trafo.apply(image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

	String tmp_file_pairs;
	NEW_TMP_FILE(tmp_file_pairs);
	trafo.clear();
	trafo.setName("pairs");
	TransformationDescription::PairVector pairs;
	pairs.push_back(make_pair(1.2f,5.2f));
	pairs.push_back(make_pair(2.2f,6.25f));
	pairs.push_back(make_pair(3.2f,7.3f));
	trafo.setPairs(pairs);
	trafo_xml.store(tmp_file_pairs,trafo);
	trafo_xml.load(tmp_file_pairs,trafo2);
	TEST_STRING_EQUAL(trafo2.getName(),"pairs");
	TEST_EQUAL(trafo2.getParameters().size(),0);
	TEST_EQUAL(trafo2.getPairs().size(),3);
	TEST_REAL_SIMILAR(trafo2.getPairs()[0].first,1.2);
	TEST_REAL_SIMILAR(trafo2.getPairs()[1].first,2.2);
	TEST_REAL_SIMILAR(trafo2.getPairs()[2].first,3.2);
	TEST_REAL_SIMILAR(trafo2.getPairs()[0].second,5.2);
	TEST_REAL_SIMILAR(trafo2.getPairs()[1].second,6.25);
	TEST_REAL_SIMILAR(trafo2.getPairs()[2].second,7.3);

	trafo.setName("interpolated_linear");
  {
    // The actual transformation will be constructed when it is applied for the first time, so let us try this out.
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = pre_image;
    trafo.apply(image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

  trafo.setName("mumble_pfrwoarpfz");
  {
    // The actual transformation will be constructed when it is applied for the first time, so let us try this out.
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = pre_image;
    TEST_EXCEPTION(Exception::IllegalArgument,trafo.apply(image));
  }

	String tmp_file_bspline;
	NEW_TMP_FILE(tmp_file_bspline);
	trafo.clear();
	trafo.setName("b_spline");
  trafo.setParam("num_breakpoints", 4);
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
  trafo.setPairs(pairs);
  trafo_xml.store(tmp_file_pairs,trafo);
  trafo_xml.load(tmp_file_pairs,trafo2);
  TEST_STRING_EQUAL(trafo2.getName(),"b_spline");
  TEST_EQUAL(trafo2.getParam("num_breakpoints"),4);
  TEST_EQUAL(trafo2.getParameters().size(),1);
  TEST_EQUAL(trafo2.getPairs().size(),11);
  TEST_REAL_SIMILAR(trafo2.getPairs()[0].first,1.2);
  TEST_REAL_SIMILAR(trafo2.getPairs()[0].second,5.2);
  TEST_REAL_SIMILAR(trafo2.getPairs()[10].first,3.7);
  TEST_REAL_SIMILAR(trafo2.getPairs()[10].second,-2.4);
  for ( Int breaks = 0; breaks < 10; ++breaks)
  {
    if ( breaks == 1) continue;
    trafo.setParam("num_breakpoints", breaks);
    // The actual transformation will be constructed when it is applied for the first time, so let us try this out.
    DoubleReal pre_image = 234255132.43212;
    DoubleReal image = pre_image;
    STATUS("breaks: " << breaks);
    trafo.apply(image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
