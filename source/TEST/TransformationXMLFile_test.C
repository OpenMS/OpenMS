// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

TransformationXMLFile* ptr;
CHECK((TransformationXMLFile()))
{
	ptr = new TransformationXMLFile();
	TEST_NOT_EQUAL(ptr,0);
}
RESULT

CHECK([EXTRA] static bool isValid(const String& filename))
{
	TransformationXMLFile f;
	TEST_EQUAL(f.isValid("data/TransformationXMLFile_1.xml"),true);	
	TEST_EQUAL(f.isValid("data/TransformationXMLFile_2.xml"),true);	
	TEST_EQUAL(f.isValid("data/TransformationXMLFile_3.xml"),false);	
	TEST_EQUAL(f.isValid("data/TransformationXMLFile_4.xml"),true);	
}
RESULT


CHECK(void load(const String& filename, TransformationDescription& transformation))
{
	TransformationDescription trafo;
	TransformationXMLFile trafo_xml;

	trafo_xml.load("data/TransformationXMLFile_1.xml",trafo);
	TEST_STRING_EQUAL(trafo.getName(),"none");
	TEST_EQUAL(trafo.getParameters().empty(),true);

	trafo_xml.load("data/TransformationXMLFile_2.xml",trafo);
	TEST_STRING_EQUAL(trafo.getName(),"linear");
	TEST_EQUAL(trafo.getParameters().size(),2);
	TEST_REAL_EQUAL(trafo.getParam("slope"),3.141592653589793238);
	TEST_REAL_EQUAL(trafo.getParam("intercept"),2.718281828459045235);

	trafo_xml.load("data/TransformationXMLFile_4.xml",trafo);
	TEST_STRING_EQUAL(trafo.getName(),"pairs");
	TEST_EQUAL(trafo.getParameters().size(),0);
	TEST_EQUAL(trafo.getPairs().size(),3);
	TEST_REAL_EQUAL(trafo.getPairs()[0].first,1.2);
	TEST_REAL_EQUAL(trafo.getPairs()[1].first,2.2);
	TEST_REAL_EQUAL(trafo.getPairs()[2].first,3.2);
	TEST_REAL_EQUAL(trafo.getPairs()[0].second,5.2);
	TEST_REAL_EQUAL(trafo.getPairs()[1].second,6.25);
	TEST_REAL_EQUAL(trafo.getPairs()[2].second,7.3);
}
RESULT

CHECK(void store(String filename, const TransformationDescription& transformation))
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
	TEST_REAL_EQUAL(trafo2.getParam("slope"),3.141592653589793238);
	TEST_REAL_EQUAL(trafo2.getParam("intercept"),2.718281828459045235);

	String tmp_file_pairs;
	NEW_TMP_FILE(tmp_file_pairs);
	trafo.clear();
	trafo.setName("pairs");
	TransformationDescription::PairVector pairs;
	pairs.push_back(make_pair(1.2,5.2));
	pairs.push_back(make_pair(2.2,6.25));
	pairs.push_back(make_pair(3.2,7.3));
	trafo.setPairs(pairs);
	trafo_xml.store(tmp_file_pairs,trafo);
	trafo_xml.load(tmp_file_pairs,trafo2);
	TEST_STRING_EQUAL(trafo2.getName(),"pairs");
	TEST_EQUAL(trafo2.getParameters().size(),0);
	TEST_EQUAL(trafo2.getPairs().size(),3);
	TEST_REAL_EQUAL(trafo2.getPairs()[0].first,1.2);
	TEST_REAL_EQUAL(trafo2.getPairs()[1].first,2.2);
	TEST_REAL_EQUAL(trafo2.getPairs()[2].first,3.2);
	TEST_REAL_EQUAL(trafo2.getPairs()[0].second,5.2);
	TEST_REAL_EQUAL(trafo2.getPairs()[1].second,6.25);
	TEST_REAL_EQUAL(trafo2.getPairs()[2].second,7.3);
}
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
