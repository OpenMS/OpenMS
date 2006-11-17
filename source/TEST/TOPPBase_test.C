// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm, Clemens Gröpl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/TOPPBase2.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TOPPBase, "$Id$")

class TOPPBaseTest
	: public TOPPBase2
{
	public:
		TOPPBaseTest()
			: TOPPBase2("TOPPBaseTest", "A test class")
		{
			main(0,0);
		}
		
		TOPPBaseTest(int argc , char** argv)
			: TOPPBase2("TOPPBaseTest", "A test class")
		{
			main(argc,argv);	
		}
		
		virtual void registerOptionsAndFlags_()
		{
			registerStringOption_("stringoption","<string>","string default","string description");
			registerIntOption_("intoption","<int>",4711,"int description");
			registerDoubleOption_("doubleoption","<double>",0.4711,"double description");
			registerFlag_("flag","flag description");
		}

		String getStringOption(const String& name) const
		{
			return getStringOption_(name);
		}

		double getDoubleOption(const String& name) const
		{
			return getDoubleOption_(name);
		}

		SignedInt getIntOption(const String& name) const
		{
			return getIntOption_(name);
		}

		bool getFlag(const String& name) const
		{
			return getFlag_(name);
		}
			
		virtual ExitCodes main_(int /*argc*/ , char** /*argv*/)
		{
			return EXECUTION_OK;
		}

// COMMENTED OUT AS THESE METHODS ARE private!
// DO NOT DELETE AS THESE TESTS CAN STILL BE USED WHEN CHANGING TOPPBase!
// SIMPLY MAKE THESE METHODS protected FOR A WHILE.
//		String const& getIniLocation() const
//		{
//			return getIniLocation_();
//		}
//
//		String getParamAsString(const String& key, const String& default_value="") const
//		{
//			return getParamAsString_(key,default_value);
//		}
//		
//		SignedInt getParamAsInt(const String& key, SignedInt default_value=0) const
//		{
//			return getParamAsInt_(key,default_value);
//		}
//		
//		double getParamAsDouble(const String& key, double default_value=0) const
//		{
//			return getParamAsDouble_(key,default_value);
//		}
//		
//		bool getParamAsBool(const String& key, bool default_value=false) const
//		{
//			return getParamAsBool_(key,default_value);
//		}
//		
//		DataValue const& getParam(const String& key) const
//		{
//			return getParam_(key);
//		}
//		
//		Param const& getParam() const
//		{
//			return getParam_();
//		}
//		
//		Param getParamCopy( const std::string& prefix ) const
//		{
//			return getParamCopy_(prefix);
//		}
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOPPBaseTest* ptr = 0;
CHECK(TOPPBase())
	ptr = new TOPPBaseTest();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~TOPPBase())
	delete ptr;
RESULT

//parts to build command lines
char* a1 ="TOPPBaseTest";
char* a3 ="-ini";
char* a5 ="-instance";
char* a6 ="6";
char* a7 ="data/TOPPBase_toolcommon.ini";
char* a8 ="data/TOPPBase_common.ini";
char* a9 ="5";
char* a10 ="-stringoption";
char* a11 ="-flag";
char* a12 ="commandline";
char* a13 ="4.5";
char* a14 ="-intoption";
char* a15 ="-doubleoption";

// COMMENTED OUT AS THESE METHODS ARE private!
// DO NOT DELETE AS THESE TESTS CAN STILL BE USED WHEN CHANGING TOPPBase!
// SIMPLY MAKE THESE METHODS protected FOR A WHILE.
//
//CHECK(String const& getIniLocation_() const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_EQUAL(tmp.getIniLocation(),"TOPPBaseTest:1:")
//	//command line
//	char* instance_cl[3] = {a1, a5, a9}; //command line: "TOPPTOPPBaseTest -instance 5"
//	TOPPBaseTest tmp2(3,instance_cl);
//	TEST_EQUAL(tmp2.getIniLocation(),"TOPPBaseTest:5:")
//RESULT
//
//CHECK(bool getParamAsBool_(const String& key, bool default_value=false) const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_EQUAL(tmp.getParamAsBool("flag",false),false);
//	TEST_EQUAL(tmp.getParamAsBool("flag",true),true);
//	//command line
//	char* flag_cl[2] = {a1, a11}; //command line: "TOPPTOPPBaseTest -flag"
//	TOPPBaseTest tmp2(2,flag_cl);
//	TEST_EQUAL(tmp2.getParamAsBool("flag",false),true);
//RESULT
//
//CHECK(SignedInt getParamAsInt_(const String& key, SignedInt default_value=0) const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_EQUAL(tmp.getParamAsInt("stringoption",17),17);
//	//command line
//	char* int_cl[3] = {a1, a10, a9}; //command line: "TOPPTOPPBaseTest -stringoption 5"
//	TOPPBaseTest tmp2(3,int_cl);
//	TEST_EQUAL(tmp2.getParamAsInt("stringoption",17),5);
//RESULT
//
//CHECK(double getParamAsDouble_(const String& key, double default_value=0) const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_REAL_EQUAL(tmp.getParamAsDouble("stringoption",17.0),17.0);
//	//command line
//	char* double_cl[3] = {a1, a10, a13}; //command line: "TOPPTOPPBaseTest -stringoption 4.5"
//	TOPPBaseTest tmp2(3,double_cl);
//	TEST_REAL_EQUAL(tmp2.getParamAsDouble("stringoption",17.0),4.5);
//RESULT
//
//CHECK(String getParamAsString_(const String& key, const String& default_value="") const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_EQUAL(tmp.getParamAsString("stringoption","leer"),"leer");
//	//command line
//	char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPTOPPBaseTest -stringoption commandline"
//	TOPPBaseTest tmp2(3,string_cl);
//	TEST_EQUAL(tmp2.getParamAsString("stringoption","leer"),"commandline");
//RESULT
//
//CHECK(DataValue const& getParam_(const String& key) const)
//	//default 
//	TOPPBaseTest tmp;
//	TEST_EQUAL(tmp.getParam("stringoption"),DataValue::EMPTY);
//	
//	//command line
//	char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPTOPPBaseTest -stringoption commandline"
//	TOPPBaseTest tmp2(3,string_cl);
//	TEST_EQUAL(tmp2.getParam("stringoption"),DataValue("commandline"));
//	
//	//command line (when there is a ini file value too)
//	char* both_cl[5] = {a1, a10, a12, a3, a7}; //command line: "TOPPTOPPBaseTest -stringoption commandline -ini data/TOPPBase_toolcommon.ini"
//	TOPPBaseTest tmp3(5,both_cl);
//	TEST_EQUAL(tmp3.getParam("stringoption"),DataValue("commandline"));
//	
//	//ini file: instance section
//	char* common_cl[3] = {a1, a3, a7}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
//	TOPPBaseTest tmp4(3,common_cl);
//	TEST_EQUAL(tmp4.getParam("stringoption"),DataValue("instance1"));
//	char* common5_cl[5] = {a1, a3, a7, a5, a9}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 5"
//	TOPPBaseTest tmp5(5,common5_cl);
//	TEST_EQUAL(tmp5.getParam("stringoption"),DataValue("instance5"));
//	
//	//ini file: tool common section
//	char* common6_cl[5] = {a1, a3, a7, a5, a6}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 6"
//	TOPPBaseTest tmp6(5,common6_cl);
//	TEST_EQUAL(tmp6.getParam("stringoption"),DataValue("toolcommon"));
//
//	//ini file: common section
//	char* common7_cl[5] = {a1, a3, a8, a5, a6}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_common.ini -instance 6"
//	TOPPBaseTest tmp7(5,common7_cl);
//	TEST_EQUAL(tmp7.getParam("stringoption"),DataValue("common"));
//
//	//ini file: inheritence
//	//TODO
//RESULT
//
//CHECK(Param const& getParam_() const)
//	//TODO
//RESULT
//
//CHECK(Param getParamCopy_( const std::string& prefix ) const)
//	//TODO
//RESULT


CHECK(String getStringOption_(const String& name) const)
	//default 
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getStringOption("stringoption"),"string default");
	//command line
	char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPTOPPBaseTest -stringoption commandline"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getStringOption("stringoption"),"commandline");
	
	//command line (when there is a ini file value too)
	char* both_cl[5] = {a1, a10, a12, a3, a7}; //command line: "TOPPTOPPBaseTest -stringoption commandline -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp3(5,both_cl);
	TEST_EQUAL(tmp3.getStringOption("stringoption"),DataValue("commandline"));
	
	//ini file: instance section
	char* common_cl[3] = {a1, a3, a7}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp4(3,common_cl);
	TEST_EQUAL(tmp4.getStringOption("stringoption"),DataValue("instance1"));
	char* common5_cl[5] = {a1, a3, a7, a5, a9}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 5"
	TOPPBaseTest tmp5(5,common5_cl);
	TEST_EQUAL(tmp5.getStringOption("stringoption"),DataValue("instance5"));
	
	//ini file: tool common section
	char* common6_cl[5] = {a1, a3, a7, a5, a6}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 6"
	TOPPBaseTest tmp6(5,common6_cl);
	TEST_EQUAL(tmp6.getStringOption("stringoption"),DataValue("toolcommon"));

	//ini file: common section
	char* common7_cl[5] = {a1, a3, a8, a5, a6}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_common.ini -instance 6"
	TOPPBaseTest tmp7(5,common7_cl);
	TEST_EQUAL(tmp7.getStringOption("stringoption"),DataValue("common"));
RESULT

CHECK(String getIntOption_(const String& name) const)
	//default 
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIntOption("intoption"),4711);
	//command line
	char* string_cl[3] = {a1, a14, a6}; //command line: "TOPPTOPPBaseTest -intoption 6"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getIntOption("intoption"),6);
RESULT

CHECK(String getDoubleOption_(const String& name) const)
	//default 
	TOPPBaseTest tmp;
	TEST_REAL_EQUAL(tmp.getDoubleOption("doubleoption"),0.4711);
	//command line
	char* string_cl[3] = {a1, a15, a13}; //command line: "TOPPTOPPBaseTest -doubleoption 4.5"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_REAL_EQUAL(tmp2.getDoubleOption("doubleoption"),4.5);
RESULT

CHECK(bool getFlag_(const String& name) const)
	//default 
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getFlag("flag"),false);
	//command line
	char* flag_cl[2] = {a1, a11}; //command line: "TOPPTOPPBaseTest -flag"
	TOPPBaseTest tmp2(2,flag_cl);
	TEST_EQUAL(tmp2.getFlag("flag"),true);
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



