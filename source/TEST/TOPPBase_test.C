// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/TOPPBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TOPPBase, "$Id$");

//test class with optional parameters
class TOPPBaseTest
	: public TOPPBase
{
	public:
		TOPPBaseTest()
			: TOPPBase("TOPPBaseTest", "A test class")
		{
			main(0,0);
		}

		TOPPBaseTest(int argc , char** argv)
			: TOPPBase("TOPPBaseTest", "A test class")
		{
			main(argc,argv);
		}

		virtual void registerOptionsAndFlags_()
		{
			registerStringOption_("stringoption","<string>","string default","string description",false);
			registerIntOption_("intoption","<int>",4711,"int description",false);
			registerDoubleOption_("doubleoption","<double>",0.4711,"double description",false);
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

		Int getIntOption(const String& name) const
		{
			return getIntOption_(name);
		}

		Param const& getParam() const
		{
			return getParam_();
		}

		bool getFlag(const String& name) const
		{
			return getFlag_(name);
		}

		bool setByUser(const String& name) const
		{
			return setByUser_(name);
		}

		virtual ExitCodes main_(int /*argc*/ , char** /*argv*/)
		{
			return EXECUTION_OK;
		}

		String const& getIniLocation() const
		{
			return getIniLocation_();
		}

		void inputFileReadable(const String& filename) const
		{
			inputFileReadable_(filename);
		}

		void outputFileWritable(const String& filename) const
		{
			outputFileWritable_(filename);
		}

		void parseRange(const String& text, double& low, double& high) const
		{
			parseRange_(text, low, high);
		}

};

// Test class for no-optional parameters
class TOPPBaseTestNOP
	: public TOPPBase
{
	public:
		TOPPBaseTestNOP()
			: TOPPBase("TOPPBaseTest", "A test class with non-optional parameters")
		{
			main(0,0);
		}

		TOPPBaseTestNOP(int argc , char** argv)
			: TOPPBase("TOPPBaseTestNOP", "A test class with non-optional parameters")
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

		Int getIntOption(const String& name) const
		{
			return getIntOption_(name);
		}

		virtual ExitCodes main_(int /*argc*/ , char** /*argv*/)
		{
			return EXECUTION_OK;
		}
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOPPBaseTest* ptr = 0;
CHECK((TOPPBase(const String& tool_name, const String& tool_description)))
	ptr = new TOPPBaseTest();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~TOPPBase()))
	delete ptr;
RESULT

CHECK(ExitCodes main(int argc, char **argv))
	// is tested implicitly in all tests
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
char* a16 ="4711";

CHECK(([EXTRA]String const& getIniLocation_() const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIniLocation(),"TOPPBaseTest:1:")
	//command line
	char* instance_cl[3] = {a1, a5, a9}; //command line: "TOPPTOPPBaseTest -instance 5"
	TOPPBaseTest tmp2(3,instance_cl);
	TEST_EQUAL(tmp2.getIniLocation(),"TOPPBaseTest:5:")
RESULT

CHECK([EXTRA] bool setByUser_(const String& name) const)
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.setByUser("intoption"),false);

	//command line
	char* string_cl[3] = {a1, a14, a16}; //command line: "TOPPTOPPBaseTest -intoption 4711"
	TOPPBaseTest tmp2(3,string_cl);

	TEST_EQUAL(tmp2.setByUser("intoption"),true);
	TEST_EQUAL(tmp2.setByUser("stringoption"),false);
	TEST_EQUAL(tmp2.setByUser("doubleoption"),false);

	//ini file
	char* both_cl[3] = {a1, a3, a7}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp3(3,both_cl);

	TEST_EQUAL(tmp3.setByUser("intoption"),false);
	TEST_EQUAL(tmp3.setByUser("stringoption"),true);
	TEST_EQUAL(tmp3.setByUser("doubleoption"),false);
RESULT

CHECK(([EXTRA]String getStringOption_(const String& name) const))
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

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getStringOption("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getStringOption("imleeewenit"));

	//missing required parameters
	char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp8(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp8.getStringOption("stringoption"));
RESULT

CHECK(([EXTRA]String getIntOption_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIntOption("intoption"),4711);
	//command line
	char* string_cl[3] = {a1, a14, a6}; //command line: "TOPPTOPPBaseTest -intoption 6"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getIntOption("intoption"),6);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getIntOption("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getIntOption("imleeewenit"));

	//missing required parameters
	char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp3(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp3.getIntOption("intoption"));
RESULT

CHECK(([EXTRA]String getDoubleOption_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_REAL_EQUAL(tmp.getDoubleOption("doubleoption"),0.4711);
	//command line
	char* string_cl[3] = {a1, a15, a13}; //command line: "TOPPTOPPBaseTest -doubleoption 4.5"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_REAL_EQUAL(tmp2.getDoubleOption("doubleoption"),4.5);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getDoubleOption("intoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getDoubleOption("imleeewenit"));

	//missing required parameters
	char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp3(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp3.getDoubleOption("doubleoption"));
RESULT

CHECK(([EXTRA]bool getFlag_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getFlag("flag"),false);
	//command line
	char* flag_cl[2] = {a1, a11}; //command line: "TOPPTOPPBaseTest -flag"
	TOPPBaseTest tmp2(2,flag_cl);
	TEST_EQUAL(tmp2.getFlag("flag"),true);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getFlag("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getFlag("imleeewenit"));
RESULT

CHECK(([EXTRA]void inputFileReadable_(const String& filename) const))
	TOPPBaseTest tmp;
	TEST_EXCEPTION(Exception::FileNotFound,tmp.inputFileReadable("/this/file/does/not/exist.txt"));
	TEST_EXCEPTION(Exception::FileEmpty,tmp.inputFileReadable("data/TOPPBase_empty.txt"));
	tmp.inputFileReadable("data/TOPPBase_common.ini");
RESULT

CHECK(([EXTRA]void outputFileWritable_(const String& filename) const))
	TEST_EXCEPTION(Exception::UnableToCreateFile,TOPPBaseTest().outputFileWritable("/this/file/cannot/be/written/does_not_exists.txt"));

	String filename;
	NEW_TMP_FILE(filename);
	TOPPBaseTest().outputFileWritable(filename);
RESULT

CHECK(([EXTRA]void parseRange_(const String& text, double& low, double& high) const))
	TOPPBaseTest topp;
	double a = -1.0;
	double b = -1.0;
	String s;

	s = ":";
	topp.parseRange(s,a,b);
	TEST_REAL_EQUAL(a,-1.0);
	TEST_REAL_EQUAL(b,-1.0);

	s = "4.5:";
	topp.parseRange(s,a,b);
	TEST_REAL_EQUAL(a,4.5);
	TEST_REAL_EQUAL(b,-1.0);

	s = ":5.5";
	topp.parseRange(s,a,b);
	TEST_REAL_EQUAL(a,4.5);
	TEST_REAL_EQUAL(b,5.5);

	s = "6.5:7.5";
	topp.parseRange(s,a,b);
	TEST_REAL_EQUAL(a,6.5);
	TEST_REAL_EQUAL(b,7.5);
RESULT

CHECK(([EXTRA]Param getParam_( const std::string& prefix ) const))
{
	//ini file
	char* tmp_argv[] = {a1, a3, a7}; //command line: "TOPPTOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp_topp(sizeof(tmp_argv)/sizeof(*tmp_argv),tmp_argv);

	Param good_params;
	good_params.setValue( "TOPPBaseTest:stringoption", "toolcommon" );
	good_params.setValue( "ini", "data/TOPPBase_toolcommon.ini" );
	good_params.setValue( "stringoption", "instance1" );

	TEST_EQUAL(tmp_topp.getParam(), good_params);
}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



