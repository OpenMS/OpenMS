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
#include <OpenMS/APPLICATIONS/TOPPBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TOPPBase, "$Id$")

class Test
	: public TOPPBase
{
	public:
		Test()
			: TOPPBase("TOPPTest")
		{
			main(0,0);
		}
		
		Test(int argc , char** argv)
			: TOPPBase("TOPPTest")
		{
			main(argc,argv);	
		}
		
		virtual void printToolUsage_() const
		{
			
		}
		
		virtual void printToolHelpOpt_() const
		{
			
		}	
	
		virtual void setOptionsAndFlags_()
		{
			options_["-option"] = "option_internal";
			flags_["-flag"] = "flag_internal";
		}

		String getParamAsString(const String& key, const String& default_value="") const
		{
			return getParamAsString_(key,default_value);
		}
		
		SignedInt getParamAsInt(const String& key, SignedInt default_value=0) const
		{
			return getParamAsInt_(key,default_value);
		}
		
		double getParamAsDouble(const String& key, double default_value=0) const
		{
			return getParamAsDouble_(key,default_value);
		}
		
		bool getParamAsBool(const String& key, bool default_value=false) const
		{
			return getParamAsBool_(key,default_value);
		}
		
		DataValue const& getParam(const String& key) const
		{
			return getParam_(key);
		}
		
		Param const& getParam() const
		{
			return getParam_();
		}
		
		Param getParamCopy( const std::string& prefix ) const
		{
			return getParamCopy_(prefix);
		}
			
		virtual ExitCodes main_(int /*argc*/ , char** /*argv*/)
		{
			return EXECUTION_OK;
		}
};

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Test* ptr = 0;
CHECK(TOPPBase())
	ptr = new Test();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~TOPPBase())
	delete ptr;
RESULT

//parts to build command lines
char* a1 ="TOPPTest";
char* a3 ="-ini";
char* a5 ="-n";
char* a6 ="6";
char* a7 ="data/TOPPBase_toolcommon.ini";
char* a8 ="data/TOPPBase_common.ini";
char* a9 ="5";
char* a10 ="-option";
char* a11 ="-flag";
char* a12 ="commandline";
char* a13 ="4.5";


CHECK(String const& getToolName() const)
	Test tmp;
	TEST_EQUAL(tmp.getToolName(),"TOPPTest")
RESULT

CHECK(String const& getIniLocation() const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getIniLocation(),"TOPPTest:1:")
	//command line
	char* instance_cl[3] = {a1, a5, a9}; //command line: "TOPPTest -n 5"
	Test tmp2(3,instance_cl);
	TEST_EQUAL(tmp2.getIniLocation(),"TOPPTest:5:")
RESULT

CHECK(bool getParamAsBool_(const String& key, bool default_value=false) const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getParamAsBool("flag_internal",false),false);
	TEST_EQUAL(tmp.getParamAsBool("flag_internal",true),true);
	//command line
	char* flag_cl[2] = {a1, a11}; //command line: "TOPPTest -flag"
	Test tmp2(2,flag_cl);
	TEST_EQUAL(tmp2.getParamAsBool("flag_internal",false),true);
RESULT

CHECK(SignedInt getParamAsInt_(const String& key, SignedInt default_value=0) const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getParamAsInt("option_internal",17),17);
	//command line
	char* int_cl[3] = {a1, a10, a9}; //command line: "TOPPTest -option 5"
	Test tmp2(3,int_cl);
	TEST_EQUAL(tmp2.getParamAsInt("option_internal",17),5);
RESULT

CHECK(double getParamAsDouble_(const String& key, double default_value=0) const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getParamAsDouble("option_internal",17.0),17.0);
	//command line
	char* double_cl[3] = {a1, a10, a13}; //command line: "TOPPTest -option 4.5"
	Test tmp2(3,double_cl);
	TEST_EQUAL(tmp2.getParamAsDouble("option_internal",17.0),4.5);
RESULT

CHECK(String getParamAsString_(const String& key, const String& default_value="") const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getParamAsString("option_internal","leer"),"leer");
	//command line
	char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPTest -option commandline"
	Test tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getParamAsString("option_internal","leer"),"commandline");
RESULT

CHECK(DataValue const& getParam_(const String& key) const)
	//default 
	Test tmp;
	TEST_EQUAL(tmp.getParam("option_internal"),DataValue::EMPTY);
	
	//command line
	char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPTest -option commandline"
	Test tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getParam("option_internal"),DataValue("commandline"));
	
	//command line (when there is a ini file value too)
	char* both_cl[5] = {a1, a10, a12, a3, a7}; //command line: "TOPPTest -option commandline -ini data/TOPPBase_toolcommon.ini"
	Test tmp3(5,both_cl);
	TEST_EQUAL(tmp3.getParam("option_internal"),DataValue("commandline"));
	
	//ini file: instance section
	char* common_cl[3] = {a1, a3, a7}; //command line: "TOPPTest -ini data/TOPPBase_toolcommon.ini"
	Test tmp4(3,common_cl);
	TEST_EQUAL(tmp4.getParam("option_internal"),DataValue("instance1"));
	char* common5_cl[5] = {a1, a3, a7, a5, a9}; //command line: "TOPPTest -ini data/TOPPBase_toolcommon.ini -n 5"
	Test tmp5(5,common5_cl);
	TEST_EQUAL(tmp5.getParam("option_internal"),DataValue("instance5"));
	
	//ini file: tool common section
	char* common6_cl[5] = {a1, a3, a7, a5, a6}; //command line: "TOPPTest -ini data/TOPPBase_toolcommon.ini -n 6"
	Test tmp6(5,common6_cl);
	TEST_EQUAL(tmp6.getParam("option_internal"),DataValue("toolcommon"));

	//ini file: common section
	char* common7_cl[5] = {a1, a3, a8, a5, a6}; //command line: "TOPPTest -ini data/TOPPBase_common.ini -n 6"
	Test tmp7(5,common7_cl);
	TEST_EQUAL(tmp7.getParam("option_internal"),DataValue("common"));

	//ini file: inheritence
	//TODO
RESULT

CHECK(Param const& getParam_() const)
	//TODO
RESULT

CHECK(Param getParamCopy_( const std::string& prefix ) const)
	//TODO
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



