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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

//test class with optional parameters
class TOPPBaseTest
  : public TOPPBase
{
  public:
    TOPPBaseTest()
      : TOPPBase("TOPPBaseTest", "A test class", false)
    {
      main(0,0);
    }

    TOPPBaseTest(int argc ,const char** argv)
      : TOPPBase("TOPPBaseTest", "A test class", false)
    {
      main(argc,argv);
    }

    virtual void registerOptionsAndFlags_()
    {
      registerStringOption_("stringoption","<string>","string default","string description",false);
      registerIntOption_("intoption","<int>",4711,"int description",false);
      registerDoubleOption_("doubleoption","<double>",0.4711,"double description",false);
      registerIntList_("intlist","<intlist>",IntList::create("1,2,3,4"),"intlist description",false);
      registerDoubleList_("doublelist","<doublelist>",DoubleList::create("0.4711,1.022,4.0"),"doubelist description",false);
      registerStringList_("stringlist","<stringlist>",StringList::create("abc,def,ghi,jkl"),"stringlist description",false);
      registerFlag_("flag","flag description");

      //for testing write_ini parameter (and with it setDefaults)
      registerStringList_("stringlist2","<stringlist>",StringList::create("hopla,dude"),"stringlist with restrictions",false);
      vector<String> rest;
      rest.push_back("hopla");
      rest.push_back("dude");
      setValidStrings_("stringlist2",rest);

      registerIntList_("intlist2","<int>",IntList::create("3,4,5"),"intlist with restrictions",false);
      setMinInt_("intlist2",2);
      setMaxInt_("intlist2",6);

      registerDoubleList_("doublelist2","<double>",DoubleList::create("1.2,2.33"),"doublelist with restrictions",false);
      setMinFloat_("doublelist2",0.2);
      setMaxFloat_("doublelist2",5.4);

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

    StringList getStringList(const String& name) const
    {
      return getStringList_(name);
    }

    IntList getIntList(const String& name) const
    {
      return getIntList_(name);
    }

    DoubleList getDoubleList(const String& name) const
    {
      return getDoubleList_(name);
    }

    Param const& getParam() const
    {
      return getParam_();
    }

    bool getFlag(const String& name) const
    {
      return getFlag_(name);
    }

    virtual ExitCodes main_(int /*argc*/ , const char** /*argv*/)
    {
      return EXECUTION_OK;
    }

    String const& getIniLocation() const
    {
      return getIniLocation_();
    }

    void inputFileReadable(const String& filename, const String& param_name) const
    {
      inputFileReadable_(filename, param_name);
    }

    void outputFileWritable(const String& filename, const String& param_name) const
    {
      outputFileWritable_(filename, param_name);
    }

    void addDataProcessing(MSExperiment<>& map, DataProcessing::ProcessingAction action)
    {
    	DataProcessing dp = getProcessingInfo_(action);
    	
      addDataProcessing_(map, dp);
      
      //additionally test FeatureMap and ConsensusMap
      FeatureMap<> f_map;
      addDataProcessing_(f_map, dp);
      
      ConsensusMap c_map;
      addDataProcessing_(c_map, dp);
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
      : TOPPBase("TOPPBaseTest", "A test class with non-optional parameters", false)
    {
      main(0,0);
    }

    TOPPBaseTestNOP(int argc , const char** argv)
      : TOPPBase("TOPPBaseTestNOP", "A test class with non-optional parameters", false)
    {
      main(argc,argv);
    }

    virtual void registerOptionsAndFlags_()
    {
      registerStringOption_("stringoption","<string>","string default","string description");
      registerIntOption_("intoption","<int>",4711,"int description");
      registerDoubleOption_("doubleoption","<double>",0.4711,"double description");
      registerFlag_("flag","flag description");
      registerStringList_("stringlist","<stringlist>",StringList::create("abc,def,ghi,jkl"),"stringlist description");
      registerIntList_("intlist","<intlist>",IntList::create("1,2,3,4"),"intlist description");
      registerDoubleList_("doublelist","<doublelist>",DoubleList::create("0.4711,1.022,4.0"),"doubelist description");
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

    StringList getStringList(const String& name) const
    {
      return getStringList_(name);
    }
        IntList getIntList(const String& name) const
    {
      return getIntList_(name);
    }

    DoubleList getDoubleList(const String& name) const
    {
      return getDoubleList_(name);
    }

    virtual ExitCodes main_(int /*argc*/ , const char** /*argv*/)
    {
      return EXECUTION_OK;
    }
};

/////////////////////////////////////////////////////////////

  START_TEST(TOPPBase, "$Id$");

/////////////////////////////////////////////////////////////

TOPPBaseTest* ptr = 0;
START_SECTION((TOPPBase(const String &name, const String &description, bool official=true, bool id_tag_support=false, const String &version="")))
	ptr = new TOPPBaseTest();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~TOPPBase()))
	delete ptr;
END_SECTION

START_SECTION((static Map<String,StringList> getToolList()))
	TEST_EQUAL(TOPPBaseTest::getToolList().has("FileInfo"),true)
	TEST_EQUAL(TOPPBaseTest::getToolList().has("ImaginaryTool"),false)
	TEST_EQUAL(TOPPBaseTest::getToolList()["FileInfo"].empty(),true)
	TEST_EQUAL(TOPPBaseTest::getToolList()["FeatureFinder"].empty(),false)

END_SECTION

START_SECTION((ExitCodes main(int argc, const char**argv)))
	NOT_TESTABLE
	// is tested implicitly in all tests
END_SECTION

//parts to build command lines
const char* a1 ="TOPPBaseTest";
const char* a3 ="-ini";
const char* a5 ="-instance";
const char* a6 ="6";
// needed to get the correct pathes
char* a7;
std::string temp_a7(OPENMS_GET_TEST_DATA_PATH("TOPPBase_toolcommon.ini"));
a7 = new char[temp_a7.size() + 1];
strcpy(a7, temp_a7.c_str());
//
char* a8;
std::string temp_a8(OPENMS_GET_TEST_DATA_PATH("TOPPBase_common.ini"));
a8 = new char[temp_a8.size() + 1];
strcpy(a8, temp_a8.c_str());
//
const char* a9 ="5";
const char* a10 ="-stringoption";
const char* a11 ="-flag";
const char* a12 ="commandline";
const char* a13 ="4.5";
const char* a14 ="-intoption";
const char* a15 ="-doubleoption";
const char* a16 ="4711";
const char* a17 ="-stringlist";
const char* a18 ="-intlist";
const char* a19 ="-doublelist";
const char* a20 ="0.411";
const char* a21 = "-write_ini";
START_SECTION(([EXTRA]String const& getIniLocation_() const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIniLocation(),"TOPPBaseTest:1:")
	//command line
	const char* instance_cl[3] = {a1, a5, a9}; //command line: "TOPPBaseTest -instance 5"
	TOPPBaseTest tmp2(3,instance_cl);
	TEST_EQUAL(tmp2.getIniLocation(),"TOPPBaseTest:5:")
END_SECTION


START_SECTION(([EXTRA]String getStringOption_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getStringOption("stringoption"),"string default");
	//command line
	const char* string_cl[3] = {a1, a10, a12}; //command line: "TOPPBaseTest -stringoption commandline"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getStringOption("stringoption"),"commandline");

	//command line (when there is a ini file value too)
	const char* both_cl[5] = {a1, a10, a12, a3, a7}; //command line: "TOPPBaseTest -stringoption commandline -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp3(5,both_cl);
	TEST_EQUAL(tmp3.getStringOption("stringoption"),DataValue("commandline"));

	//ini file: instance section
	const char* common_cl[3] = {a1, a3, a7}; //command line: "TOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp4(3,common_cl);
	TEST_EQUAL(tmp4.getStringOption("stringoption"),DataValue("instance1"));
	const char* common5_cl[5] = {a1, a3, a7, a5, a9}; //command line: "TOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 5"
	TOPPBaseTest tmp5(5,common5_cl);
	TEST_EQUAL(tmp5.getStringOption("stringoption"),DataValue("instance5"));

	//ini file: tool common section
	const char* common6_cl[5] = {a1, a3, a7, a5, a6}; //command line: "TOPPBaseTest -ini data/TOPPBase_toolcommon.ini -instance 6"
	TOPPBaseTest tmp6(5,common6_cl);
	TEST_EQUAL(tmp6.getStringOption("stringoption"),DataValue("toolcommon"));

	//ini file: common section
	const char* common7_cl[5] = {a1, a3, a8, a5, a6}; //command line: "TOPPBaseTest -ini data/TOPPBase_common.ini -instance 6"
	TOPPBaseTest tmp7(5,common7_cl);
	TEST_EQUAL(tmp7.getStringOption("stringoption"),DataValue("common"));

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getStringOption("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getStringOption("imleeewenit"));

	//missing required parameters
	const char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp8(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp8.getStringOption("stringoption"));

	//test option write_ini
	String filename;
	NEW_TMP_FILE(filename);
	const char* f_name = filename.c_str();
	const char* write_ini[3]={a1, a21, f_name};

	TOPPBaseTest tmp9(3, write_ini);
	Param p1, p2;
	p1.load(filename);
	//remove id pool (the path is dependent on the installation path)
	p1.remove("TOPPBaseTest:1:id_pool");

	//every parameter except for help,ini.instance, write_ini and write_wsdl
	//toolname : TOPPBaseTest
	p2.setValue("TOPPBaseTest:version",VersionInfo::getVersion());
	p2.setValue("TOPPBaseTest:1:stringoption","string default","string description");
	p2.setValue("TOPPBaseTest:1:intoption",4711,"int description");
	p2.setValue("TOPPBaseTest:1:doubleoption",0.4711,"double description");
	p2.setValue("TOPPBaseTest:1:intlist",IntList::create("1,2,3,4"),"intlist description");
	p2.setValue("TOPPBaseTest:1:doublelist",DoubleList::create("0.4711,1.022,4.0"),"doubelist description");
	p2.setValue("TOPPBaseTest:1:stringlist",StringList::create("abc,def,ghi,jkl"),"stringlist description");
	p2.setValue("TOPPBaseTest:1:flag","false","flag description");
	p2.setValue("TOPPBaseTest:1:log","TOPP.log","Location of the log file");
	p2.setValue("TOPPBaseTest:1:debug",0,"Sets the debug level");
	p2.setValue("TOPPBaseTest:1:threads",1, "Sets the number of threads allowed to be used by the TOPP tool");
	p2.setValue("TOPPBaseTest:1:no_progress","false","Disables progress logging to command line");
	p2.setValue("TOPPBaseTest:1:test","false","Enables the test mode (needed for software testing only)");
	//with restriction
  p2.setValue("TOPPBaseTest:1:stringlist2",StringList::create("hopla,dude"),"stringlist with restrictions");
	vector<String> rest;
	rest.push_back("hopla");
	rest.push_back("dude");
	String stringlist2 = "TOPPBaseTest:1:stringlist2";
	p2.setValidStrings(stringlist2,rest);
	String intlist2 = "TOPPBaseTest:1:intlist2";
	String doublelist2 = "TOPPBaseTest:1:doublelist2";
	p2.setValue(intlist2,IntList::create("3,4,5"),"intlist with restriction");
	p2.setMinInt(intlist2,2);
	p2.setMaxInt(intlist2,6);
	p2.setValue(doublelist2,DoubleList::create("1.2,2.33"),"doubelist with restrictions");
	p2.setMinFloat(doublelist2,0.2);
	p2.setMaxFloat(doublelist2,5.4);	
	TEST_EQUAL(p1,p2)
END_SECTION

START_SECTION(([EXTRA]String getIntOption_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIntOption("intoption"),4711);
	//command line
	const char* string_cl[3] = {a1, a14, a6}; //command line: "TOPPBaseTest -intoption 6"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_EQUAL(tmp2.getIntOption("intoption"),6);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getIntOption("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getIntOption("imleeewenit"));

	//missing required parameters
	const char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp3(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp3.getIntOption("intoption"));
END_SECTION

START_SECTION(([EXTRA]String getDoubleOption_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_REAL_SIMILAR(tmp.getDoubleOption("doubleoption"),0.4711);
	//command line
	const char* string_cl[3] = {a1, a15, a13}; //command line: "TOPPBaseTest -doubleoption 4.5"
	TOPPBaseTest tmp2(3,string_cl);
	TEST_REAL_SIMILAR(tmp2.getDoubleOption("doubleoption"),4.5);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getDoubleOption("intoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getDoubleOption("imleeewenit"));

	//missing required parameters
	const char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp3(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp3.getDoubleOption("doubleoption"));
END_SECTION

START_SECTION(([EXTRA] String getIntList_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIntList("intlist"),IntList::create("1,2,3,4"))
	//command line
	const char* string_cl[5]={a1, a18, a6 ,a9 ,a16}; //commandline: "TOPPBaseTest -intlist 6 5 4711"
	TOPPBaseTest tmp2(5, string_cl);
	TEST_EQUAL(tmp2.getIntList("intlist"),IntList::create("6,5,4711"))

	const char* string_cl1[3]={a1, a18, a6}; //commandline: "TOPPBaseTest -intlist 6"
	TOPPBaseTest tmp3(3, string_cl1);
	TEST_EQUAL(tmp3.getIntList("intlist"),IntList::create("6"))

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getIntList("intoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getIntList("imleeewenit"));

	//missing required parameters
	const char* string_cl2[2] = {a1, a11};
	TOPPBaseTestNOP tmp4(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp4.getIntList("intlist"));
END_SECTION

START_SECTION(([EXTRA] String getDoubleList_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getDoubleList("doublelist"),DoubleList::create("0.4711,1.022,4.0"));
	//command line
	const char* string_cl[3]={a1, a19, a20}; //commandline:"TOPPBaseTest -doublelist 0.411"
	TOPPBaseTest tmp2(3, string_cl);
	TEST_EQUAL(tmp2.getDoubleList("doublelist"),DoubleList::create("0.411"));
	const char* a21 = "4.0";
	const char* string_cl2[5]={a1,a19,a20,a13,a21};//commandline :"TOPPBaseTest -doublelist 0.411 4.5 4.0
	TOPPBaseTest tmp3(5,string_cl2);
	TEST_EQUAL(tmp3.getDoubleList("doublelist"),DoubleList::create("0.411,4.5,4.0"));

	const char* string_cl21[4]={a1,a19,a20,a13};//commandline :"TOPPBaseTest -doublelist 0.411 4.5
	TOPPBaseTest tmp31(4,string_cl21);
	TEST_EQUAL(tmp31.getDoubleList("doublelist"),DoubleList::create("0.411,4.5"));

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getDoubleList("intoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getDoubleList("imleeewenit"));
	//missing required parameters
	const char* string_cl3[2] = {a1, a11};
	TOPPBaseTestNOP tmp4(2,string_cl3);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp4.getDoubleList("doublelist"));
END_SECTION

START_SECTION(([EXTRA] String getStringList_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getStringList("stringlist"),StringList::create("abc,def,ghi,jkl"));
	//command line
	const char* string_cl[3]={a1,a17,a12};	//commandline: "TOPPBaseTest -stringlist conmandline"
	TOPPBaseTest tmp2(3, string_cl);
	TEST_EQUAL(tmp2.getStringList("stringlist"),StringList::create("commandline"))

	const char* string_cl2[5]={a1,a17,a12,a7, a8};	//commandline: "TOPPBaseTest -stringlist conmandline data/TOPPBase_toolcommon.ini data/TOPPBase_common.ini"
	TOPPBaseTest tmp3(5, string_cl2);
	StringList tmp_stringlist;
	tmp_stringlist << "commandline";
	tmp_stringlist << OPENMS_GET_TEST_DATA_PATH("TOPPBase_toolcommon.ini");
	tmp_stringlist << OPENMS_GET_TEST_DATA_PATH("TOPPBase_common.ini");
	TEST_EQUAL(tmp3.getStringList("stringlist"),tmp_stringlist);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getStringList("intoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getStringList("imleeewenit"));

	//missing required parameters
	const char* string_cl3[2] = {a1, a11};
	TOPPBaseTestNOP tmp4(2,string_cl3);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp4.getStringList("stringlist"));

END_SECTION

START_SECTION(([EXTRA]bool getFlag_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getFlag("flag"),false);
	//command line
	const char* flag_cl[2] = {a1, a11}; //command line: "TOPPBaseTest -flag"
	TOPPBaseTest tmp2(2,flag_cl);
	TEST_EQUAL(tmp2.getFlag("flag"),true);

	TEST_EXCEPTION(Exception::WrongParameterType,tmp2.getFlag("doubleoption"));
	TEST_EXCEPTION(Exception::UnregisteredParameter,tmp2.getFlag("imleeewenit"));
END_SECTION

START_SECTION(([EXTRA]void inputFileReadable_(const String& filename, const String& param_name) const))
	TOPPBaseTest tmp;
	TEST_EXCEPTION(Exception::FileNotFound,tmp.inputFileReadable("/this/file/does/not/exist.txt","someparam"));
	TEST_EXCEPTION(Exception::FileEmpty,tmp.inputFileReadable(OPENMS_GET_TEST_DATA_PATH("TOPPBase_empty.txt"), "someparam"));
	tmp.inputFileReadable(OPENMS_GET_TEST_DATA_PATH("TOPPBase_common.ini"),"ini");
END_SECTION

START_SECTION(([EXTRA]void outputFileWritable_(const String& filename, const String& param_name) const))
	TEST_EXCEPTION(Exception::UnableToCreateFile,TOPPBaseTest().outputFileWritable("/this/file/cannot/be/written/does_not_exists.txt","someparam"));

	String filename;
	NEW_TMP_FILE(filename);
	TOPPBaseTest().outputFileWritable(filename, "");
	//Actually writing something to the file is not necessary, but on Mac all tmp files are called 'source_<line>.tmp'.
	//So we have to make sure the file is empty. Otherwise the test might fail...
	TextFile dummy;
	dummy.resize(5);
	dummy.store(filename);
END_SECTION

START_SECTION(([EXTRA]void parseRange_(const String& text, double& low, double& high) const))
	TOPPBaseTest topp;
	double a = -1.0;
	double b = -1.0;
	String s;

	s = ":";
	topp.parseRange(s,a,b);
	TEST_REAL_SIMILAR(a,-1.0);
	TEST_REAL_SIMILAR(b,-1.0);

	s = "4.5:";
	topp.parseRange(s,a,b);
	TEST_REAL_SIMILAR(a,4.5);
	TEST_REAL_SIMILAR(b,-1.0);

	s = ":5.5";
	topp.parseRange(s,a,b);
	TEST_REAL_SIMILAR(a,4.5);
	TEST_REAL_SIMILAR(b,5.5);

	s = "6.5:7.5";
	topp.parseRange(s,a,b);
	TEST_REAL_SIMILAR(a,6.5);
	TEST_REAL_SIMILAR(b,7.5);
END_SECTION

START_SECTION(([EXTRA]Param getParam_( const std::string& prefix ) const))
{
	//ini file
	const char* tmp_argv[] = {a1, a3, a7}; //command line: "TOPPBaseTest -ini data/TOPPBase_toolcommon.ini"
	TOPPBaseTest tmp_topp(sizeof(tmp_argv)/sizeof(*tmp_argv),tmp_argv);

	Param good_params;
	good_params.setValue( "TOPPBaseTest:stringoption", "toolcommon" );
	good_params.setValue( "ini", OPENMS_GET_TEST_DATA_PATH("TOPPBase_toolcommon.ini") );
	good_params.setValue( "stringoption", "instance1" );

	TEST_EQUAL(tmp_topp.getParam(), good_params);
}
END_SECTION

START_SECTION(([EXTRA] data processing methods))
	MSExperiment<> exp;
	exp.resize(2);

	TOPPBaseTest topp;
	topp.addDataProcessing(exp, DataProcessing::ALIGNMENT);
	
	for (Size i=0; i<exp.size(); ++i)
	{
		TEST_EQUAL(exp[i].getDataProcessing().size(),1)
		TEST_EQUAL(exp[i].getDataProcessing()[0].getSoftware().getName(),"TOPPBaseTest")
		TEST_NOT_EQUAL(exp[i].getDataProcessing()[0].getSoftware().getVersion(),"1.1")
		TEST_EQUAL(exp[i].getDataProcessing()[0].getCompletionTime().isValid(),true)
		TEST_EQUAL(exp[i].getDataProcessing()[0].getProcessingActions().size(),1)
		TEST_EQUAL(*(exp[i].getDataProcessing()[0].getProcessingActions().begin()),DataProcessing::ALIGNMENT)
	}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



