// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <cstdlib>

#include <QStringList>
///////////////////////////

using namespace OpenMS;
using namespace std;

//test class with optional parameters
class TOPPBaseTest
  : public TOPPBase
{
  public:
    TOPPBaseTest()
      : TOPPBase("TOPPBaseTest", "A test class", false, {}, false)
    {
      char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      main(0,nullptr);
    }

    TOPPBaseTest(int argc ,const char** argv)
      : TOPPBase("TOPPBaseTest", "A test class", false, {}, false)
    {
      char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      main(argc,argv);
    }

    void registerOptionsAndFlags_() override
    {
      registerStringOption_("stringoption","<string>","string default","string description",false);
      registerIntOption_("intoption","<int>",4711,"int description",false);
      registerDoubleOption_("doubleoption","<double>",0.4711,"double description",false);
      registerIntList_("intlist","<intlist>",ListUtils::create<Int>("1,2,3,4"),"intlist description",false);
      registerDoubleList_("doublelist","<doublelist>", ListUtils::create<double>("0.4711,1.022,4.0"),"doubelist description",false);
      registerStringList_("stringlist","<stringlist>", ListUtils::create<String>("abc,def,ghi,jkl"),"stringlist description",false);
      registerFlag_("flag","flag description");

      //for testing write_ini parameter (and with it setDefaults)
      registerStringList_("stringlist2","<stringlist>", ListUtils::create<String>("hopla,dude"),"stringlist with restrictions",false);
      vector<String> rest;
      rest.push_back("hopla");
      rest.push_back("dude");
      setValidStrings_("stringlist2",rest);

      registerIntList_("intlist2","<int>",ListUtils::create<Int>("3,4,5"),"intlist with restrictions",false);
      setMinInt_("intlist2",2);
      setMaxInt_("intlist2",6);

      registerDoubleList_("doublelist2","<double>", ListUtils::create<double>("1.2,2.33"),"doublelist with restrictions",false);
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

    ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
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

    void addDataProcessing(PeakMap& map, DataProcessing::ProcessingAction action)
    {
    	DataProcessing dp = getProcessingInfo_(action);

      addDataProcessing_(map, dp);

      //additionally test FeatureMap and ConsensusMap
      FeatureMap f_map;
      addDataProcessing_(f_map, dp);

      ConsensusMap c_map;
      addDataProcessing_(c_map, dp);
    }

    bool parseRange(const String& text, double& low, double& high) const
    {
      return parseRange_(text, low, high);
    }

    TOPPBase::ExitCodes runExternalProcess(const QString& executable, const QStringList& arguments, const QString& workdir) const
    {
      return runExternalProcess_(executable, arguments, workdir);
    }

};

// Test class for no-optional parameters
class TOPPBaseTestNOP
  : public TOPPBase
{
  public:
    TOPPBaseTestNOP()
      : TOPPBase("TOPPBaseTestNOP", "A test class with non-optional parameters", false, {}, false)
    {
      char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      main(0,nullptr);
    }

    TOPPBaseTestNOP(int argc , const char** argv)
      : TOPPBase("TOPPBaseTestNOP", "A test class with non-optional parameters", false, {}, false)
    {
      char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      main(argc,argv);
    }

    void registerOptionsAndFlags_() override
    {
      registerStringOption_("stringoption","<string>","","string description");
      registerIntOption_("intoption","<int>",0,"int description",false);
      registerDoubleOption_("doubleoption","<double>", -1.0,"double description", false);
      registerFlag_("flag","flag description");
      registerStringList_("stringlist","<stringlist>", ListUtils::create<String>(""),"stringlist description");
      registerIntList_("intlist","<intlist>",ListUtils::create<Int>(""),"intlist description");
      registerDoubleList_("doublelist","<doublelist>", ListUtils::create<double>(""),"doubelist description");
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

    ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
    {
      return EXECUTION_OK;
    }
};

// Test class for parameters derived from a Param object
class TOPPBaseTestParam: public TOPPBase
{
  public:
    TOPPBaseTestParam(const Param& param):
			TOPPBase("TOPPBaseTestParam", "A test class with parameters derived from Param", false, {}, false), test_param_(param)
    {
      static char* var = (char *)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      main(0, nullptr);
    }

    void registerOptionsAndFlags_() override
    {
      registerFullParam_(test_param_);
    }

    ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
    {
      return EXECUTION_OK;
    }

    const Param& getParam() const
    {
      return getParam_();
    }

  private:
    Param test_param_;
};

//test class with optional parameters
class TOPPBaseCmdParseTest
  : public TOPPBase
{

public:
  TOPPBaseCmdParseTest()
    : TOPPBase("TOPPBaseCmdParseTest", "A test class to test parts of the cmd parser functionality", false, {}, false)
  {}

  void registerOptionsAndFlags_() override
  {
  }

  ExitCodes run(int argc , const char** argv)
  {
    static char* var = (char *)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
    return main(argc, argv);
  }

  ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
  {
    static char* var = (char *)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
    return EXECUTION_OK;
  }
};

//test class with optional parameters
class TOPPBaseCmdParseSubsectionsTest
: public TOPPBase
{

public:
  TOPPBaseCmdParseSubsectionsTest()
  : TOPPBase("TOPPBaseCmdParseSubsectionsTest", "A test class to test parts of the cmd parser functionality", false, {}, false)
  {}

  void registerOptionsAndFlags_() override
  {
    registerStringOption_("stringoption","<string>","","string description");
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("other", "Other parameters section");
  }

  Param getSubsectionDefaults_(const String & section) const override
  {
    Param p;
    if (section == "algorithm")
    {
      p.setValue("param1", "param1_value", "param1_description");
      p.setValue("param2", "param2_value", "param2_description");
    }
    else // "other"
    {
      p.setValue("param3", "param3_value", "param3_description");
      p.setValue("param4", "param4_value", "param4_description");
      p.setValue("flagparam", "false", "this will be a flag");
      p.setValidStrings("flagparam", {"true","false"});
      p.setValue("nonflagparam", "true", "this will be a string param with true/false");
      p.setValidStrings("nonflagparam", {"true","false"});
    }
    return p;
  }

  ExitCodes run(int argc , const char** argv)
  {
    static char* var = (char *)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
    return main(argc, argv);
  }

  ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
  {
    return EXECUTION_OK;
  }

  String getStringOption(String name)
  {
    return getStringOption_(name);
  }

  Param getParam()
  {
    return getParam_();
  }
};

/////////////////////////////////////////////////////////////

  START_TEST(TOPPBase, "$Id$");

/////////////////////////////////////////////////////////////

TOPPBaseTest* ptr = nullptr;
TOPPBaseTest* nullPointer = nullptr;
START_SECTION(TOPPBase(const String& name, const String& description, bool official = true, const std::vector<Citation>& citations = {}, bool toolhandler_test = true))
	ptr = new TOPPBaseTest();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~TOPPBase()))
	delete ptr;
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
const char* test = "-test";

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
	const char* string_cl2[2] = {a1, a11}; // TOPPBaseTest -flag
	TOPPBaseTestNOP tmp8(2,string_cl2);
	TEST_EXCEPTION(Exception::RequiredParameterNotGiven,tmp8.getStringOption("stringoption"));

	//test option write_ini
	String filename;
	NEW_TMP_FILE(filename);
	const char* f_name = filename.c_str();
	const char* write_ini[3]={a1, a21, f_name};

	TOPPBaseTest tmp9(3, write_ini);
	Param p1, p2;
  ParamXMLFile paramFile;
	paramFile.load(filename, p1);
	//remove id pool (the path is dependent on the installation path)
	p1.remove("TOPPBaseTest:1:id_pool");

	//every parameter except for help,ini.instance, write_ini and write_wsdl
	//toolname : TOPPBaseTest
	p2.setValue("TOPPBaseTest:version",VersionInfo::getVersion());
	p2.setValue("TOPPBaseTest:1:stringoption","string default","string description");
	p2.setValue("TOPPBaseTest:1:intoption",4711,"int description");
	p2.setValue("TOPPBaseTest:1:doubleoption",0.4711,"double description");
	p2.setValue("TOPPBaseTest:1:intlist",ListUtils::create<Int>("1,2,3,4"),"intlist description");
	p2.setValue("TOPPBaseTest:1:doublelist", ListUtils::create<double>("0.4711,1.022,4.0"),"doubelist description");
	p2.setValue("TOPPBaseTest:1:stringlist", std::vector<std::string>{"abc","def","ghi","jkl"},"stringlist description");
	p2.setValue("TOPPBaseTest:1:flag","false","flag description");
  p2.setValue("TOPPBaseTest:1:log","","Name of log file (created only when specified)");
	p2.setValue("TOPPBaseTest:1:debug",0,"Sets the debug level");
	p2.setValue("TOPPBaseTest:1:threads",1, "Sets the number of threads allowed to be used by the TOPP tool");
	p2.setValue("TOPPBaseTest:1:no_progress","false","Disables progress logging to command line");
	p2.setValue("TOPPBaseTest:1:force","false","Overwrite tool specific checks.");
	p2.setValue("TOPPBaseTest:1:test","false","Enables the test mode (needed for software testing only)");
	//with restriction
  p2.setValue("TOPPBaseTest:1:stringlist2", std::vector<std::string>{"hopla","dude"},"stringlist with restrictions");
	vector<std::string> rest = {"hopla","dude"};
	String stringlist2 = "TOPPBaseTest:1:stringlist2";
	p2.setValidStrings(stringlist2,rest);
	String intlist2 = "TOPPBaseTest:1:intlist2";
	String doublelist2 = "TOPPBaseTest:1:doublelist2";
	p2.setValue(intlist2,ListUtils::create<Int>("3,4,5"),"intlist with restriction");
	p2.setMinInt(intlist2,2);
	p2.setMaxInt(intlist2,6);
	p2.setValue(doublelist2, ListUtils::create<double>("1.2,2.33"),"doubelist with restrictions");
	p2.setMinFloat(doublelist2,0.2);
	p2.setMaxFloat(doublelist2,5.4);
	TEST_EQUAL(p1,p2)
	WHITELIST("version")
	TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("TOPPBase_test_write_ini_out.ini"))

  String filename2;
  NEW_TMP_FILE(filename2);
  const char* f_name2 = filename2.c_str();
  const char* b = "TOPPBaseCmdParseSubsectionsTest";
  const char* write_ini2[3]={b, a21, f_name2};
  TOPPBaseCmdParseSubsectionsTest tmp10{};
	tmp10.run(3, write_ini2);
  TEST_FILE_SIMILAR(filename2, OPENMS_GET_TEST_DATA_PATH("TOPPBase_test_write_ini_subsec_out.ini"))
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
	//-> not testable, as ints cannot be made 'required' (no NAN supported)
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
END_SECTION

START_SECTION(([EXTRA] String getIntList_(const String& name) const))
	//default
	TOPPBaseTest tmp;
	TEST_EQUAL(tmp.getIntList("intlist") == ListUtils::create<Int>("1,2,3,4"), true)
	//command line
	const char* string_cl[5]={a1, a18, a6 ,a9 ,a16}; //commandline: "TOPPBaseTest -intlist 6 5 4711"
	TOPPBaseTest tmp2(5, string_cl);
	TEST_EQUAL(tmp2.getIntList("intlist") == ListUtils::create<Int>("6,5,4711"), true)

	const char* string_cl1[3]={a1, a18, a6}; //commandline: "TOPPBaseTest -intlist 6"
	TOPPBaseTest tmp3(3, string_cl1);
	TEST_EQUAL(tmp3.getIntList("intlist") == ListUtils::create<Int>("6"), true)

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
	TEST_EQUAL(tmp.getDoubleList("doublelist") == ListUtils::create<double>("0.4711,1.022,4.0"), true)
	//command line
	const char* string_cl[3]={a1, a19, a20}; //commandline:"TOPPBaseTest -doublelist 0.411"
	TOPPBaseTest tmp2(3, string_cl);
	TEST_EQUAL(tmp2.getDoubleList("doublelist") == ListUtils::create<double>("0.411"), true)
	const char* a210 = "4.0";
	const char* string_cl2[5]={a1, a19, a20, a13, a210};//commandline :"TOPPBaseTest -doublelist 0.411 4.5 4.0
	TOPPBaseTest tmp3(5,string_cl2);
	TEST_EQUAL(tmp3.getDoubleList("doublelist") == ListUtils::create<double>("0.411,4.5,4.0"), true)

	const char* string_cl21[4]={a1,a19,a20,a13};//commandline :"TOPPBaseTest -doublelist 0.411 4.5
	TOPPBaseTest tmp31(4,string_cl21);
	TEST_EQUAL(tmp31.getDoubleList("doublelist") == ListUtils::create<double>("0.411,4.5"), true)

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
	TEST_EQUAL(tmp.getStringList("stringlist") == ListUtils::create<String>("abc,def,ghi,jkl"), true)
	//command line
	const char* string_cl[3]={a1,a17,a12};	//commandline: "TOPPBaseTest -stringlist conmandline"
	TOPPBaseTest tmp2(3, string_cl);
	TEST_EQUAL(tmp2.getStringList("stringlist") == ListUtils::create<String>("commandline"), true)

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
	TEST_EXCEPTION(Exception::FileNotFound, tmp.inputFileReadable("/this/file/does/not/exist.txt", "someparam"));
	TEST_EXCEPTION(Exception::FileEmpty, tmp.inputFileReadable(OPENMS_GET_TEST_DATA_PATH("TOPPBase_empty.txt"), "someparam"));
	tmp.inputFileReadable(OPENMS_GET_TEST_DATA_PATH("TOPPBase_common.ini"), "ini");
END_SECTION

START_SECTION(([EXTRA]void outputFileWritable_(const String& filename, const String& param_name) const))
	TEST_EXCEPTION(Exception::UnableToCreateFile,TOPPBaseTest().outputFileWritable("/this/file/cannot/be/written/does_not_exists.txt","someparam"));

	String filename;
	NEW_TMP_FILE(filename);
	TOPPBaseTest().outputFileWritable(filename, "");
	//Actually writing something to the file is not necessary, but on Mac all tmp files are called 'source_<line>.tmp'.
	//So we have to make sure the file is empty. Otherwise the test might fail...
	TextFile dummy;
  dummy.addLine("");dummy.addLine("");dummy.addLine("");dummy.addLine("");dummy.addLine("");
	dummy.store(filename);
END_SECTION

START_SECTION(([EXTRA]void parseRange_(const String& text, double& low, double& high) const))
{
	TOPPBaseTest topp;
	double a = -1.0;
	double b = -1.0;

	String s = ":";
	bool result = topp.parseRange(s, a, b);
	TEST_REAL_SIMILAR(a, -1.0);
	TEST_REAL_SIMILAR(b, -1.0);
  TEST_EQUAL(result, false);

	s = "4.5:";
	result = topp.parseRange(s, a, b);
	TEST_REAL_SIMILAR(a, 4.5);
	TEST_REAL_SIMILAR(b, -1.0);
  TEST_EQUAL(result, true);

	s = ":5.5";
	result = topp.parseRange(s, a, b);
	TEST_REAL_SIMILAR(a, 4.5);
	TEST_REAL_SIMILAR(b, 5.5);
  TEST_EQUAL(result, true);

	s = "6.5:7.5";
	result = topp.parseRange(s, a, b);
	TEST_REAL_SIMILAR(a, 6.5);
	TEST_REAL_SIMILAR(b, 7.5);
  TEST_EQUAL(result, true);
}
END_SECTION

START_SECTION(([EXTRA] TOPPBase::ExitCodes TOPPBase::runExternalProcess_(const QString& executable, const QStringList& arguments, const QString& workdir) const))
{

// we just need ANY commandline tool available on (hopefully) all boxes.
// note that commands like "dir" or "type" are only known within cmd.exe and are not actual executables (unlike on Linux)
#ifdef OPENMS_WINDOWSPLATFORM
  const QString exe = "cmd";
  const QStringList args = QStringList() << "/C" << "echo hi";
  const QStringList args_broken = QStringList() << "/C" << "doesnotexist";
#else
  const QString exe = "ls";
  const QStringList args("-l");
  const QStringList args_broken = QStringList() << "-0";
#endif // 

  TOPPBaseTest topp;
  auto result = topp.runExternalProcess("/path/does/not/exists.exe", QStringList(), "");
  TEST_EQUAL(result, TOPPBase::EXTERNAL_PROGRAM_NOTFOUND);

  result = topp.runExternalProcess(exe, args_broken, "");
  TEST_EQUAL(result, TOPPBase::EXTERNAL_PROGRAM_ERROR);
  
  result = topp.runExternalProcess(exe, args, "");
  TEST_EQUAL(result, TOPPBase::EXECUTION_OK);
}
END_SECTION

START_SECTION(([EXTRA] data processing methods))
	PeakMap exp;
	exp.resize(2);

	TOPPBaseTest topp;
	topp.addDataProcessing(exp, DataProcessing::ALIGNMENT);

	for (Size i=0; i<exp.size(); ++i)
	{
		TEST_EQUAL(exp[i].getDataProcessing().size(),1)
		TEST_EQUAL(exp[i].getDataProcessing()[0]->getSoftware().getName(),"TOPPBaseTest")
		TEST_NOT_EQUAL(exp[i].getDataProcessing()[0]->getSoftware().getVersion(),"1.1")
		TEST_EQUAL(exp[i].getDataProcessing()[0]->getCompletionTime().isValid(),true)
		TEST_EQUAL(exp[i].getDataProcessing()[0]->getProcessingActions().size(),1)
		TEST_EQUAL(*(exp[i].getDataProcessing()[0]->getProcessingActions().begin()),DataProcessing::ALIGNMENT)
	}
END_SECTION

START_SECTION(([EXTRA] const Param& getParam_()))
{
	Param test_param;
	test_param.setValue("param_int", 123, "param int description");
	test_param.setValue("param_double", -4.56, "param double description");
	test_param.setValue("param_string", "test", "param string description");
	test_param.setValue("param_stringlist", std::vector<std::string>{"this","is","a","test"}, "param stringlist description");
	test_param.setValue("param_intlist", ListUtils::create<Int>("7,-8,9"), "param intlist description");
	test_param.setValue("param_doublelist", ListUtils::create<double>("123,-4.56,0.789"), "param doublelist description");
	test_param.setValue("param_flag", "true", "param flag description");
	test_param.setValidStrings("param_flag", {"true","false"});

	TOPPBaseTestParam temp(test_param);
	Param result = temp.getParam(); // contains "test_param" + some default stuff
	for (Param::ParamIterator it = test_param.begin(); it != test_param.end(); ++it)
	{
		TEST_EQUAL(*it == result.getEntry(it.getName()), true);
	}
}
END_SECTION

START_SECTION((static void setMaxNumberOfThreads(int num_threads)))
{
  // this is a helper function that is only working if openmp is active
  // due to bugs in the different OpenMP implementations it is not realy
  // testable
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([EXTRA] misc options on command line))
{
  // misc text option
  const char* string_cl[3] = {a1, a12, test}; //command line: "TOPPBaseTest commandline"
  TOPPBaseCmdParseTest tmp1;
  TOPPBase::ExitCodes ec1 = tmp1.run(3,string_cl);
  TEST_EQUAL(ec1, TOPPBase::ILLEGAL_PARAMETERS)

  // unknown option
  TOPPBaseCmdParseTest tmp2;
  const char* string_cl_2[4] = {a1, a10, a12, test}; //command line: "TOPPBaseTest -stringoption commandline"
  TOPPBase::ExitCodes ec2 = tmp1.run(4,string_cl_2);
  TEST_EQUAL(ec2, TOPPBase::ILLEGAL_PARAMETERS)
}
END_SECTION

const char* a22 = "-algorithm:param1";
const char* a23 = "-algorithm:param2";
const char* a24 = "-other:param3";
const char* a25 = "-other:param4";
const char* a26 = "val1";
const char* a27 = "val2";
const char* a28 = "val3";
const char* a29 = "val4";
std::string temp_a30(OPENMS_GET_TEST_DATA_PATH("TOPPBaseCmdParseSubsectionsTest.ini"));
const char* a30 = temp_a30.c_str();

START_SECTION(([EXTRA] test subsection parameters))
{
  const char* string_cl_1[4] = {a1, a10, a12, test}; //command line: "TOPPBaseTest -stringoption commandline"
  TOPPBaseCmdParseSubsectionsTest tmp1;
  TOPPBase::ExitCodes ec1 = tmp1.run(4, string_cl_1);
  TEST_EQUAL(ec1, TOPPBase::EXECUTION_OK)
  TEST_EQUAL(tmp1.getStringOption("stringoption"), "commandline");
  TEST_EQUAL(tmp1.getParam().getValue("algorithm:param1"), "param1_value");
  TEST_EQUAL(tmp1.getParam().getValue("algorithm:param2"), "param2_value");
  TEST_EQUAL(tmp1.getParam().getValue("other:param3"), "param3_value");
  TEST_EQUAL(tmp1.getParam().getValue("other:param4"), "param4_value");

  // overwrite from cmd
  const char* string_cl_2[12] = {a1, a10, a12, a22, a26, a23, a27, a24, a28, a25, a29, test}; //command line: "TOPPBaseTest -algorithm:param1 val1 -algorithm:param2 val2 -algorithm:param3 val3 -algorithm:param4 val4 -stringoption commandline"
  TOPPBaseCmdParseSubsectionsTest tmp2;
  TOPPBase::ExitCodes ec2 = tmp2.run(12, string_cl_2);
  TEST_EQUAL(ec2, TOPPBase::EXECUTION_OK)
  TEST_EQUAL(tmp2.getStringOption("stringoption"), "commandline");
  TEST_EQUAL(tmp2.getParam().getValue("algorithm:param1"), "val1");
  TEST_EQUAL(tmp2.getParam().getValue("algorithm:param2"), "val2");
  TEST_EQUAL(tmp2.getParam().getValue("other:param3"), "val3");
  TEST_EQUAL(tmp2.getParam().getValue("other:param4"), "val4");

  // overwrite ini values from cmd
  const char* string_cl_3[10] = {a1, a3, a30, a22, a26, a25, a29, a10, a12, test }; //command line: "TOPPBaseTest -ini TOPPBaseCmdParseSubsectionsTest.ini -algorithm:param1 val1 -algorithm:param4 val4"
  TOPPBaseCmdParseSubsectionsTest tmp3;
  TOPPBase::ExitCodes ec3 = tmp3.run(10, string_cl_3);
  TEST_EQUAL(ec3, TOPPBase::EXECUTION_OK)
  TEST_EQUAL(tmp3.getStringOption("stringoption"), "commandline");
  TEST_EQUAL(tmp3.getParam().getValue("algorithm:param1"), "val1");
  TEST_EQUAL(tmp3.getParam().getValue("algorithm:param2"), "param2_ini_value");
  TEST_EQUAL(tmp3.getParam().getValue("other:param3"), "param3_ini_value");
  TEST_EQUAL(tmp3.getParam().getValue("other:param4"), "val4");
}
END_SECTION

delete [] a7;
delete [] a8;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



