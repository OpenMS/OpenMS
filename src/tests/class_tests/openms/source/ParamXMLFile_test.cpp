// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/FORMAT/TextFile.h>

///////////////////////////

#include <fstream>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;
using namespace std;

START_TEST(ParamXMLFile, "$Id")

ParamXMLFile* ptr = nullptr;
ParamXMLFile* nullPtr = nullptr;

START_SECTION((ParamXMLFile()))
{
  ptr = new ParamXMLFile();
  TEST_NOT_EQUAL(ptr, nullPtr)
  delete ptr;
}
END_SECTION


START_SECTION((void load(const String& filename, Param& param)))
	Param p2;
  ParamXMLFile paramFile;
	TEST_EXCEPTION(Exception::FileNotFound, paramFile.load("FileDoesNotExist.xml",p2))
END_SECTION

Param p;
p.setValue("test:float",17.4f,"floatdesc");
p.setValue("test:string","test,test,test","stringdesc");
p.setValue("test:int",17,"intdesc");
p.setValue("test2:float",17.5f);
p.setValue("test2:string","test2");
p.setValue("test2:int",18);
p.setSectionDescription("test","sectiondesc");
p.addTags("test:float", {"a","b","c"});

START_SECTION((void store(const String& filename, const Param& param) const))
  ParamXMLFile paramFile;

	Param p2(p);
	p2.setValue("test:a:a1", 47.1,"a1desc\"<>\nnewline");
	p2.setValue("test:b:b1", 47.1);
	p2.setSectionDescription("test:b","bdesc\"<>\nnewline");
	p2.setValue("test2:a:a1", 47.1);
	p2.setValue("test2:b:b1", 47.1,"",{"advanced"});
	p2.setSectionDescription("test2:a","adesc");

	//exception
	Param p300;
	TEST_EXCEPTION(Exception::UnableToCreateFile, paramFile.store("/does/not/exist/FileDoesNotExist.xml",p300))

	String filename;
	NEW_TMP_FILE(filename);
	paramFile.store(filename,p2);
	Param p3;
	paramFile.load(filename,p3);
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), float(p3.getValue("test:float")))
	TEST_EQUAL(p2.getValue("test:string"), p3.getValue("test:string"))
	TEST_EQUAL(p2.getValue("test:int"), p3.getValue("test:int"))
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), float(p3.getValue("test2:float")))
	TEST_EQUAL(p2.getValue("test2:string"), p3.getValue("test2:string"))
	TEST_EQUAL(p2.getValue("test2:int"), p3.getValue("test2:int"))

	TEST_STRING_EQUAL(p2.getDescription("test:float"), p3.getDescription("test:float"))
	TEST_STRING_EQUAL(p2.getDescription("test:string"), p3.getDescription("test:string"))
	TEST_STRING_EQUAL(p2.getDescription("test:int"), p3.getDescription("test:int"))
	TEST_EQUAL(p3.getSectionDescription("test"),"sectiondesc")
	TEST_EQUAL(p3.getDescription("test:a:a1"),"a1desc\"<>\nnewline")
	TEST_EQUAL(p3.getSectionDescription("test:b"),"bdesc\"<>\nnewline")
	TEST_EQUAL(p3.getSectionDescription("test2:a"),"adesc")
	TEST_EQUAL(p3.hasTag("test2:b:b1","advanced"),true)
	TEST_EQUAL(p3.hasTag("test2:a:a1","advanced"),false)
	TEST_EQUAL(ParamXMLFile().isValid(filename, std::cerr),true)

	// advanced
	NEW_TMP_FILE(filename);
	Param p7;
	p7.setValue("true",5,"",{"advanced"});
	p7.setValue("false",5,"");

	paramFile.store(filename,p7);
	TEST_EQUAL(ParamXMLFile().isValid(filename, std::cerr),true)
	Param p8;
	paramFile.load(filename,p8);

	TEST_EQUAL(p8.getEntry("true").tags.count("advanced")==1, true)
	TEST_EQUAL(p8.getEntry("false").tags.count("advanced")==1, false)

	//restrictions
	NEW_TMP_FILE(filename);
	Param p5;
	p5.setValue("int",5);
	p5.setValue("int_min",5);
	p5.setMinInt("int_min",4);
	p5.setValue("int_max",5);
	p5.setMaxInt("int_max",6);
	p5.setValue("int_min_max",5);
	p5.setMinInt("int_min_max",0);
	p5.setMaxInt("int_min_max",10);

	p5.setValue("float",5.1);
	p5.setValue("float_min",5.1);
	p5.setMinFloat("float_min",4.1);
	p5.setValue("float_max",5.1);
	p5.setMaxFloat("float_max",6.1);
	p5.setValue("float_min_max",5.1);
	p5.setMinFloat("float_min_max",0.1);
	p5.setMaxFloat("float_min_max",10.1);

	vector<std::string> strings;
	p5.setValue("string","bli");
	strings.push_back("bla");
	strings.push_back("bluff");
	p5.setValue("string_2","bla");
	p5.setValidStrings("string_2",strings);

		//list restrictions
	vector<std::string> strings2 = {"xml", "txt"};
	p5.setValue("stringlist2",std::vector<std::string>{"a.txt","b.xml","c.pdf"});
	p5.setValue("stringlist",std::vector<std::string>{"aa.C","bb.h","c.doxygen"});
	p5.setValidStrings("stringlist2",strings2);

	p5.setValue("intlist",ListUtils::create<Int>("2,5,10"));
	p5.setValue("intlist2",ListUtils::create<Int>("2,5,10"));
	p5.setValue("intlist3",ListUtils::create<Int>("2,5,10"));
	p5.setValue("intlist4",ListUtils::create<Int>("2,5,10"));
	p5.setMinInt("intlist2",1);
	p5.setMaxInt("intlist3",11);
	p5.setMinInt("intlist4",0);
	p5.setMaxInt("intlist4",15);

	p5.setValue("doublelist",ListUtils::create<double>("1.2,3.33,4.44"));
	p5.setValue("doublelist2",ListUtils::create<double>("1.2,3.33,4.44"));
	p5.setValue("doublelist3",ListUtils::create<double>("1.2,3.33,4.44"));
	p5.setValue("doublelist4",ListUtils::create<double>("1.2,3.33,4.44"));

	p5.setMinFloat("doublelist2",1.1);
	p5.setMaxFloat("doublelist3",4.45);
	p5.setMinFloat("doublelist4",0.1);
	p5.setMaxFloat("doublelist4",5.8);


	paramFile.store(filename,p5);
	TEST_EQUAL(paramFile.isValid(filename, std::cerr),true)
	Param p6;
	paramFile.load(filename,p6);


	TEST_EQUAL(p6.getEntry("int").min_int, -numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("int").max_int, numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("int_min").min_int, 4)
  TEST_EQUAL(p6.getEntry("int_min").max_int, numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("int_max").min_int, -numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("int_max").max_int, 6)
  TEST_EQUAL(p6.getEntry("int_min_max").min_int, 0)
  TEST_EQUAL(p6.getEntry("int_min_max").max_int, 10)

  TEST_REAL_SIMILAR(p6.getEntry("float").min_float, -numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("float").max_float, numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("float_min").min_float, 4.1)
  TEST_REAL_SIMILAR(p6.getEntry("float_min").max_float, numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("float_max").min_float, -numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("float_max").max_float, 6.1)
  TEST_REAL_SIMILAR(p6.getEntry("float_min_max").min_float, 0.1)
  TEST_REAL_SIMILAR(p6.getEntry("float_min_max").max_float, 10.1)

  TEST_EQUAL(p6.getEntry("string").valid_strings.size(), 0)
  TEST_EQUAL(p6.getEntry("string_2").valid_strings.size(), 2)
  TEST_EQUAL(p6.getEntry("string_2").valid_strings[0], "bla")
  TEST_EQUAL(p6.getEntry("string_2").valid_strings[1], "bluff")


  TEST_EQUAL(p6.getEntry("stringlist").valid_strings.size(), 0)
  TEST_EQUAL(p6.getEntry("stringlist2").valid_strings.size(), 2)
  TEST_EQUAL(p6.getEntry("stringlist2").valid_strings[0], "xml")
  TEST_EQUAL(p6.getEntry("stringlist2").valid_strings[1], "txt")

  TEST_EQUAL(p6.getEntry("intlist").min_int, -numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("intlist").max_int, numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("intlist2").min_int, 1)
  TEST_EQUAL(p6.getEntry("intlist2").max_int, numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("intlist3").min_int, -numeric_limits<Int>::max())
  TEST_EQUAL(p6.getEntry("intlist3").max_int, 11)
  TEST_EQUAL(p6.getEntry("intlist4").min_int, 0)
  TEST_EQUAL(p6.getEntry("intlist4").max_int, 15)

  TEST_REAL_SIMILAR(p6.getEntry("doublelist").min_float, -numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("doublelist").max_float, numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("doublelist2").min_float, 1.1)
  TEST_REAL_SIMILAR(p6.getEntry("doublelist2").max_float, numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("doublelist3").min_float, -numeric_limits<double>::max())
  TEST_REAL_SIMILAR(p6.getEntry("doublelist3").max_float, 4.45)
  TEST_REAL_SIMILAR(p6.getEntry("doublelist4").min_float, 0.1)
  TEST_REAL_SIMILAR(p6.getEntry("doublelist4").max_float, 5.8)

  // NaN for float/double
  Param p_nan;
  p_nan.setValue("float_nan", std::numeric_limits<float>::quiet_NaN());
  p_nan.setValue("double_nan", std::numeric_limits<double>::quiet_NaN());
  NEW_TMP_FILE(filename);
  paramFile.store(filename, p_nan);
  Param p_nan2;
  paramFile.load(filename, p_nan2);
  TEST_TRUE(std::isnan((double)p_nan2.getValue("float_nan")))
  TEST_TRUE(std::isnan((double)p_nan2.getValue("double_nan")))
  // ... test the actual values written to INI (should be 'NaN', not 'nan', for compatibility with downstream tools, like Java's double)
  TextFile tf;
  tf.load(filename);
  TEST_TRUE((tf.begin() + 2)->hasSubstring("value=\"NaN\""))
  TEST_TRUE((tf.begin() + 3)->hasSubstring("value=\"NaN\""))

	//Test if an empty Param written to a file validates against the schema
	NEW_TMP_FILE(filename);
	Param p4;
	paramFile.store(filename,p4);
	TEST_EQUAL(paramFile.isValid(filename, std::cerr),true)
END_SECTION

START_SECTION((void writeXMLToStream(std::ostream *os_ptr, const Param &param) const ))
{
	Param p;
	p.setValue("stringlist", std::vector<std::string>{"a","bb","ccc"}, "StringList Description");
	p.setValue("intlist", ListUtils::create<Int>("1,22,333"));
	p.setValue("item", String("bla"));
	p.setValue("stringlist2", std::vector<std::string>());
	p.setValue("intlist2", ListUtils::create<Int>(""));
	p.setValue("item1", 7);
	p.setValue("intlist3", ListUtils::create<Int>("1"));
	p.setValue("stringlist3", std::vector<std::string>{"1"});
	p.setValue("item3", 7.6);
	p.setValue("doublelist", ListUtils::create<double>("1.22,2.33,4.55"));
	p.setValue("doublelist3", ListUtils::create<double>("1.4"));
  p.setValue("file_parameter", "", "This is a file parameter.");
  p.addTag("file_parameter", "input file");
  p.setValidStrings("file_parameter", std::vector<std::string>{"*.mzML","*.mzXML"});
  p.setValue("advanced_parameter", "", "This is an advanced parameter.", {"advanced"});
  p.setValue("flag", "false", "This is a flag i.e. in a command line input it does not need a value.");
  p.setValidStrings("flag",{"true","false"});
  p.setValue("noflagJustTrueFalse", "true", "This is not a flag but has a boolean meaning.");
  p.setValidStrings("noflagJustTrueFalse", {"true","false"});

  String filename;
  NEW_TMP_FILE(filename)
  std::ofstream s(filename.c_str(), std::ios::out);
  ParamXMLFile paramFile;
  paramFile.writeXMLToStream(&s,p);
  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("ParamXMLFile_test_writeXMLToStream.xml"))
}
END_SECTION

START_SECTION([EXTRA] loading and storing of lists)
  ParamXMLFile paramFile;

	Param p;
	p.setValue("stringlist", std::vector<std::string>{"a","bb","ccc"});
	p.setValue("intlist", ListUtils::create<Int>("1,22,333"));
	p.setValue("item", String("bla"));
	p.setValue("stringlist2", std::vector<std::string>());
	p.setValue("intlist2", ListUtils::create<Int>(""));
	p.setValue("item1", 7);
	p.setValue("intlist3", ListUtils::create<Int>("1"));
	p.setValue("stringlist3", std::vector<std::string>{"1"});
	p.setValue("item3", 7.6);
	p.setValue("doublelist", ListUtils::create<double>("1.22,2.33,4.55"));
	p.setValue("doublelist2", ListUtils::create<double>(""));
	p.setValue("doublelist3", ListUtils::create<double>("1.4"));
	//store
	String filename;
	NEW_TMP_FILE(filename);
	paramFile.store(filename,p);
	//load
	Param p2;
	paramFile.load(filename,p2);

	TEST_EQUAL(p2.size(),12);

	TEST_EQUAL(p2.getValue("stringlist").valueType(), ParamValue::STRING_LIST)
	std::vector<std::string> list = p2.getValue("stringlist");
	TEST_EQUAL(list.size(),3)
	TEST_EQUAL(list[0],"a")
	TEST_EQUAL(list[1],"bb")
	TEST_EQUAL(list[2],"ccc")

	TEST_EQUAL(p2.getValue("stringlist2").valueType(), ParamValue::STRING_LIST)
	list = p2.getValue("stringlist2");
	TEST_EQUAL(list.size(),0)

	TEST_EQUAL(p2.getValue("stringlist").valueType(), ParamValue::STRING_LIST)
	list = p2.getValue("stringlist3");
	TEST_EQUAL(list.size(),1)
	TEST_EQUAL(list[0],"1")

	TEST_EQUAL(p2.getValue("intlist").valueType(), ParamValue::INT_LIST)
	IntList intlist = p2.getValue("intlist");
	TEST_EQUAL(intlist.size(),3);
	TEST_EQUAL(intlist[0], 1)
	TEST_EQUAL(intlist[1], 22)
	TEST_EQUAL(intlist[2], 333)

	TEST_EQUAL(p2.getValue("intlist2").valueType(),ParamValue::INT_LIST)
	intlist = p2.getValue("intlist2");
	TEST_EQUAL(intlist.size(),0)

	TEST_EQUAL(p2.getValue("intlist3").valueType(),ParamValue::INT_LIST)
	intlist = p2.getValue("intlist3");
	TEST_EQUAL(intlist.size(),1)
	TEST_EQUAL(intlist[0],1)

	TEST_EQUAL(p2.getValue("doublelist").valueType(), ParamValue::DOUBLE_LIST)
	DoubleList doublelist = p2.getValue("doublelist");
	TEST_EQUAL(doublelist.size(),3);
	TEST_EQUAL(doublelist[0], 1.22)
	TEST_EQUAL(doublelist[1], 2.33)
	TEST_EQUAL(doublelist[2], 4.55)

	TEST_EQUAL(p2.getValue("doublelist2").valueType(),ParamValue::DOUBLE_LIST)
	doublelist = p2.getValue("doublelist2");
	TEST_EQUAL(doublelist.size(),0)

	TEST_EQUAL(p2.getValue("doublelist3").valueType(),ParamValue::DOUBLE_LIST)
	doublelist = p2.getValue("doublelist3");
	TEST_EQUAL(doublelist.size(),1)
	TEST_EQUAL(doublelist[0],1.4)

END_SECTION


START_SECTION(([EXTRA] Escaping of characters))
	Param p;
  ParamXMLFile paramFile;

	p.setValue("string",String("bla"),"string");
	p.setValue("string_with_ampersand", String("bla2&blubb"), "string with ampersand");
	p.setValue("string_with_ampersand_in_descr", String("blaxx"), "String with & in description");
	p.setValue("string_with_single_quote", String("bla'xxx"), "String with single quotes");
	p.setValue("string_with_single_quote_in_descr", String("blaxxx"), "String with ' quote in description");
	p.setValue("string_with_double_quote", String("bla\"xxx"), "String with double quote");
	p.setValue("string_with_double_quote_in_descr", String("bla\"xxx"), "String with \" description");
	p.setValue("string_with_greater_sign", String("bla>xxx"), "String with greater sign");
	p.setValue("string_with_greater_sign_in_descr", String("bla greater xxx"), "String with >");
	p.setValue("string_with_less_sign", String("bla<xxx"), "String with less sign");
	p.setValue("string_with_less_sign_in_descr", String("bla less sign_xxx"), "String with less sign <");


	String filename;
	NEW_TMP_FILE(filename)
	paramFile.store(filename,p);

	Param p2;
	paramFile.load(filename,p2);

	TEST_STRING_EQUAL(p2.getDescription("string"), "string")

  TEST_STRING_EQUAL(p.getValue("string_with_ampersand"), String("bla2&blubb"))
  TEST_STRING_EQUAL(p.getDescription("string_with_ampersand_in_descr"), "String with & in description")
  TEST_STRING_EQUAL(p.getValue("string_with_single_quote"), String("bla'xxx"))
  TEST_STRING_EQUAL(p.getDescription("string_with_single_quote_in_descr"), "String with ' quote in description")
  TEST_STRING_EQUAL(p.getValue("string_with_double_quote"), String("bla\"xxx"))
  TEST_STRING_EQUAL(p.getDescription("string_with_double_quote_in_descr"), "String with \" description")
  TEST_STRING_EQUAL(p.getValue("string_with_greater_sign"), String("bla>xxx"))
  TEST_STRING_EQUAL(p.getDescription("string_with_greater_sign_in_descr"), "String with >")
  TEST_STRING_EQUAL(p.getValue("string_with_less_sign"), String("bla<xxx"))
  TEST_STRING_EQUAL(p.getDescription("string_with_less_sign_in_descr"), "String with less sign <")
END_SECTION

START_SECTION([EXTRA] loading pre 1.6.2 files and storing them in 1.6.2 format)
{
	Param p;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("Param_pre16_update.ini"), p);

  // test some of the former tags if they were loaded correctly
  TEST_EQUAL(p.getValue("SpectraFilterMarkerMower:version"), "1.11.0")
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "input file"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "required"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "advanced"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:out", "output file"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "required"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "advanced"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:log", "advanced"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:log", "required"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:no_progress", "advanced"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:no_progress", "required"), false)

  // write as 1.6.2 ini and check if the output is as expected
	String filename;
	NEW_TMP_FILE(filename)
	paramFile.store(filename,p);

  TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("Param_post16_update.ini"))
}
END_SECTION

START_SECTION([EXTRA] loading 1.6.2 files)
{
	Param p;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("Param_post16_update.ini"), p);

  // test some of the former tags if they were loaded correctly
  TEST_EQUAL(p.getValue("SpectraFilterMarkerMower:version"), "1.11.0")
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "input file"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "required"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "advanced"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:out", "output file"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "required"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:in", "advanced"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:log", "advanced"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:log", "required"), false)

  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:no_progress", "advanced"), true)
  TEST_EQUAL(p.hasTag("SpectraFilterMarkerMower:1:no_progress", "required"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

