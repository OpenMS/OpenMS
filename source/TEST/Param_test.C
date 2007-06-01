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

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Param* d10_ptr = 0;
CHECK((Param()))
	d10_ptr = new Param();
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~Param()))
	delete d10_ptr;
RESULT

CHECK((const DataValue& getValue(const String& key) const))
	Param p;
	TEST_EQUAL(p.getValue("key")==DataValue::EMPTY, true)
RESULT

CHECK((const String& getDescription(const String& key) const))
	Param p;
	TEST_EQUAL(p.getDescription("key"),String::EMPTY)
RESULT

CHECK((void setValue(const String& key, const String& value, String description)))
	Param p;
	p.setValue("key","value");
	TEST_EQUAL(p.getValue("key"), "value")
	TEST_STRING_EQUAL(p.getDescription("key"), String::EMPTY)
	p.setValue("key","value","thisvalue");
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
	p.setValue("key","value");
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
RESULT

CHECK((void setValue(const String& key, Int value, String description)))
	Param p;
	p.setValue("key",17);
	TEST_EQUAL(Int(p.getValue("key")), 17)
	TEST_STRING_EQUAL(p.getDescription("key"), String::EMPTY)
	p.setValue("key",17,"thisvalue");
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
	p.setValue("key",17);
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
RESULT

CHECK((void setValue(const String& key, float value, String description)))
	Param p;
	p.setValue("key",17.4f);
	TEST_REAL_EQUAL(float(p.getValue("key")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("key"), String::EMPTY)
	p.setValue("key",17.4f,"thisvalue");
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
	p.setValue("key",17.4f);
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
RESULT

CHECK((void setValue(const String& key, double value, String description)))
	Param p;
	p.setValue("key",17.4);
	TEST_REAL_EQUAL(double(p.getValue("key")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("key"), String::EMPTY)
	p.setValue("key",17.4,"thisvalue");
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
	p.setValue("key",17.4);
	TEST_STRING_EQUAL(p.getDescription("key"), "thisvalue")
RESULT

CHECK((bool empty() const))
	Param p;
	TEST_EQUAL(p.empty(), true)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.empty(), false)
RESULT

CHECK((void clear()))
	Param p;
	p.setValue("key",17.4,"keydesc");
	p.clear();
	TEST_EQUAL(p.empty(), true)
	TEST_EQUAL(p.getDescription("key"),"")
	
RESULT

CHECK((UInt size() const))
	Param p;
	TEST_EQUAL(p.size(), 0)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.size(), 1)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.size(), 1)
RESULT

CHECK((void setDescription(const String& location, const String& description)))
	Param p;
	//no value -> no desctription
	p.setDescription("test","sectiondesc");
	p.setDescription("test:float","valuedesc");
	TEST_EQUAL(p.getDescription("test"), String::EMPTY)
	TEST_EQUAL(p.getDescription("test:float"), String::EMPTY)
	
	//value -> description	
	p.setValue("test:float",47.1);
	p.setDescription("test","sectiondesc");
	p.setDescription("test:float","valuedesc");
	TEST_EQUAL(p.getDescription("test"), "sectiondesc")
	TEST_EQUAL(p.getDescription("test:float"), "valuedesc")
RESULT


Param p;
p.setValue("test:float",17.4f,"floatdesc");
p.setValue("test:string","test,test,test","stringdesc");
p.setValue("test:int",17,"intdesc");
p.setValue("test2:float",17.5f);
p.setValue("test2:string","test2");
p.setValue("test2:int",18);
p.setDescription("test","sectiondesc");

CHECK((Param(const Param& rhs)))
	Param p2(p);
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p.getDescription("test:int"), "intdesc")
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")
RESULT

CHECK((Param& operator = (const Param& rhs)))
	Param p2;
	p2=p;
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")
RESULT


CHECK((void remove(const String& prefix)))
	Param p2(p);
	
	p2.remove("test:float");
	TEST_EQUAL(p2.getValue("test:float"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")

	p2.remove("test:");
	TEST_EQUAL(p2.getValue("test:float"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test:string"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test:int"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), String::EMPTY)
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")

	p2.remove("test");
	TEST_EQUAL(p2.getValue("test:float"), DataValue::EMPTY)
	TEST_EQUAL(p2.getValue("test:string"), DataValue::EMPTY)
	TEST_EQUAL(p2.getValue("test:int"), DataValue::EMPTY)
	TEST_EQUAL(p2.getValue("test2:float"), DataValue::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), DataValue::EMPTY)
	TEST_EQUAL(p2.getValue("test2:int"), DataValue::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"")
	TEST_STRING_EQUAL(p2.getDescription("test:float"), String::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), String::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), String::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), String::EMPTY)
RESULT


CHECK((bool operator == (const Param& rhs) const))
	Param p2(p);
	TEST_EQUAL(p==p2, true)
	p2.setValue("test:float",17.5f);
	TEST_EQUAL(p==p2, false)
	p2 = p;
	p2.setValue("test:float3",17.4f);
	TEST_EQUAL(p==p2, false)
	p2 = p;
	p2.remove("test:float");
	TEST_EQUAL(p==p2, false)
RESULT

CHECK((void load(const String& filename) throw(Exception::FileNotFound, Exception::ParseError)))
	Param p2;
	TEST_EXCEPTION(Exception::FileNotFound, p2.load("FileDoesNotExist.xml"))	
RESULT

CHECK((void store(const String& filename) const throw(Exception::UnableToCreateFile)))
	Param p2(p);
	p2.setValue("test:a:a1", 47.1,"a1desc");
	p2.setValue("test:b:b1", 47.1);
	p2.setDescription("test:b","bdesc");
	p2.setValue("test2:a:a1", 47.1);
	p2.setValue("test2:b:b1", 47.1);
	p2.setDescription("test2:a","adesc");
	
	//exception
	Param p300;
	TEST_EXCEPTION(Exception::UnableToCreateFile, p300.store("/does/not/exist/FileDoesNotExist.xml"))	
	
	String filename;
	NEW_TMP_FILE(filename);
	p2.store(filename);
	Param p3;
	p3.load(filename);
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), float(p3.getValue("test:float")))
	TEST_EQUAL(p2.getValue("test:string"), p3.getValue("test:string"))
	TEST_EQUAL(p2.getValue("test:int"), p3.getValue("test:int"))
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), float(p3.getValue("test2:float")))
	TEST_EQUAL(p2.getValue("test2:string"), p3.getValue("test2:string"))
	TEST_EQUAL(p2.getValue("test2:int"), p3.getValue("test2:int"))	
	
	TEST_STRING_EQUAL(p2.getDescription("test:float"), p3.getDescription("test:float"))
	TEST_STRING_EQUAL(p2.getDescription("test:string"), p3.getDescription("test:string"))
	TEST_STRING_EQUAL(p2.getDescription("test:int"), p3.getDescription("test:int"))
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), p3.getDescription(""))
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), p3.getDescription(""))
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), p3.getDescription(""))
	TEST_EQUAL(p3.getDescription("test"),"sectiondesc")
	TEST_EQUAL(p3.getDescription("test:a"),"")
	TEST_EQUAL(p3.getDescription("test:a:a1"),"a1desc")
	TEST_EQUAL(p3.getDescription("test:b"),"bdesc")
	TEST_EQUAL(p3.getDescription("test:b:b1"),"")
	TEST_EQUAL(p3.getDescription("test2:a"),"adesc")
RESULT

CHECK((void insert(String prefix, const Param& para)))
	Param p2;
	p2.insert("test3",p);
	TEST_REAL_EQUAL(float(p2.getValue("test3:test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test3:test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test3:test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test3:test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test3:test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test3:test:int"), "intdesc")
	TEST_REAL_EQUAL(float(p2.getValue("test3:test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test3:test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test3:test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test3:test"),"sectiondesc")
		
	p2.insert("",p);
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)	
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")
RESULT

CHECK((Param copy(const String& prefix, bool remove_prefix=false, String new_prefix="") const))
	Param p2;

	p2 = p.copy("notthere:");
	TEST_EQUAL((p2==Param()),true)

	p2 = p.copy("test:");
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(p2.getValue("test2:float"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:int"), DataValue::EMPTY)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"")

	p2 = p.copy("test:",true);
	TEST_REAL_EQUAL(float(p2.getValue("float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("float"), "floatdesc")
	TEST_EQUAL(p2.getValue("string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("string"), "stringdesc")
	TEST_EQUAL(p2.getDescription("test"),"")
	
	p2 = p.copy("test:",true,"tttest");
	TEST_REAL_EQUAL(float(p2.getValue("tttest:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("tttest:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("tttest:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("tttest:string"), "stringdesc")
	TEST_EQUAL(p2.getDescription("tttest:test"),"")

	p2 = p.copy("test:",false,"tttest");
	TEST_REAL_EQUAL(float(p2.getValue("tttest:test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("tttest:test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("tttest:test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("tttest:test:string"), "stringdesc")
	TEST_EQUAL(p2.getDescription("tttest:test"),"")
	
	p2 = p.copy("test");
	TEST_REAL_EQUAL(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_EQUAL(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getDescription("test"),"sectiondesc")
RESULT

CHECK((Param copyWithInherit(const String& old_prefix, const String& new_prefix="") const))
{
	Param p0;
	p0.setValue("test:float",17.4f,"test:float");
	p0.setValue("test:inherit","test2","test:inherit");
	p0.setValue("test:int",17,"test:int");
	p0.setValue("test:string","test,test,test","test:string");

	p0.setValue("test2:double",18.2,"test2:double");
	p0.setValue("test2:float",17.5f,"test2:float");
	p0.setValue("test2:inherit","test3:test3a","test2:inherit");
	p0.setValue("test2:string","test2","test2:string");

	p0.setValue("test3:bla","wrong","test3:bla");
	p0.setValue("test3:test3a:anotherint",99,"test3:test3a:anotherint");
	p0.setValue("test3:test3a:bla","blubber","test3:test3a:bla");
	p0.setValue("test3:test3a:inherit","non-existent:location","test3:test3a:inherit");

	Param p2;

	p2 = p0.copyWithInherit("notthere:");
	TEST_EQUAL((p2==Param()),true)

	p2 = p0.copyWithInherit("test:","new_prefix");
	TEST_EQUAL(float(p2.getValue("new_prefix:float")), 17.4f);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:float"), "test:float");
	TEST_STRING_EQUAL(p2.getValue("new_prefix:string"), "test,test,test");
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:string"), "test:string");
	TEST_EQUAL(int(p2.getValue("new_prefix:int")), 17);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:int"), "test:int");
	TEST_EQUAL(double(p2.getValue("new_prefix:double")), 18.2);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:double"), "test2:double");
	TEST_EQUAL(p2.getValue("new_prefix:nostring"), DataValue::EMPTY);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:nostring"), String::EMPTY);
	TEST_EQUAL(int(p2.getValue("new_prefix:anotherint")), 99);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:anotherint"), "test3:test3a:anotherint");
	TEST_STRING_EQUAL(p2.getValue("new_prefix:bla"), "blubber");
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:bla"), "test3:test3a:bla");
	TEST_EQUAL(p2.getValue("new_prefix:inherit"), DataValue::EMPTY);
	TEST_STRING_EQUAL(p2.getDescription("new_prefix:inherit"), String::EMPTY);

	Param p3;

	p3.setValue("circle1:inherit","circle2");
	p3.setValue("circle1:iwashere1","incircle1");
	p3.setValue("circle2:inherit","circle3");
	p3.setValue("circle2:iwashere2","incircle2");
	p3.setValue("circle3:inherit","circle4");
	p3.setValue("circle3:iwashere3","incircle3");
	p3.setValue("circle4:inherit","circle1");
	p3.setValue("circle4:iwashere4","incircle4");
	STATUS(p3);

	Param p4;
	// p3.copyWithInherit("circle1:"); // debugging
	TEST_EXCEPTION(Exception::ParseError, p4 = p3.copyWithInherit("circle1:"));
	STATUS(p4);

	p3.remove("circle4:inherit");

	// without new_prefix
	{
		p4 = p3.copyWithInherit("circle1:");
		
		Param p5;
		
		p5.setValue("iwashere1","incircle1");
		p5.setValue("iwashere2","incircle2");
		p5.setValue("iwashere3","incircle3");
		p5.setValue("iwashere4","incircle4");
		STATUS(p5);

		TEST_EQUAL(p4==p5,true);
	}

	// with new_prefix
	{
		p4 = p3.copyWithInherit("circle1:","new_prefix");
		
		Param p5;
		
		p5.setValue("new_prefix:iwashere1","incircle1");
		p5.setValue("new_prefix:iwashere2","incircle2");
		p5.setValue("new_prefix:iwashere3","incircle3");
		p5.setValue("new_prefix:iwashere4","incircle4");
		STATUS(p5);

		TEST_EQUAL(p4==p5,true);
	}
	
}
RESULT

CHECK((void setDefaults(const Param& defaults, String prefix="", bool showMessage=false)))
	Param defaults;
	defaults.setValue("float",1.0f,"float");	
	defaults.setValue("float2",2.0f,"float2");
	defaults.setValue("string","default string1","string");
	defaults.setValue("string2","default string2","string2");
	defaults.setValue("PATH:onlyfordescription",45.2);
	defaults.setDescription("PATH","PATHdesc");
	Param p2;
	p2.setValue("PATH:float",-1.0f,"PATH:float");
	p2.setValue("PATH:string","some string","PATH:string");
	p2.setValue("float",-2.0f,"float");
	p2.setValue("string","other string","string");
	
	TEST_EQUAL(p2.size(),4);
	
	p2.setDefaults(defaults);
	TEST_EQUAL(p2.size(),7);
	TEST_REAL_EQUAL(float(p2.getValue("float")),-2.0);
	TEST_STRING_EQUAL(p2.getDescription("float"),"float");
	TEST_REAL_EQUAL(float(p2.getValue("float2")),2.0);
	TEST_STRING_EQUAL(p2.getDescription("float2"),"float2");
	TEST_EQUAL(string(p2.getValue("string")),"other string");
	TEST_STRING_EQUAL(p2.getDescription("string"),"string");
	TEST_EQUAL(string(p2.getValue("string2")),"default string2");
	TEST_STRING_EQUAL(p2.getDescription("string2"),"string2");
	TEST_STRING_EQUAL(p2.getDescription("PATH"),"PATHdesc");

	p2.setDefaults(defaults,"PATH");
	TEST_EQUAL(p2.size(),10);
	TEST_REAL_EQUAL(float(p2.getValue("PATH:float")),-1.0);
	TEST_STRING_EQUAL(p2.getDescription("PATH:float"),"PATH:float");
	TEST_REAL_EQUAL(float(p2.getValue("PATH:float2")),2.0);
	TEST_STRING_EQUAL(p2.getDescription("PATH:float2"),"float2");
	TEST_EQUAL(string(p2.getValue("PATH:string")),"some string");
	TEST_STRING_EQUAL(p2.getDescription("PATH:string"),"PATH:string");
	TEST_EQUAL(string(p2.getValue("PATH:string2")),"default string2");
	TEST_STRING_EQUAL(p2.getDescription("PATH:string2"),"string2");
	TEST_STRING_EQUAL(p2.getDescription("PATH"),"PATHdesc");
RESULT

char* a1 ="executable";
char* a2 ="-a";
char* a3 ="av";
char* a4 ="-b";
char* a5 ="bv";
char* a6 ="-c";
char* a7 ="cv";
char* a8 ="rv1";
char* a9 ="rv2";
char* command_line[9]; // "executable -a av -b bv -c cv rv1 rv2"
command_line[0] = a1;
command_line[1] = a2;
command_line[2] = a3;
command_line[3] = a4;
command_line[4] = a5;
command_line[5] = a6;
command_line[6] = a7;
command_line[7] = a8;
command_line[8] = a9;

char* command_line2[6]; // "executable -a av -b -c cv"
command_line2[0] = a1;
command_line2[1] = a2;
command_line2[2] = a3;
command_line2[3] = a4;
command_line2[4] = a6;
command_line2[5] = a7;

char* command_line3[6]; // "executable -a -b -c cv rv1"
command_line3[0] = a1;
command_line3[1] = a2;
command_line3[2] = a4;
command_line3[3] = a6;
command_line3[4] = a7;
command_line3[5] = a8;


CHECK((void parseCommandLine(const int argc, char** argv, String prefix = "")))
	Param p2,p3;
	p2.parseCommandLine(9,command_line,"test4");
	p3.setValue("test4:-a","av");
	p3.setValue("test4:-b","bv");
	p3.setValue("test4:-c","cv");
	p3.setValue("test4:misc","rv1 rv2");
	TEST_EQUAL(p2==p3,true)

	Param p20,p30;
	p20.parseCommandLine(6,command_line2);
	p30.setValue("-a","av");
	p30.setValue("-b","");
	p30.setValue("-c","cv");
	TEST_EQUAL(p20==p30,true)
RESULT

CHECK((void parseCommandLine(const int argc, char** argv, const std::map<String, String>& options_with_argument, const std::map<String, String>& options_without_argument, const String& misc="misc", const String& unknown="unknown")))
	map<String,String> with,without;
	with["-a"]="a";
	with["-b"]="b";
	with["-c"]="c";
	
	Param p4,p5;
	p4.parseCommandLine(9,command_line,with,without,"misc_","unknown_");
	p5.setValue("a","av");
	p5.setValue("b","bv");
	p5.setValue("c","cv");
	p5.setValue("misc_","rv1 rv2");
	TEST_EQUAL(p4==p5,true)

	with.clear();
	with["-a"]="a";
	without["-b"]="b";
	
	Param p40,p50;
	p40.parseCommandLine(9,command_line,with,without,"misc__","unknown__");
	p50.setValue("a","av");
	p50.setValue("b","true");
	p50.setValue("misc__","bv cv rv1 rv2");
	p50.setValue("unknown__","-c");
	TEST_EQUAL(p40==p50,true)

	//"executable -a av -b -c cv"	
	Param p400,p500;
	p400.parseCommandLine(6,command_line2,with,without,"misc__","unknown__");
	p500.setValue("a","av");
	p500.setValue("b","true");
	p500.setValue("misc__","cv");
	p500.setValue("unknown__","-c");
	TEST_EQUAL(p400==p500,true)

	//"executable -a -b -c cv rv1"
	Param p4000,p5000;
	p4000.parseCommandLine(6,command_line3,with,without,"misc__","unknown__");
	p5000.setValue("a","");
	p5000.setValue("b","true");
	p5000.setValue("misc__","cv rv1");
	p5000.setValue("unknown__","-c");
	TEST_EQUAL(p4000==p5000,true)
RESULT


CHECK(([EXTRA] ConstIterator begin() const))
	TEST_EQUAL("test2:float", p.begin()->first)
	TEST_EQUAL(p.getValue("test2:float"), p.begin()->second)
RESULT

CHECK(([EXTRA] ConstIterator end() const))
	Param::ConstIterator it = p.end();
	it--;
	TEST_EQUAL("test:string", it->first)
	TEST_EQUAL(p.getValue("test:string"), it->second)
RESULT

CHECK([EXTRA](friend std::ostream& operator << (std::ostream& os, const Param& param)))
	Param p;
	p.setValue("key",17.4);
	stringstream ss;
	ss << p;
	TEST_EQUAL(ss.str(), "\"key\"  ->  \"17.4\"\n")
	p.clear();
	ss.str("");
	p.setValue("key",17.4, "thiskey");
	ss<<p;
	TEST_EQUAL(ss.str(), "\"key\"  ->  \"17.4\" :thiskey\n")
	
RESULT

CHECK((ConstIterator begin() const))
	Param p;
	p.setValue("key",17.4);
	TEST_EQUAL(p.begin()->first, "key")
	TEST_EQUAL(double(p.begin()->second), 17.4)
RESULT

CHECK((ConstIterator end() const))
	Param p;
	TEST_EQUAL(p.end()==p.begin(),true)
	p.setValue("key",17.4);
	TEST_EQUAL((--p.end())==p.begin(),true)
RESULT

CHECK((void checkDefaults(const String &name, const Param &defaults, String prefix="", std::ostream &os=std::cout) const))
	ostringstream os;
	Param p,d;
	p.setValue("string",String("bla"),"string");
	p.setValue("int",5,"int");
	p.setValue("double",47.11,"double");
		
	p.checkDefaults("Test",d,"",os);
	TEST_EQUAL(os.str()=="Warning: Test received the unknown parameter 'double'!\nWarning: Test received the unknown parameter 'int'!\nWarning: Test received the unknown parameter 'string'!\n",true);
	
	d.setValue("int",5,"int");
	d.setValue("double",47.11,"double");
	os.str("");
	p.checkDefaults("Test",d,"",os);
	TEST_EQUAL(os.str()=="Warning: Test received the unknown parameter 'string'!\n",true);
	
	p.clear();
	p.setValue("pref:string",String("bla"),"pref:string");
	p.setValue("pref:int",5,"pref:int");
	p.setValue("pref:double",47.11,"pref:double");
	os.str("");
	p.checkDefaults("Test",d,"pref",os);
	TEST_EQUAL(os.str()=="Warning: Test received the unknown parameter 'string' in 'pref:'!\n",true);

	os.str("");
	p.checkDefaults("Test",d,"pref:",os);
	TEST_EQUAL(os.str()=="Warning: Test received the unknown parameter 'string' in 'pref:'!\n",true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



