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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfo.h>

///////////////////////////

START_TEST(Example, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;

MetaInfo* test;

CHECK((MetaInfo()))
	test = new MetaInfo;
	TEST_NOT_EQUAL(test, 0)
RESULT

CHECK((~MetaInfo()))
	delete test;
RESULT

MetaInfo mi;

CHECK((static MetaInfoRegistry& registry()))
	MetaInfo mi2;
	mi2.registry().registerName("testname","testdesc","testunit");
	TEST_EQUAL (mi2.registry().getIndex("testname"),1024);
	TEST_EQUAL (mi.registry().getIndex("testname"),1024);
RESULT

CHECK((void setValue(UnsignedInt index, const std::string& value)))
	string tmp;
	mi.setValue(1024,string("testtesttest"));
RESULT

CHECK((const DataValue& getValue(UnsignedInt index) const))
	string tmp;
	mi.setValue(1024,string("testtesttest"));
	tmp = string(mi.getValue(1024));
	TEST_EQUAL("testtesttest",tmp)
RESULT

CHECK((void setValue(const std::string& name, const std::string& value)))
	string tmp;
	mi.setValue("testname",string("testtesttest2"));
RESULT

CHECK((const DataValue& getValue(const std::string& name) const))
	string tmp;
	mi.setValue("testname",string("testtesttest2"));
	tmp = string(mi.getValue("testname"));
	TEST_EQUAL("testtesttest2",tmp)
RESULT

CHECK((void setValue(const std::string& name, SignedInt value)))
	SignedInt tmp;
	mi.setValue("cluster_id",4711);
	tmp = SignedInt(mi.getValue("cluster_id"));
	TEST_EQUAL(tmp,4711)
RESULT

CHECK((void setValue(const std::string& name, double value)))
	double tmp;
	mi.setValue("cluster_id",4711.1234);
	tmp = double(mi.getValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4711.1234)
RESULT

CHECK((void setValue(UnsignedInt index, SignedInt value)))
	SignedInt tmp;
	mi.setValue(2,4712);
	tmp = SignedInt(mi.getValue("cluster_id"));
	TEST_EQUAL(tmp,4712)
RESULT

CHECK((void setValue(UnsignedInt index, double value)))
	double tmp;
	mi.setValue(2,4712.1234);
	tmp = double(mi.getValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4712.1234)
RESULT

CHECK((bool empty() const))
	MetaInfo tmp;
	TEST_EQUAL(tmp.empty(),true)
	tmp.setValue(1024,string("testtesttest"));
	TEST_EQUAL(tmp.empty(),false)
RESULT

CHECK((MetaInfo(const MetaInfo& rhs)))
	MetaInfo mi3(mi);
	TEST_REAL_EQUAL(double(mi3.getValue("cluster_id")),4712.1234)
	TEST_EQUAL("testtesttest2",string(mi3.getValue("testname")))
RESULT

CHECK((MetaInfo& operator = (const MetaInfo& rhs)))
	MetaInfo mi3;
	mi3 = mi;
	TEST_REAL_EQUAL(double(mi3.getValue("cluster_id")),4712.1234)
	TEST_EQUAL("testtesttest2",string(mi3.getValue("testname")))
RESULT

CHECK((void setValue(const std::string& name, const DataValue& value)))
	DataValue tmp("testtesttest3");
	mi.setValue("testname",tmp);
	tmp = string(mi.getValue("testname"));
	TEST_EQUAL("testtesttest3",tmp)
RESULT

CHECK((void setValue(UnsignedInt index, const DataValue& value)))
	DataValue tmp("testtesttest3");
	mi.setValue(2,tmp);
	tmp = string(mi.getValue(2));
	TEST_EQUAL("testtesttest3",tmp)
RESULT

CHECK((void getKeys(std::vector<std::string>& keys) const))
	vector<string> tmp,tmp2;
	tmp.push_back("cluster_id");
	tmp.push_back("testname");
	mi.getKeys(tmp2);
	TEST_EQUAL(tmp2.size(),tmp.size())
	TEST_EQUAL(tmp2[0],tmp[0])
	TEST_EQUAL(tmp2[1],tmp[1])
	
	MetaInfo mi2(mi);
	mi2.getKeys(tmp2);	
	TEST_EQUAL(tmp2.size(),tmp.size())
	TEST_EQUAL(tmp2[0],tmp[0])
	TEST_EQUAL(tmp2[1],tmp[1])

	mi2.registry().registerName("a","test");
	mi2.registry().registerName("d","test");
	mi2.registry().registerName("x","test");
	mi2.setValue("a",1);
	mi2.setValue("d",1);
	mi2.setValue("x",1);
	mi2.getKeys(tmp2);
	tmp.clear();
	tmp.push_back("cluster_id");
	tmp.push_back("testname");
	tmp.push_back("a");
	tmp.push_back("d");
	tmp.push_back("x");	
	
	TEST_EQUAL(tmp2.size(),tmp.size())
	TEST_EQUAL(tmp2[0],tmp[0])
	TEST_EQUAL(tmp2[1],tmp[1])
	TEST_EQUAL(tmp2[2],tmp[2])
	TEST_EQUAL(tmp2[3],tmp[3])
	TEST_EQUAL(tmp2[4],tmp[4])
RESULT

CHECK(void getKeys(std::vector< UnsignedInt > &keys) const)
	MetaInfo mi;
	mi.setValue("label",String("tag"));
	mi.setValue("icon",String("kreis"));
	vector<UnsignedInt> vec;
	mi.getKeys(vec);
	TEST_EQUAL(vec.size(),2)
	TEST_EQUAL(vec[0],3)
	TEST_EQUAL(vec[1],4)
	
	mi.registry().registerName("a","test");
	mi.registry().registerName("d","test");
	mi.registry().registerName("x","test");
	mi.setValue("a",1);
	mi.setValue("d",1);
	mi.setValue("x",1);
	mi.getKeys(vec);
	
	TEST_EQUAL(vec.size(),5)
	TEST_EQUAL(vec[0],3)
	TEST_EQUAL(vec[1],4)
	TEST_EQUAL(vec[2],1025)
	TEST_EQUAL(vec[3],1026)
	TEST_EQUAL(vec[4],1027)
RESULT

CHECK((bool exists(const std::string& name) const))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists("cluster_id"),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists("cluster_id"),true)
RESULT

CHECK((bool exists(UnsignedInt index) const))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists(2),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists(2),true)
RESULT

CHECK((void clear()))
	MetaInfo i;
	TEST_EQUAL(i.empty(),true)
	i.setValue("label",String("test"));
	TEST_EQUAL(i.empty(),false)
	i.clear();
	TEST_EQUAL(i.empty(),true)
RESULT

CHECK((bool operator== (const MetaInfo& rhs) const))
	MetaInfo i,i2;
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
	i.setValue("label",String("test"));
	TEST_EQUAL(i==i2,false)
	TEST_EQUAL(i2==i,false)
	i2.setValue("label",String("test"));
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
RESULT

CHECK((bool operator!= (const MetaInfo& rhs) const))
	MetaInfo i,i2;
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
	i.setValue("label",String("test"));
	TEST_EQUAL(i!=i2,true)
	TEST_EQUAL(i2!=i,true)
	i2.setValue("label",String("test"));
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
RESULT

CHECK((void removeValue(UnsignedInt index)))
	MetaInfo i,i2;
	
	i.setValue(1,String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeValue(1);
	TEST_EQUAL(i==i2,true)
	
	//try if removing a non-existing value works as well
	i.removeValue(1234);
RESULT

CHECK((void removeValue(const std::string& name)))
	MetaInfo i,i2;
	
	i.setValue("label",String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeValue("label");
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeValue("icon");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

