// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfoInterface.h>

///////////////////////////

START_TEST(Example, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;

MetaInfoInterface* test;
CHECK((MetaInfoInterface()))
	test = new MetaInfoInterface;
	TEST_NOT_EQUAL(test, 0)
RESULT


CHECK((~MetaInfoInterface()))
	delete test;
RESULT

MetaInfoInterface mi;

CHECK((MetaInfoRegistry& metaRegistry() const))
	MetaInfoInterface mi2;
	mi2.metaRegistry().registerName("testname","testdesc","testunit");
	TEST_EQUAL(mi2.metaRegistry().getIndex("testname"),1024);
	TEST_EQUAL(mi.metaRegistry().getIndex("testname"),1024);
RESULT

CHECK((void setMetaValue(UInt index, const String& value)))
	NOT_TESTABLE //tested in the get method
RESULT

CHECK((void setMetaValue(const String& name, const String& value)))
	NOT_TESTABLE //tested in the get method
RESULT

CHECK((const DataValue& getMetaValue(UInt index) const))
	string tmp;
	mi.setMetaValue(1024,String("testtesttest"));
	tmp = String(mi.getMetaValue(1024));
	TEST_EQUAL("testtesttest",tmp)
RESULT

CHECK((const DataValue& getMetaValue(const String& name) const))
	string tmp;
	mi.setMetaValue("testname",String("testtesttest2"));
	tmp = String(mi.getMetaValue("testname"));
	TEST_EQUAL("testtesttest2",tmp)
RESULT

CHECK((void setMetaValue(const String& name, Int value)))
	Int tmp;
	mi.setMetaValue("cluster_id",-4711);
	tmp = Int(mi.getMetaValue("cluster_id"));
	TEST_EQUAL(tmp,-4711)
RESULT

CHECK((void setMetaValue(const String& name, DoubleReal value)))
	double tmp;
	mi.setMetaValue("cluster_id",4711.1234);
	tmp = double(mi.getMetaValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4711.1234)
RESULT

CHECK((void setMetaValue(UInt index, Int value)))
	Int tmp;
	mi.setMetaValue(2,-4712);
	tmp = Int(mi.getMetaValue("cluster_id"));
	TEST_EQUAL(tmp,-4712)
RESULT

CHECK((void setMetaValue(UInt index, DoubleReal value)))
	double tmp;
	mi.setMetaValue(2,4712.1234);
	tmp = double(mi.getMetaValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4712.1234)
RESULT

CHECK((void setMetaValue(const String& name, UInt value)))
	Int tmp;
	mi.setMetaValue("cluster_id",4711u);
	tmp = Int(mi.getMetaValue("cluster_id"));
	TEST_EQUAL(tmp,4711u)
RESULT

CHECK((void setMetaValue(const String& name, Real value)))
	double tmp;
	mi.setMetaValue("cluster_id",4711.12f);
	tmp = double(mi.getMetaValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4711.12f)
RESULT

CHECK((void setMetaValue(UInt index, UInt value)))
	Int tmp;
	mi.setMetaValue(2,4712u);
	tmp = Int(mi.getMetaValue("cluster_id"));
	TEST_EQUAL(tmp,4712u)
RESULT

CHECK((void setMetaValue(UInt index, Real value)))
	double tmp;
	mi.setMetaValue(2,4712.12f);
	tmp = double(mi.getMetaValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4712.12f)
RESULT

CHECK((bool isMetaEmpty() const))
	MetaInfoInterface tmp;
	TEST_EQUAL(tmp.isMetaEmpty(),true)
	tmp.setMetaValue(1024,String("testtesttest"));
	TEST_EQUAL(tmp.isMetaEmpty(),false)
RESULT

PRECISION(0.001)

CHECK((MetaInfoInterface(const MetaInfoInterface& rhs)))
	//test if copy worked
	MetaInfoInterface mi3(mi);
	TEST_REAL_EQUAL(double(mi.getMetaValue("cluster_id")),double(mi3.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_EQUAL(double(mi.getMetaValue("cluster_id")),4712.12)
	TEST_REAL_EQUAL(double(mi3.getMetaValue("cluster_id")),11.9)
RESULT

CHECK((MetaInfoInterface& operator = (const MetaInfoInterface& rhs)))
	//test if copy worked
	MetaInfoInterface mi3,mi4;
	mi3 = mi;
	TEST_REAL_EQUAL(double(mi3.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_EQUAL(double(mi.getMetaValue("cluster_id")),4712.12)
	TEST_REAL_EQUAL(double(mi3.getMetaValue("cluster_id")),11.9)
	//test what happens when left side is not empty
	mi3 = mi;
	TEST_REAL_EQUAL(double(mi3.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_EQUAL(double(mi.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	TEST_REAL_EQUAL(double(mi3.getMetaValue("cluster_id")),11.9)	
	//test what happens when source is empty
	mi3 = mi4;
	TEST_REAL_EQUAL(mi3.isMetaEmpty(),true)	
RESULT

CHECK((void setMetaValue(const String& name, const DataValue& value)))
	DataValue tmp("testtesttest3");
	mi.setMetaValue("testname",tmp);
	tmp = String(mi.getMetaValue("testname"));
	TEST_EQUAL("testtesttest3",tmp)
RESULT

CHECK((void setMetaValue(UInt index, const DataValue& value)))
	DataValue tmp("testtesttest3");
	mi.setMetaValue(2,tmp);
	tmp = String(mi.getMetaValue(2));
	TEST_EQUAL("testtesttest3",tmp)
RESULT

CHECK((void getKeys(std::vector<String>& keys) const))
	vector<String> tmp,tmp2;
	tmp.push_back("cluster_id");
	tmp.push_back("testname");
	mi.getKeys(tmp2);
	TEST_EQUAL(tmp2.size(),tmp.size())
	TEST_EQUAL(tmp2[0],tmp[0])
	TEST_EQUAL(tmp2[1],tmp[1])
	
	MetaInfoInterface mi2(mi);
	mi2.getKeys(tmp2);	
	TEST_EQUAL(tmp2.size(),tmp.size())
	TEST_EQUAL(tmp2[0],tmp[0])
	TEST_EQUAL(tmp2[1],tmp[1])

	mi2.metaRegistry().registerName("a","test");
	mi2.metaRegistry().registerName("d","test");
	mi2.metaRegistry().registerName("x","test");
	mi2.setMetaValue("a",1);
	mi2.setMetaValue("d",1);
	mi2.setMetaValue("x",1);
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

CHECK((void getKeys(std::vector< UInt > &keys) const))
	MetaInfoInterface mi;
	mi.setMetaValue("label",String("tag"));
	mi.setMetaValue("icon",String("kreis"));
	vector<UInt> vec;
	mi.getKeys(vec);
	TEST_EQUAL(vec.size(),2)
	TEST_EQUAL(vec[0],3)
	TEST_EQUAL(vec[1],4)
	
	mi.metaRegistry().registerName("a","test");
	mi.metaRegistry().registerName("d","test");
	mi.metaRegistry().registerName("x","test");
	mi.setMetaValue("a",1);
	mi.setMetaValue("d",1);
	mi.setMetaValue("x",1);
	mi.getKeys(vec);
	
	TEST_EQUAL(vec.size(),5)
	TEST_EQUAL(vec[0],3)
	TEST_EQUAL(vec[1],4)
	TEST_EQUAL(vec[2],1025)
	TEST_EQUAL(vec[3],1026)
	TEST_EQUAL(vec[4],1027)
RESULT

CHECK((bool metaValueExists(const String& name) const))
	MetaInfoInterface mi4;
	TEST_EQUAL(mi4.metaValueExists("cluster_id"),false)
	mi4.setMetaValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.metaValueExists("cluster_id"),true)
RESULT

CHECK((bool metaValueExists(UInt index) const))
	MetaInfoInterface mi4;
	TEST_EQUAL(mi4.metaValueExists(2),false)
	mi4.setMetaValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.metaValueExists(2),true)
RESULT

CHECK(([EXTRA] void getKeys(std::vector<String>& keys) const))
	std::vector<String> keys;
	mi.getKeys(keys);	
	TEST_EQUAL(keys.size(),2)
	TEST_EQUAL(keys[0],"cluster_id")
	TEST_EQUAL(keys[1],"testname")
RESULT

CHECK((void clearMetaInfo()))
	MetaInfoInterface i;
	TEST_EQUAL(i.isMetaEmpty(),true)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i.isMetaEmpty(),false)
	i.clearMetaInfo();
	TEST_EQUAL(i.isMetaEmpty(),true)
RESULT

CHECK((bool operator== (const MetaInfoInterface& rhs) const))
	MetaInfoInterface i,i2;
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i==i2,false)
	TEST_EQUAL(i2==i,false)
	i2.setMetaValue("label",String("test"));
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
RESULT

CHECK((bool operator!= (const MetaInfoInterface& rhs) const))
	MetaInfoInterface i,i2;
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i!=i2,true)
	TEST_EQUAL(i2!=i,true)
	i2.setMetaValue("label",String("test"));
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
RESULT

CHECK((void removeMetaValue(UInt index)))
	MetaInfoInterface i,i2;
	
	i.setMetaValue(1,String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeMetaValue(1);
	TEST_EQUAL(i==i2,true)
	
	//try if removing a non-existing value works as well
	i.removeMetaValue(1234);
RESULT

CHECK((void removeMetaValue(const String& name)))
	MetaInfoInterface i,i2;
	
	i.setMetaValue("label",String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeMetaValue("label");
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeMetaValue("icon");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

