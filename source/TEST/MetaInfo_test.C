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
// $Id: MetaInfo_test.C,v 1.6 2006/05/30 15:46:43 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfo.h>

///////////////////////////

START_TEST(Example, "$Id: MetaInfo_test.C,v 1.6 2006/05/30 15:46:43 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;

MetaInfo* test;

CHECK(MetaInfo::MetaInfo())
	test = new MetaInfo;
	TEST_NOT_EQUAL(test, 0)
RESULT

CHECK(~MetaInfo())
	delete test;
RESULT

MetaInfo mi;

CHECK(MetaInfoRegistry& registry())
	MetaInfo mi2;
	mi2.registry().registerName("testname","testdesc","testunit");
	TEST_EQUAL (mi2.registry().getIndex("testname"),1024);
	TEST_EQUAL (mi.registry().getIndex("testname"),1024);
RESULT

CHECK (setValue(UnsignedInt index, const std::string& value) / DataValue getValue(UnsignedInt index) const)
	string tmp;
	mi.setValue(1024,string("testtesttest"));
	tmp = string(mi.getValue(1024));
	TEST_EQUAL("testtesttest",tmp)
RESULT

CHECK (setValue(const std::string& name, const std::string& value) / DataValue getValue(const std::string& name) const)
	string tmp;
	mi.setValue("testname",string("testtesttest2"));
	tmp = string(mi.getValue("testname"));
	TEST_EQUAL("testtesttest2",tmp)
RESULT

CHECK (void setValue(const std::string& name, SignedInt))
	SignedInt tmp;
	mi.setValue("cluster_id",4711);
	tmp = SignedInt(mi.getValue("cluster_id"));
	TEST_EQUAL(tmp,4711)
RESULT

CHECK (void setValue(const std::string& name, double))
	double tmp;
	mi.setValue("cluster_id",4711.1234);
	tmp = double(mi.getValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4711.1234)
RESULT

CHECK (void setValue(const std::string& name, SignedInt))
	SignedInt tmp;
	mi.setValue(2,4712);
	tmp = SignedInt(mi.getValue("cluster_id"));
	TEST_EQUAL(tmp,4712)
RESULT

CHECK (void setValue(const std::string& name, double))
	double tmp;
	mi.setValue(2,4712.1234);
	tmp = double(mi.getValue("cluster_id"));
	TEST_REAL_EQUAL(tmp,4712.1234)
RESULT

CHECK (bool empty())
	MetaInfo tmp;
	TEST_EQUAL(tmp.empty(),true)
	tmp.setValue(1024,string("testtesttest"));
	TEST_EQUAL(tmp.empty(),false)
RESULT

CHECK(MetaInfo(const MetaInfo& rhs))
	MetaInfo mi3(mi);
	TEST_REAL_EQUAL(double(mi3.getValue("cluster_id")),4712.1234)
	TEST_EQUAL("testtesttest2",string(mi3.getValue("testname")))
RESULT

CHECK(MetaInfo(const MetaInfo& rhs))
	MetaInfo mi4;
	mi4 = mi;
	TEST_REAL_EQUAL(double(mi4.getValue("cluster_id")),4712.1234)
	TEST_EQUAL("testtesttest2",string(mi4.getValue("testname")))
RESULT

CHECK (setValue(const std::string& name, const DataValue& value) / DataValue getValue(const std::string& name) const)
	DataValue tmp("testtesttest3");
	mi.setValue("testname",tmp);
	tmp = string(mi.getValue("testname"));
	TEST_EQUAL("testtesttest3",tmp)
RESULT

CHECK (void getKeys(std::vector<std::string>& keys) const)
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

CHECK (bool exists(const std::string& name))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists("cluster_id"),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists("cluster_id"),true)
RESULT

CHECK (bool exists(UnsignedInt index))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists(2),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists(2),true)
RESULT

CHECK (void clear())
	MetaInfo i;
	TEST_EQUAL(i.empty(),true)
	i.setValue("label",String("test"));
	TEST_EQUAL(i.empty(),false)
	i.clear();
	TEST_EQUAL(i.empty(),true)
RESULT

CHECK (bool operator== (const MetaInfo& rhs) const)
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

CHECK (bool operator!= (const MetaInfo& rhs) const)
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

