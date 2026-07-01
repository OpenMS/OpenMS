// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfo.h>

///////////////////////////

START_TEST(Example, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;

MetaInfo* test = nullptr;
MetaInfo* nullPointer = nullptr;

START_SECTION((MetaInfo()))
	test = new MetaInfo;
  TEST_NOT_EQUAL(test, nullPointer)
END_SECTION

START_SECTION((~MetaInfo()))
	delete test;
END_SECTION

MetaInfo mi;

START_SECTION((static MetaInfoRegistry& registry()))
	MetaInfo mi2;
	mi2.registry().registerName("testname", "testdesc", "testunit");
	TEST_EQUAL(mi2.registry().getIndex("testname"), 1024);
	TEST_EQUAL(mi.registry().getIndex("testname"), 1024);
END_SECTION

START_SECTION((void setValue(const std::string& name, const DataValue& value)))
	NOT_TESTABLE //tested in the get method
END_SECTION

START_SECTION((void setValue(UInt index, const DataValue& value)))
	NOT_TESTABLE //tested in the get method
END_SECTION

START_SECTION((const DataValue& getValue(UInt index, const DataValue& default_value = DataValue::EMPTY) const))
{
	string tmp;
	mi.setValue(1024,std::string("testtesttest"));
	tmp =StringUtils::toStr(mi.getValue(1024));
  TEST_EQUAL(tmp, "testtesttest");
	TEST_EQUAL(mi.getValue(1025) == DataValue::EMPTY, true);
	TEST_EQUAL(mi.getValue(1025, 10) == DataValue(10), true);
}
END_SECTION

START_SECTION((const DataValue& getValue(const std::string& name, const DataValue& default_value = DataValue::EMPTY) const))
{
	string tmp;
	mi.setValue("testname",std::string("testtesttest2"));
	tmp =StringUtils::toStr(mi.getValue("testname"));
	TEST_EQUAL(tmp, "testtesttest2");
	TEST_EQUAL(mi.getValue("notdefined") == DataValue::EMPTY, true);
	TEST_EQUAL(mi.getValue("notdefined", 10) == DataValue(10), true);
}
END_SECTION

mi.setValue("cluster_id",4711.12f);
mi.setValue(2,4712.12f);

START_SECTION((bool empty() const))
	MetaInfo tmp;
	TEST_EQUAL(tmp.empty(),true)
	tmp.setValue(1024,std::string("testtesttest"));
	TEST_EQUAL(tmp.empty(),false)
END_SECTION

START_SECTION((MetaInfo(const MetaInfo& rhs)))
	MetaInfo mi3(mi);
	TEST_REAL_SIMILAR(double(mi3.getValue("cluster_id")),double(mi.getValue("cluster_id")))
	TEST_STRING_EQUAL(mi3.getValue("testname"),"testtesttest2")
END_SECTION

START_SECTION((MetaInfo& operator = (const MetaInfo& rhs)))
	MetaInfo mi3;
	mi3 = mi;
	TEST_REAL_SIMILAR(double(mi3.getValue("cluster_id")),double(mi.getValue("cluster_id")))
	TEST_STRING_EQUAL(mi3.getValue("testname"),"testtesttest2")
END_SECTION

START_SECTION((void getKeys(std::vector<std::string>& keys) const))
	vector<std::string> tmp,tmp2;
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
END_SECTION

START_SECTION((void getKeys(std::vector< UInt > &keys) const))
	MetaInfo mi;
	mi.setValue("label",std::string("tag"));
	mi.setValue("icon",std::string("kreis"));
	vector<UInt> vec;
	mi.getKeys(vec);
	TEST_EQUAL(vec.size(),2)
	TEST_EQUAL(vec[0],3)
	TEST_EQUAL(vec[1],4)

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
END_SECTION

START_SECTION((bool exists(const std::string& name) const))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists("cluster_id"),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists("cluster_id"),true)
END_SECTION

START_SECTION((bool exists(UInt index) const))
	MetaInfo mi4;
	TEST_EQUAL(mi4.exists(2),false)
	mi4.setValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.exists(2),true)
END_SECTION

START_SECTION((void clear()))
	MetaInfo i;
	TEST_EQUAL(i.empty(),true)
	i.setValue("label",std::string("test"));
	TEST_EQUAL(i.empty(),false)
	i.clear();
	TEST_EQUAL(i.empty(),true)
END_SECTION

START_SECTION((bool operator== (const MetaInfo& rhs) const))
	MetaInfo i,i2;
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
	i.setValue("label",std::string("test"));
	TEST_EQUAL(i==i2,false)
	TEST_EQUAL(i2==i,false)
	i2.setValue("label",std::string("test"));
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
END_SECTION

START_SECTION((bool operator!= (const MetaInfo& rhs) const))
	MetaInfo i,i2;
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
	i.setValue("label",std::string("test"));
	TEST_EQUAL(i!=i2,true)
	TEST_EQUAL(i2!=i,true)
	i2.setValue("label",std::string("test"));
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
END_SECTION

START_SECTION((MetaInfo & operator+=(const MetaInfo& rhs)))
  MetaInfo m_new, m_base;

  m_base.setValue("label",std::string("old"));
  m_base.setValue("exists_no_overwrite", 5.2);
  
  m_new.setValue("label",std::string("new")); // will be overwritten
  m_new.setValue("icon", 4.3);            // will be added
  
  m_base += m_new;

  TEST_EQUAL(m_base.getValue("label"),std::string("new"));
  TEST_EQUAL(m_base.getValue("icon"), 4.3);
  TEST_EQUAL(m_base.getValue("exists_no_overwrite"), 5.2);
END_SECTION

START_SECTION((void removeValue(UInt index)))
	MetaInfo i,i2;

	i.setValue(1,std::string("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeValue(1);
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeValue(1234);
END_SECTION

START_SECTION((void removeValue(const std::string& name)))
	MetaInfo i,i2;

	i.setValue("label",std::string("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeValue("label");
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeValue("icon");
END_SECTION

START_SECTION((Size size() const))
	MetaInfo mi_size;
	TEST_EQUAL(mi_size.size(), 0)
	mi_size.setValue("key1", 1);
	TEST_EQUAL(mi_size.size(), 1)
	mi_size.setValue("key2", 2);
	TEST_EQUAL(mi_size.size(), 2)
	mi_size.removeValue("key1");
	TEST_EQUAL(mi_size.size(), 1)
	mi_size.clear();
	TEST_EQUAL(mi_size.size(), 0)
END_SECTION

START_SECTION((const_iterator begin() const))
	MetaInfo mi_iter;
	TEST_EQUAL(mi_iter.begin() == mi_iter.end(), true) // empty
	mi_iter.setValue("test_key", 42);
	TEST_EQUAL(mi_iter.begin() != mi_iter.end(), true)
	TEST_EQUAL(mi_iter.begin()->second, DataValue(42))
END_SECTION

START_SECTION((const_iterator end() const))
	NOT_TESTABLE // tested with begin()
END_SECTION

START_SECTION((const_iterator cbegin() const))
	MetaInfo mi_citer;
	mi_citer.setValue("key",std::string("value"));
	auto it = mi_citer.cbegin();
	TEST_EQUAL(it->second, DataValue(std::string("value")))
END_SECTION

START_SECTION((const_iterator cend() const))
	NOT_TESTABLE // tested with cbegin()
END_SECTION

START_SECTION((iterator begin()))
	MetaInfo mi_mut;
	mi_mut.setValue("mutable_key", 100);
	auto it = mi_mut.begin();
	it->second = DataValue(200);
	TEST_EQUAL(mi_mut.getValue("mutable_key"), DataValue(200))
END_SECTION

START_SECTION((iterator end()))
	NOT_TESTABLE // tested with mutable begin()
END_SECTION

START_SECTION([EXTRA] Range-based for loop iteration)
	MetaInfo mi_range;
	mi_range.setValue("a", 1);
	mi_range.setValue("b", 2);
	mi_range.setValue("c", 3);

	int count = 0;
	for (const auto& kv : mi_range)
	{
		++count;
		TEST_EQUAL(kv.second.valueType() == DataValue::INT_VALUE, true)
	}
	TEST_EQUAL(count, 3)
END_SECTION

START_SECTION([EXTRA] operator+= with empty MetaInfo)
	MetaInfo m1, m2;
	m1.setValue("key", 1);

	// Adding empty to non-empty
	m1 += m2;
	TEST_EQUAL(m1.size(), 1)
	TEST_EQUAL(m1.getValue("key"), DataValue(1))

	// Adding non-empty to empty
	MetaInfo m3;
	m3 += m1;
	TEST_EQUAL(m3.size(), 1)
	TEST_EQUAL(m3.getValue("key"), DataValue(1))
END_SECTION

START_SECTION([EXTRA] operator+= with disjoint keys)
	MetaInfo m1, m2;
	m1.setValue("a", 1);
	m1.setValue("c", 3);
	m2.setValue("b", 2);
	m2.setValue("d", 4);

	m1 += m2;
	TEST_EQUAL(m1.size(), 4)
	TEST_EQUAL(m1.getValue("a"), DataValue(1))
	TEST_EQUAL(m1.getValue("b"), DataValue(2))
	TEST_EQUAL(m1.getValue("c"), DataValue(3))
	TEST_EQUAL(m1.getValue("d"), DataValue(4))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

