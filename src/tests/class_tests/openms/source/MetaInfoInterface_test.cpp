// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfoInterface.h>

///////////////////////////

START_TEST(Example, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace std;
using namespace OpenMS;

MetaInfoInterface* test = nullptr;
MetaInfoInterface* nullPointer = nullptr;
START_SECTION((MetaInfoInterface()))
	test = new MetaInfoInterface;
  TEST_NOT_EQUAL(test, nullPointer)
END_SECTION


START_SECTION((~MetaInfoInterface()))
	delete test;
END_SECTION

MetaInfoInterface mi;

START_SECTION((static MetaInfoRegistry& metaRegistry()))
	MetaInfoInterface mi2;
	mi2.metaRegistry().registerName("testname","testdesc","testunit");
	TEST_EQUAL(mi2.metaRegistry().getIndex("testname"),1024);
	TEST_EQUAL(mi.metaRegistry().getIndex("testname"),1024);
END_SECTION

START_SECTION((void setMetaValue(const String& name, const DataValue& value)))
	NOT_TESTABLE //tested in the get method
END_SECTION

START_SECTION((void setMetaValue(UInt index, const DataValue& value)))
	NOT_TESTABLE //tested in the get method
END_SECTION

START_SECTION((const DataValue& getMetaValue(UInt index) const))
	mi.setMetaValue(1024,String("testtesttest"));
	TEST_STRING_EQUAL(mi.getMetaValue(1024),"testtesttest");
END_SECTION

START_SECTION((const DataValue& getMetaValue(const String& name) const))
	mi.setMetaValue("testname",String("testtesttest2"));
	TEST_STRING_EQUAL(mi.getMetaValue("testname"),"testtesttest2");
END_SECTION

	mi.setMetaValue("cluster_id",4711.12f);
	mi.setMetaValue(2,4712.12f);

START_SECTION((bool isMetaEmpty() const))
	MetaInfoInterface tmp;
	TEST_EQUAL(tmp.isMetaEmpty(),true)
	tmp.setMetaValue(1024,String("testtesttest"));
	TEST_EQUAL(tmp.isMetaEmpty(),false)
END_SECTION

TOLERANCE_ABSOLUTE(0.001)

START_SECTION((MetaInfoInterface(const MetaInfoInterface& rhs)))
	//test if copy worked
	MetaInfoInterface mi3(mi);
	TEST_REAL_SIMILAR(double(mi.getMetaValue("cluster_id")),double(mi3.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_SIMILAR(double(mi.getMetaValue("cluster_id")),4712.12)
	TEST_REAL_SIMILAR(double(mi3.getMetaValue("cluster_id")),11.9)
END_SECTION

START_SECTION((MetaInfoInterface& operator = (const MetaInfoInterface& rhs)))
	//test if copy worked
	MetaInfoInterface mi3,mi4;
	mi3 = mi;
	TEST_REAL_SIMILAR(double(mi3.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_SIMILAR(double(mi.getMetaValue("cluster_id")),4712.12)
	TEST_REAL_SIMILAR(double(mi3.getMetaValue("cluster_id")),11.9)
	//test what happens when left side is not empty
	mi3 = mi;
	TEST_REAL_SIMILAR(double(mi3.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	//test if a deep copy was done
	mi3.setMetaValue("cluster_id",11.9);
	TEST_REAL_SIMILAR(double(mi.getMetaValue("cluster_id")),double(mi.getMetaValue("cluster_id")))
	TEST_REAL_SIMILAR(double(mi3.getMetaValue("cluster_id")),11.9)
	//test what happens when source is empty
	mi3 = mi4;
	TEST_EQUAL(mi3.isMetaEmpty(),true)
END_SECTION

START_SECTION((void getKeys(std::vector<String>& keys) const))
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
END_SECTION

START_SECTION((void getKeys(std::vector< UInt > &keys) const))
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
END_SECTION

START_SECTION((bool metaValueExists(const String& name) const))
	MetaInfoInterface mi4;
	TEST_EQUAL(mi4.metaValueExists("cluster_id"),false)
	mi4.setMetaValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.metaValueExists("cluster_id"),true)
END_SECTION

START_SECTION((bool metaValueExists(UInt index) const))
	MetaInfoInterface mi4;
	TEST_EQUAL(mi4.metaValueExists(2),false)
	mi4.setMetaValue("cluster_id",4712.1234);
	TEST_EQUAL(mi4.metaValueExists(2),true)
END_SECTION

START_SECTION(([EXTRA] void getKeys(std::vector<String>& keys) const))
	std::vector<String> keys;
	mi.getKeys(keys);
	TEST_EQUAL(keys.size(),2)
	TEST_EQUAL(keys[0],"cluster_id")
	TEST_EQUAL(keys[1],"testname")
END_SECTION

START_SECTION((void clearMetaInfo()))
	MetaInfoInterface i;
	TEST_EQUAL(i.isMetaEmpty(),true)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i.isMetaEmpty(),false)
	i.clearMetaInfo();
	TEST_EQUAL(i.isMetaEmpty(),true)
END_SECTION

START_SECTION((bool operator== (const MetaInfoInterface& rhs) const))
	MetaInfoInterface i,i2;
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i==i2,false)
	TEST_EQUAL(i2==i,false)
	i2.setMetaValue("label",String("test"));
	TEST_EQUAL(i==i2,true)
	TEST_EQUAL(i2==i,true)
END_SECTION

START_SECTION((bool operator!= (const MetaInfoInterface& rhs) const))
	MetaInfoInterface i,i2;
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
	i.setMetaValue("label",String("test"));
	TEST_EQUAL(i!=i2,true)
	TEST_EQUAL(i2!=i,true)
	i2.setMetaValue("label",String("test"));
	TEST_EQUAL(i!=i2,false)
	TEST_EQUAL(i2!=i,false)
END_SECTION

START_SECTION((void removeMetaValue(UInt index)))
	MetaInfoInterface i,i2;

	i.setMetaValue(1,String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeMetaValue(1);
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeMetaValue(1234);
END_SECTION

START_SECTION((void removeMetaValue(const String& name)))
	MetaInfoInterface i,i2;

	i.setMetaValue("label",String("bla"));
	TEST_EQUAL(i==i2,false)
	i.removeMetaValue("label");
	TEST_EQUAL(i==i2,true)

	//try if removing a non-existing value works as well
	i.removeMetaValue("icon");
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

