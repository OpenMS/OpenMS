// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h> // for "ParameterInformation"

///////////////////////////

using namespace OpenMS;
using namespace std;

#ifdef _MSC_VER  // disable optimization in VS only for this test (as its size triggers 'heap-overflow' during compile otherwise)
#pragma warning (disable: 4748) // disable warning that occurs when switching optimzation off (as /GS is still enabled)
#pragma optimize( "", off )
#endif

START_TEST(Param, "$Id$")

//////////////////// Param::ParamEntry /////////////////////////////
////////////////////////////////////////////////////////////////////

Param::ParamEntry* pe_ptr = nullptr;
Param::ParamEntry* pe_nullPointer = nullptr;
START_SECTION(([Param::ParamEntry] ParamEntry()))
	pe_ptr = new Param::ParamEntry();
  TEST_NOT_EQUAL(pe_ptr,pe_nullPointer)
END_SECTION

START_SECTION(([Param::ParamEntry] ~ParamEntry()))
	delete pe_ptr;
END_SECTION

START_SECTION(([Param::ParamEntry] ParamEntry(const std::string &n, const DataValue &v, const std::string &d, const std::vector<std::string> &t=std::vector<std::string>())))
	Param::ParamEntry pe("n","v","d",{"advanced"});
	TEST_EQUAL(pe.name,"n")
	TEST_EQUAL(pe.description,"d")
	TEST_EQUAL(pe.value,"v")
	TEST_EQUAL(pe.tags.count("advanced")==1,true)

	 pe = Param::ParamEntry("n1","v1","d1");
	TEST_EQUAL(pe.name,"n1")
	TEST_EQUAL(pe.description,"d1")
	TEST_EQUAL(pe.value,"v1")
	TEST_EQUAL(pe.tags.count("advanced")==1,false)
END_SECTION

START_SECTION(([Param::ParamEntry] bool isValid(std::string& message) const))

	Param p;
	std::string m;
	p.setValue("int",5);
	TEST_EQUAL(p.getEntry("int").isValid(m),true);
	p.setMinInt("int",5);
	TEST_EQUAL(p.getEntry("int").isValid(m),true);
	p.setMaxInt("int",8);
	TEST_EQUAL(p.getEntry("int").isValid(m),true);
	p.setValue("int",10);
	TEST_EQUAL(p.getEntry("int").isValid(m),false);

	p.setValue("float",5.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setMinFloat("float",5.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setMaxFloat("float",8.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setValue("float",10.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),false);

	p.setValue("float",5.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setMinFloat("float",5.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setMaxFloat("float",8.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),true);
	p.setValue("float",10.1);
	TEST_EQUAL(p.getEntry("float").isValid(m),false);


	vector<std::string> strings;
	strings.push_back("bla");
	strings.push_back("bluff");
	p.setValue("string","bli");
	TEST_EQUAL(p.getEntry("string").isValid(m),true);
	p.setValidStrings("string",strings);
	TEST_EQUAL(p.getEntry("string").isValid(m),false);

	p.setValue("string_2","bla");
	TEST_EQUAL(p.getEntry("string_2").isValid(m),true);
	p.setValidStrings("string_2",strings);
	TEST_EQUAL(p.getEntry("string_2").isValid(m),true);

END_SECTION

START_SECTION(([Param::ParamEntry] bool operator==(const ParamEntry& rhs) const))
	Param::ParamEntry n1("n","d","v",{"advanced"});
	Param::ParamEntry n2("n","d","v",{"advanced"});

	TEST_EQUAL(n1==n2,true)

	n2.name = "name";
	TEST_EQUAL(n1==n2,false)
	n2 = n1;

	n2.value = "bla";
	TEST_EQUAL(n1==n2,false)
	n2 = n1;

	n2.description = "bla";
	TEST_EQUAL(n1==n2,true)

	n2.tags.clear();
	TEST_EQUAL(n1==n2,true)
END_SECTION

////////////////// Param::ParamNode ////////////////////////////////
////////////////////////////////////////////////////////////////////

Param::ParamNode* pn_ptr = nullptr;
Param::ParamNode* pn_nullPointer = nullptr;
START_SECTION(([Param::ParamNode] ParamNode()))
	pn_ptr = new Param::ParamNode();
  TEST_NOT_EQUAL(pn_ptr,pn_nullPointer)
END_SECTION

START_SECTION(([Param::ParamNode] ~ParamNode()))
	delete pn_ptr;
END_SECTION

START_SECTION(([Param::ParamNode] ParamNode(const std::string& n, const std::string& d)))
	Param::ParamNode n("n","d");
	TEST_EQUAL(n.name,"n")
	TEST_EQUAL(n.description,"d")

	n = Param::ParamNode("n1","d1");
	TEST_EQUAL(n.name,"n1")
	TEST_EQUAL(n.description,"d1")
END_SECTION

START_SECTION(([Param::ParamNode] bool operator==(const ParamNode& rhs) const))
	Param::ParamNode n1("n","d");
	Param::ParamNode n2("n","d");

	TEST_EQUAL(n1==n2,true)

	n2.name = "name";
	TEST_EQUAL(n1==n2,false)
	n2 = n1;

	n2.description = "bla";
	TEST_EQUAL(n1==n2,true)
	n2 = n1;

	n2.nodes.resize(5);
	TEST_EQUAL(n1==n2,false)
	n2 = n1;

	n2.entries.resize(5);
	TEST_EQUAL(n1==n2,false)
	n2 = n1;

	n2.entries.push_back(Param::ParamEntry("a","x",""));
	n2.entries.push_back(Param::ParamEntry("b","y",""));
	n1.entries.push_back(Param::ParamEntry("b","y",""));
	n1.entries.push_back(Param::ParamEntry("a","x",""));
	TEST_EQUAL(n1==n2,true)

	n2.nodes.push_back(Param::ParamNode("a","x"));
	n2.nodes.push_back(Param::ParamNode("b","y"));
	n1.nodes.push_back(Param::ParamNode("b","y"));
	n1.nodes.push_back(Param::ParamNode("a","x"));
	TEST_EQUAL(n1==n2,true)
END_SECTION

START_SECTION(([Param::ParamNode] std::string suffix(const std::string &key) const ))
	Param::ParamNode node;
	TEST_EQUAL(node.suffix(""),"")
	TEST_EQUAL(node.suffix("A"),"A")
	TEST_EQUAL(node.suffix("A:A"),"A")
	TEST_EQUAL(node.suffix("A:AB"),"AB")
	TEST_EQUAL(node.suffix("AB:A"),"A")
	TEST_EQUAL(node.suffix(":A"),"A")
END_SECTION

//Dummy Tree:
// A
// |-B(1)
// |-C
// | |-D(2)
// | |-E(3)
// |-B
//   |-G(4)
Param::ParamNode pn,n;
Param::ParamEntry e;
pn.name="A";
e.name="B"; e.value=1; pn.entries.push_back(e);
n.name="C"; pn.nodes.push_back(n);
e.name="D"; e.value=1; pn.nodes[0].entries.push_back(e);
e.name="E"; e.value=1; pn.nodes[0].entries.push_back(e);
n.name="B"; pn.nodes.push_back(n);
e.name="G"; e.value=1; pn.nodes[1].entries.push_back(e);


START_SECTION(([Param::ParamNode] Size size() const ))
	TEST_EQUAL(pn.size(),4)
	TEST_EQUAL(pn.nodes[0].size(),2)
	TEST_EQUAL(pn.nodes[1].size(),1)
END_SECTION

START_SECTION(([Param::ParamNode] EntryIterator findEntry(const std::string& name)))
	TEST_EQUAL(pn.findEntry("A")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("B")!=pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("C")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("D")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("E")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("F")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("G")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("H")==pn.entries.end(),true)
END_SECTION

START_SECTION(([Param::ParamNode] NodeIterator findNode(const std::string& name)))
	TEST_EQUAL(pn.findNode("A")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("B")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("C")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("D")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("E")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("F")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("G")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("H")==pn.nodes.end(),true)
END_SECTION

START_SECTION(([Param::ParamNode] ParamNode* findParentOf(const std::string &name)))
  TEST_EQUAL(pn.findParentOf("A"),pn_nullPointer)
	TEST_EQUAL(pn.findParentOf("B"),&pn)
	TEST_EQUAL(pn.findParentOf("C"),&pn)
	TEST_EQUAL(pn.findParentOf("C:D"),&(pn.nodes[0]))
	TEST_EQUAL(pn.findParentOf("C:E"),&(pn.nodes[0]))
  TEST_EQUAL(pn.findParentOf("F"),pn_nullPointer)
	TEST_EQUAL(pn.findParentOf("B:G"),&(pn.nodes[1]))
  TEST_EQUAL(pn.findParentOf("X"),pn_nullPointer)
  TEST_EQUAL(pn.findParentOf("H:X"),pn_nullPointer)
  TEST_EQUAL(pn.findParentOf("H:C:X"),pn_nullPointer)
  TEST_EQUAL(pn.findParentOf("H:C:"),pn_nullPointer)
END_SECTION

START_SECTION(([Param::ParamNode] ParamEntry* findEntryRecursive(const std::string& name)))
	TEST_EQUAL(pn.findEntryRecursive("A"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("B"),&(pn.entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("C"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("C:D"),&(pn.nodes[0].entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("C:E"),&(pn.nodes[0].entries[1]))
	TEST_EQUAL(pn.findEntryRecursive("F"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("B:G"),&(pn.nodes[1].entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("X"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("H:X"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("H:C:X"),pe_nullPointer)
	TEST_EQUAL(pn.findEntryRecursive("H:C:"),pe_nullPointer)
END_SECTION

//Dummy Tree:
// A
// |-B(1)
// |-C
// | |-D(2)
// | |-E(3)
// |-B
// | |-G(4)
// |-F
//   |-H(5)

START_SECTION(([Param::ParamNode] void insert(const ParamNode& node, const std::string& prefix = "")))
	Param::ParamNode node("","");
	node.entries.push_back(Param::ParamEntry("H",5,"",{"advanced"}));
	pn.insert(node,"F");
  TEST_NOT_EQUAL(pn.findEntryRecursive("F:H"),pe_nullPointer)

	pn.insert(node,"F:Z");
  TEST_NOT_EQUAL(pn.findEntryRecursive("F:Z:H"),pe_nullPointer)

	pn.insert(node,"F:Z:");
  TEST_NOT_EQUAL(pn.findEntryRecursive("F:Z::H"),pe_nullPointer)

	pn.insert(node,"FD:ZD:D");
  TEST_NOT_EQUAL(pn.findEntryRecursive("FD:ZD:D:H"),pe_nullPointer)

	node.name = "W";
	pn.insert(node);
  TEST_NOT_EQUAL(pn.findEntryRecursive("W:H"),pe_nullPointer)

	pn.insert(node,"Q");
  TEST_NOT_EQUAL(pn.findEntryRecursive("QW:H"),pe_nullPointer)
END_SECTION

START_SECTION(([Param::ParamNode] void insert(const ParamEntry& entry, const std::string& prefix = "")))
	Param::ParamEntry entry("H","","5",{"advanced"});

	pn.insert(entry);
  TEST_NOT_EQUAL(pn.findEntryRecursive("H"),pe_nullPointer)

	pn.insert(entry,"F");
  TEST_NOT_EQUAL(pn.findEntryRecursive("FH"),pe_nullPointer)

	pn.insert(entry,"G:");
  TEST_NOT_EQUAL(pn.findEntryRecursive("G:H"),pe_nullPointer)

	pn.insert(entry,"FD:ZD:D");
  TEST_NOT_EQUAL(pn.findEntryRecursive("FD:ZD:DH"),pe_nullPointer)
END_SECTION


////////////////// Param::ParamIterator ////////////////////////////
////////////////////////////////////////////////////////////////////


Param::ParamIterator* pi_ptr = nullptr;
Param::ParamIterator* pi_nullPointer = nullptr;
START_SECTION(([Param::ParamIterator] ParamIterator()))
	pi_ptr = new Param::ParamIterator();
  TEST_NOT_EQUAL(pi_ptr,pi_nullPointer)
END_SECTION

START_SECTION(([Param::ParamIterator] ~ParamIterator()))
	delete(pi_ptr);
END_SECTION

START_SECTION(([Param::ParamIterator] ParamIterator(const Param::ParamNode& root)))
	Param::ParamNode node;
	pi_ptr = new Param::ParamIterator(node);
  TEST_NOT_EQUAL(pi_ptr,pi_nullPointer)
  delete pi_ptr;
END_SECTION

START_SECTION(([Param::ParamIterator] const Param::ParamEntry& operator*()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",{"advanced"}));
	Param::ParamIterator it(node);
	TEST_EQUAL((*it).name,"name")
	TEST_EQUAL((*it).value,"value");
	TEST_EQUAL((*it).description,"description")
	TEST_EQUAL((*it).tags.count("advanced")==1,true)
END_SECTION

START_SECTION(([Param::ParamIterator] const Param::ParamEntry* operator->()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",{"advanced"}));
	Param::ParamIterator it(node);
	TEST_EQUAL(it->name,"name");
	TEST_EQUAL(it->value,"value");
	TEST_EQUAL(it->description,"description");
	TEST_EQUAL(it->tags.count("advanced")==1,true);
END_SECTION

//complicated subtree
// Root
//  |-A=1
//  |-R
//	| |-S
//  | | |-B=2
//  | | |-C=3
//  | | 
//  | |-U (empty)
//  |
//  |-T
//    |-D=4
Param::ParamNode root, r, s, t, u;
root.name="root";
r.name="r";
s.name="s";
t.name="t";
root.entries.push_back(Param::ParamEntry("A","1",""));
s.entries.push_back(Param::ParamEntry("B","2",""));
s.description="s_desc";
s.entries.push_back(Param::ParamEntry("C","3",""));
t.entries.push_back(Param::ParamEntry("D","4",""));
r.nodes.push_back(s);
u.description="empty";
r.nodes.push_back(u);
root.nodes.push_back(r);
root.nodes.push_back(t);

START_SECTION(([Param::ParamIterator] ParamIterator& operator++()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",{"advanced"}));
	node.entries.push_back(Param::ParamEntry("name2","value2","description2"));
	node.entries.push_back(Param::ParamEntry("name3","value3","description3",{"advanced"}));

	//linear list
	Param::ParamIterator it(node);
	++it;
	TEST_EQUAL(it->name,"name2");
	TEST_EQUAL(it->value,"value2");
	TEST_EQUAL(it->description,"description2");
	TEST_EQUAL(it->tags.count("advanced")==1,false);

	++it;
	TEST_EQUAL(it->name,"name3");
	TEST_EQUAL(it->value,"value3");
	TEST_EQUAL(it->description,"description3");
	TEST_EQUAL(it->tags.count("advanced")==1,true);

	++it;

	//subtree
	node.name = "root";
	node.nodes.push_back(node);
	node.nodes[0].name = "tree";
	node.nodes[0].entries[0].name = "name4";
	node.nodes[0].entries[1].name = "name5";
	node.nodes[0].entries[2].name = "name6";

	it = Param::ParamIterator(node);
	TEST_EQUAL(it->name,"name");
	TEST_EQUAL(it->value,"value");
	TEST_EQUAL(it->description,"description");
	TEST_EQUAL(it->tags.count("advanced")==1,true);

	++it;
	TEST_EQUAL(it->name,"name2");
	TEST_EQUAL(it->value,"value2");
	TEST_EQUAL(it->description,"description2");
	TEST_EQUAL(it->tags.count("advanced")==1,false);

	++it;
	TEST_EQUAL(it->name,"name3");
	TEST_EQUAL(it->value,"value3");
	TEST_EQUAL(it->description,"description3");
	TEST_EQUAL(it->tags.count("advanced")==1,true);

	++it;
	TEST_EQUAL(it->name,"name4");
	TEST_EQUAL(it->value,"value");
	TEST_EQUAL(it->description,"description");
	TEST_EQUAL(it->tags.count("advanced")==1,true);

	++it;
	TEST_EQUAL(it->name,"name5");
	TEST_EQUAL(it->value,"value2");
	TEST_EQUAL(it->description,"description2");
	TEST_EQUAL(it->tags.count("advanced")==1,false);

	++it;
	TEST_EQUAL(it->name,"name6");
	TEST_EQUAL(it->value,"value3");
	TEST_EQUAL(it->description,"description3");
	TEST_EQUAL(it->tags.count("advanced")==1,true);

	++it;

	//complicated subtree
	Param::ParamIterator it2(root);

	TEST_EQUAL(it2->name,"A");
	TEST_EQUAL(it2->value,"1");
	++it2;

	TEST_EQUAL(it2->name,"B");
	TEST_EQUAL(it2->value,"2");
	++it2;

	TEST_EQUAL(it2->name,"C");
	TEST_EQUAL(it2->value,"3");
	++it2;

	TEST_EQUAL(it2->name,"D");
	TEST_EQUAL(it2->value,"4");
	++it2;
END_SECTION

START_SECTION(([Param::ParamIterator] ParamIterator operator++(int)))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",{"advanced"}));
	node.entries.push_back(Param::ParamEntry("name2","value2","description2"));
	node.entries.push_back(Param::ParamEntry("name3","value3","description3",{"advanced"}));

	//linear list
	Param::ParamIterator it(node), it2(node);

	it2 = it++;
	TEST_EQUAL(it->name,"name2");
	TEST_EQUAL(it->value,"value2");
	TEST_EQUAL(it->description,"description2");
	TEST_EQUAL(it->tags.count("advanced")==1,false);
	TEST_EQUAL(it2->name,"name");
	TEST_EQUAL(it2->value,"value");
	TEST_EQUAL(it2->description,"description");
	TEST_EQUAL(it2->tags.count("advanced")==1,true);
END_SECTION

START_SECTION(([Param::ParamIterator] std::string getName() const))
	Param::ParamIterator it(root);

	TEST_EQUAL(it.getName(),"A");
	++it;

	TEST_EQUAL(it.getName(),"r:s:B");
	++it;

	TEST_EQUAL(it.getName(),"r:s:C");
	++it;

	TEST_EQUAL(it.getName(),"t:D");
	++it;
END_SECTION


START_SECTION(([Param::ParamIterator] bool operator==(const ParamIterator& rhs) const))
	Param::ParamIterator begin(root), begin2(root), end;
	TEST_EQUAL(begin==end, false)
	TEST_TRUE(begin == begin)
	TEST_TRUE(begin == begin2)
	TEST_TRUE(end == end)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin2;
	TEST_TRUE(begin == begin2)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin2;
	TEST_TRUE(begin == begin2)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin2;
	TEST_TRUE(begin == begin2)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_TRUE(begin == end)
	TEST_EQUAL(begin2==end, false)

	++begin2;
	TEST_TRUE(begin == begin2)
	TEST_TRUE(begin == end)
	TEST_TRUE(begin2 == end)
END_SECTION

START_SECTION(([Param::ParamIterator] bool operator!=(const ParamIterator& rhs) const))
	Param::ParamIterator begin(root), begin2(root), end;
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	TEST_TRUE(begin == begin2)
	TEST_TRUE(begin == begin)
	TEST_TRUE(begin2 == begin2)
	TEST_TRUE(end == end)
END_SECTION


START_SECTION(([Param::ParamIterator] const std::vector< TraceInfo>& getTrace() const))

	//Recap:
	//complicated subtree
	// Root
	//  |-A=1
	//  |-R
	//	| |-S
	//  | | |-B=2
	//  | | |-C=3
	//  | | 
	//  | |-U (empty)
	//  |
	//  |-T
	//    |-D=4
	//A
	Param::ParamIterator it(root);
	TEST_EQUAL(it.getTrace().size(),0);
	++it;

	//r:s:B
	TEST_EQUAL(it.getTrace().size(),2);
	TEST_EQUAL(it.getTrace()[0].name,"r");
	TEST_EQUAL(it.getTrace()[0].opened,true);
	TEST_EQUAL(it.getTrace()[1].name,"s");
	TEST_EQUAL(it.getTrace()[1].opened,true);
	TEST_EQUAL(it.getTrace()[1].description,"s_desc");
	++it;

	//r:s:C
	TEST_EQUAL(it.getTrace().size(),0);
	++it;

	//t:D
	TEST_EQUAL(it.getTrace().size(),3);
	TEST_EQUAL(it.getTrace()[0].name,"s");
	TEST_EQUAL(it.getTrace()[0].opened,false);
	TEST_EQUAL(it.getTrace()[1].name,"r");
	TEST_EQUAL(it.getTrace()[1].opened,false);
	TEST_EQUAL(it.getTrace()[2].name,"t");
	TEST_EQUAL(it.getTrace()[2].opened,true);
	++it;

	//end()
	TEST_EQUAL(it.getTrace().size(),1);
	TEST_EQUAL(it.getTrace()[0].name,"t");
	TEST_EQUAL(it.getTrace()[0].opened,false);
END_SECTION

///////////////////////// Param ///////////////////////////////
///////////////////////////////////////////////////////////////



Param* d10_ptr = nullptr;
Param* d10_nullPointer = nullptr;
START_SECTION((Param()))
	d10_ptr = new Param();
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~Param()))
	delete d10_ptr;
END_SECTION

START_SECTION((bool exists(const std::string& key) const))
	Param p;
	TEST_EQUAL(p.exists(""), false)
	TEST_EQUAL(p.exists("key"), false)
	TEST_EQUAL(p.exists("key:value"), false)
END_SECTION

START_SECTION((bool hasSection(const std::string& key) const))
  Param p;
	p.addSection("test", "");
	p.addSection("test:test", "");
	TEST_EQUAL(p.hasSection("test"), true)
	TEST_EQUAL(p.hasSection("test:test"), true)
	TEST_EQUAL(p.hasSection("test:"), true)
	TEST_EQUAL(p.hasSection("test:test:"), true)
	TEST_EQUAL(p.hasSection("sectionThatDoesNotExist"), false)
	TEST_EQUAL(p.hasSection("AnotherSection:"), false)
	TEST_EQUAL(p.hasSection("Section:WithinSection"), false)
	TEST_EQUAL(p.hasSection("Section:WithinSection:"), false)
END_SECTION

START_SECTION((const DataValue& getValue(const std::string &key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue("key:value"))
END_SECTION

START_SECTION((const std::string& getSectionDescription(const std::string& key) const))
	Param p;
	TEST_EQUAL(p.getSectionDescription(""),"")
	TEST_EQUAL(p.getSectionDescription("key"),"")
	TEST_EQUAL(p.getSectionDescription("key:value"),"")
END_SECTION

START_SECTION((const std::string& getDescription(const std::string &key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key:value"))
END_SECTION

START_SECTION((const ParamEntry& getEntry(const std::string &key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key:value"))
END_SECTION

START_SECTION((void setValue(const std::string &key, const DataValue& value, const std::string &description="", const std::stringList &tags=std::stringList())))
	Param p;
	p.setValue("key","value");
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), "value")
	TEST_EQUAL(p.getDescription("key"), "")
	TEST_EQUAL(p.hasTag("key","advanced"), false)

	p.setValue("key","value","description",{"advanced"});
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), "value")
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)

	p.setValue("key:key","value2","description2");
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), "value")
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
	TEST_EQUAL(p.exists("key:key"), true)
	TEST_EQUAL(p.getValue("key:key"), "value2")
	TEST_EQUAL(p.getDescription("key:key"), "description2")
	TEST_EQUAL(p.hasTag("key:key","advanced"), false)
END_SECTION

START_SECTION((std::vector<std::string> getTags(const std::string& key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getTags("key"))

	p.setValue("key","value");
	TEST_EQUAL(p.getTags("key").size(),0)
END_SECTION

START_SECTION((void addTag(const std::string& key, const std::string& tag)))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTag("key","bla"))
	std::vector<std::string> error_list;
	error_list.push_back("a,b");
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTags("key",error_list))

	p.setValue("key","value");
	TEST_EQUAL(p.getTags("key").size(),0)
	p.addTag("key","advanced");
	TEST_EQUAL(p.getTags("key").size(),1)
	p.addTag("key","advanced");
	TEST_EQUAL(p.getTags("key").size(),1)
	p.addTag("key","advanced2");
	TEST_EQUAL(p.getTags("key").size(),2)
END_SECTION

START_SECTION((bool hasTag(const std::string& key, const std::string& tag) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.hasTag("key","bla"))

	p.setValue("key","value");
	TEST_EQUAL(p.hasTag("key","advanced"),false)
	TEST_EQUAL(p.hasTag("key","advanced2"),false)
	p.addTag("key","advanced");
	TEST_EQUAL(p.hasTag("key","advanced"),true)
	TEST_EQUAL(p.hasTag("key","advanced2"),false)
	p.addTag("key","advanced2");
	TEST_EQUAL(p.hasTag("key","advanced"),true)
	TEST_EQUAL(p.hasTag("key","advanced2"),true)
END_SECTION

START_SECTION((void addTags(const std::string& key, const std::vector<std::string>& tags)))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTags("key",std::vector<std::string>()))
	std::vector<std::string> error_list;
	error_list.push_back("a,b");
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTags("key",error_list))

	p.setValue("key","value");
	TEST_EQUAL(p.hasTag("key","advanced"),false)
	TEST_EQUAL(p.hasTag("key","advanced2"),false)
	p.addTags("key",{"advanced","advanced2"});
	TEST_EQUAL(p.hasTag("key","advanced"),true)
	TEST_EQUAL(p.hasTag("key","advanced2"),true)
END_SECTION

START_SECTION((void clearTags(const std::string& key)))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.clearTags("key"))
	p.setValue("key","value");
	p.addTag("key","advanced");
	TEST_EQUAL(p.getTags("key").size(),1)
	p.clearTags("key");
	TEST_EQUAL(p.getTags("key").size(),0)
END_SECTION

START_SECTION((bool empty() const))
	Param p;
	TEST_EQUAL(p.empty(), true)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.empty(), false)

	Param p2;
	TEST_EQUAL(p2.empty(), true)
	p2.setValue("a:key",17.4f);
	TEST_EQUAL(p2.empty(), false)
END_SECTION

START_SECTION((void clear()))
	Param p;
	p.setValue("key",17.4,"keydesc");
	p.clear();
	TEST_EQUAL(p.empty(), true)

	Param p2;
	p2.setValue("a:b:key",17.4,"keydesc");
	p2.clear();
	TEST_EQUAL(p2.empty(), true)
END_SECTION

START_SECTION((Size size() const))
	Param p;
	TEST_EQUAL(p.size(), 0)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.size(), 1)
	p.setValue("key",17.4f);
	TEST_EQUAL(p.size(), 1)
	p.setValue("key:a",17.5f);
	TEST_EQUAL(p.size(), 2)
	p.setValue("key:a",18.5f);
	TEST_EQUAL(p.size(), 2)
	p.setValue("key:b",18.5f);
	TEST_EQUAL(p.size(), 3)
	p.setValue("b",18.5f);
	TEST_EQUAL(p.size(), 4)
END_SECTION

START_SECTION((void setSectionDescription(const std::string &key, const std::string &description)))
	Param p;

	p.setValue("test:test",47.1);
	p.setValue("test2:test",47.1);
	p.setValue("test:test2:test",47.1);
	p.setValue("test:test:test",47.1);
	p.setSectionDescription("test","a");
	p.setSectionDescription("test2","b");
	p.setSectionDescription("test:test","c");
	p.setSectionDescription("test:test2","d");
	TEST_EQUAL(p.getSectionDescription("test"), "a")
	TEST_EQUAL(p.getSectionDescription("test2"), "b")
	TEST_EQUAL(p.getSectionDescription("test:test"), "c")
	TEST_EQUAL(p.getSectionDescription("test:test2"), "d")
END_SECTION

START_SECTION([EXTRA](friend std::ostream& operator << (std::ostream& os, const Param& param)))
	Param p;
	p.setValue("key", 17.5);
	stringstream ss;
	ss << p;
	TEST_EQUAL(ss.str(), "\"key\" -> \"17.5\"\n")

	ss.str("");
	p.setValue("key", 17.5, "thiskey");
	ss<<p;
	TEST_EQUAL(ss.str(), "\"key\" -> \"17.5\" (thiskey)\n")

	ss.str("");
	p.clear();
	p.setValue("tree:key", 17.5);
	ss<<p;
	TEST_EQUAL(ss.str(), "\"tree|key\" -> \"17.5\"\n")
END_SECTION

START_SECTION((void insert(const std::string& prefix, const Param &param)))
	Param p;
	p.setValue("a",17,"intdesc");
	p.setValue("n1:b",17.4f,"floatdesc");
	p.setValue("n1:c","test,test,test","stringdesc");
	p.setValue("n2:d",17.5f);
	p.setSectionDescription("n1","sectiondesc");

	Param p2;

	p2.insert("prefix",p);
	TEST_EQUAL(p2.size(),4)
	TEST_EQUAL(Int(p2.getValue("prefixa")), 17)
	TEST_STRING_EQUAL(p2.getDescription("prefixa"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("prefixn1:b")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("prefixn1:b"), "floatdesc")
	TEST_EQUAL(p2.getValue("prefixn1:c"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("prefixn1:c"), "stringdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("prefixn2:d")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("prefixn2:d"), "")
	TEST_EQUAL(p2.getSectionDescription("prefixn1"),"sectiondesc")

	p2.insert("",p);
	TEST_EQUAL(p2.size(),8)
	TEST_EQUAL(Int(p2.getValue("a")), 17)
	TEST_STRING_EQUAL(p2.getDescription("a"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("n1:b")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("n1:b"), "floatdesc")
	TEST_EQUAL(p2.getValue("n1:c"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("n1:c"), "stringdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("n2:d")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("n2:d"), "")
	TEST_EQUAL(p2.getSectionDescription("n1"),"sectiondesc")

	p2.insert("n3:",p);
	TEST_EQUAL(p2.size(),12)
	TEST_EQUAL(Int(p2.getValue("n3:a")), 17)
	TEST_STRING_EQUAL(p2.getDescription("n3:a"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("n3:n1:b")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("n3:n1:b"), "floatdesc")
	TEST_EQUAL(p2.getValue("n3:n1:c"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("n3:n1:c"), "stringdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("n3:n2:d")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("n3:n2:d"), "")
	TEST_EQUAL(p2.getSectionDescription("n3:n1"),"sectiondesc")

	p.clear();
	p.setValue("a",18,"intdesc");
	p.setValue("n1:b",17.7f,"floatdesc");
	p.setValue("n1:c","test,test,test,test","stringdesc");
	p.setValue("n2:d",17.8f);

	p2.insert("",p);
	TEST_EQUAL(p2.size(),12)
	TEST_EQUAL(Int(p2.getValue("a")), 18)
	TEST_REAL_SIMILAR(float(p2.getValue("n1:b")), 17.7)
	TEST_EQUAL(p2.getValue("n1:c"), "test,test,test,test")
	TEST_REAL_SIMILAR(float(p2.getValue("n2:d")), 17.8)
END_SECTION

Param p_src;
p_src.setValue("test:float",17.4f,"floatdesc");
p_src.setValue("test:string","test,test,test","stringdesc");
p_src.setValue("test:int",17,"intdesc");
p_src.setValue("test2:float",17.5f);
p_src.setValue("test2:string","test2");
p_src.setValue("test2:int",18);
p_src.setSectionDescription("test","sectiondesc");
p_src.addTags("test:float", {"a", "b", "c"});

START_SECTION((Param(const Param& rhs)))
	Param p2(p_src);
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p_src.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p_src.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p_src.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), "")
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), "")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), "")
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
	TEST_EQUAL(p2.getTags("test:float").size(), 3)
	TEST_EQUAL(p2.getTags("test:float") == ListUtils::create<std::string>("a,b,c"), true)
END_SECTION

START_SECTION((Param& operator = (const Param& rhs)))
	Param p2;
	p2=p_src;
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p_src.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p_src.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), "")
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), "")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), "")
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
	TEST_EQUAL(p2.getTags("test:float").size(), 3)
	TEST_EQUAL(p2.getTags("test:float") == ListUtils::create<std::string>("a,b,c"), true)
END_SECTION

START_SECTION((Param copy(const std::string &prefix, bool remove_prefix=false) const))
	Param p2;

	p2 = p_src.copy("notthere:");
	TEST_EQUAL((p2.empty()),true)

	p2 = p_src.copy("test:");

	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EXCEPTION(Exception::ElementNotFound, p2.getValue("test2:float"))

	p2 = p_src.copy("test:",true);
	TEST_REAL_SIMILAR(float(p2.getValue("float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("float"), "floatdesc")
	TEST_EQUAL(p2.getValue("string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("string"), "stringdesc")

	p2 = p_src.copy("test");
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), "")
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), "")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), "")
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
END_SECTION

START_SECTION((void remove(const std::string& key)))

	Param p2(p_src);
	p2.setValue("test:string2","test,test");

	TEST_EQUAL(p2.size(),7)

	p2.remove("test");
	TEST_EQUAL(p2.size(),7)

	p2.remove("test2");
	TEST_EQUAL(p2.size(),7)

	p2.remove("test:strin");
	TEST_EQUAL(p2.size(),7)

	p2.remove("test:string");
	TEST_EQUAL(p2.size(),6)

	p2.remove("test:string2");
	TEST_EQUAL(p2.size(),5)

	p2.remove("test:float");
	TEST_EQUAL(p2.size(),4)

	p2.remove("test:int");
	TEST_EQUAL(p2.size(),3)

  // test deletion of nodes (when using a trailing ':')
  p2 = p_src;
	p2.setValue("test:string2","an entry");
  p2.setValue("test:string2:e1","subnode with entries");
  p2.setValue("test:string2:sn2","subsubnode with entries");
  p2.setValue("test:string2:sn2:e1","subsubnode with entries");
  p2.setValue("test:string2:sn2:e2","subsubnode with entries");

  Param p3 = p2;

  TEST_EQUAL(p2.size(),11)

  std::cout << "p2 is " << p2 << "\n";

  p2.remove("test:"); // test subtree removal
	TEST_EQUAL(p2.size(),3)


  p3.remove("test:string2:sn2:e2:"); // nothing should happen
  TEST_EQUAL(p3.size(),11)

  p3.remove("test:string2:sn2:e1");  // delete one, the parent node is still populated
  TEST_EQUAL(p3.size(),10)

  p3.remove("test:string2:sn2:e2");  // delete last entry in subnode sn2
  TEST_EQUAL(p3.size(),9)


END_SECTION

START_SECTION((void removeAll(const std::string& prefix)))
	Param p2(p_src);

	p2.removeAll("test:float");
	TEST_EXCEPTION(Exception::ElementNotFound, p2.getValue("test:float"))
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")

	p2.removeAll("test:");
	TEST_EXCEPTION(Exception::ElementNotFound, p2.getValue("test:string"))
	TEST_EXCEPTION(Exception::ElementNotFound, p2.getValue("test:int"))
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)

	p2.removeAll("test");
	TEST_EQUAL(p2.empty(),true)

	cout << p2;
END_SECTION


START_SECTION((bool operator == (const Param& rhs) const))
	Param p2(p_src);
	TEST_TRUE(p_src == p2)
	p2.setValue("test:float",17.5f);
	TEST_EQUAL(p_src==p2, false)
	p2 = p_src;
	p2.setValue("test:float3",17.4f);
	TEST_EQUAL(p_src==p2, false)
	p2 = p_src;
	p2.removeAll("test:float");
	TEST_EQUAL(p_src==p2, false)

	//it should be independent of entry order
	Param p3,p4;
	p3.setValue("1",1);
	p3.setValue("2",2);
	p4.setValue("2",2);
	p4.setValue("1",1);
	TEST_TRUE(p3 == p4)

	//it should be independent of node order
	Param p5,p6;
	p5.setValue("1:1",1);
	p5.setValue("2:1",1);
	p6.setValue("2:1",1);
	p6.setValue("1:1",1);
	TEST_TRUE(p5 == p6)

END_SECTION

START_SECTION((void setDefaults(const Param& defaults, const std::string& prefix="", bool showMessage=false)))
	Param defaults;
	defaults.setValue("float",1.0f,"float");
	defaults.setValue("float2",2.0f,"float2");
	defaults.setValue("string","default string1","string");
	defaults.setValue("string2","default string2","string2");
	defaults.setValue("PATH:onlyfordescription",45.2);

	defaults.setValue("stringlist",std::vector<std::string>{"a","b","c"},"stringlist");
	defaults.setValue("stringlist2",std::vector<std::string>{"d","e","f"},"stringlist2");
	defaults.setValue("intlist",ListUtils::create<Int>("1,2,3"),"intlist");
	defaults.setValue("intlist2",ListUtils::create<Int>("11,22,33"),"intlist2");
	defaults.setValue("doublelist",ListUtils::create<double>("1.2,2.3"),"doublelist");
	defaults.setValue("doublelist2",ListUtils::create<double>("11.22,22.33"),"doublelist2");
	defaults.setSectionDescription("PATH","PATHdesc");
	Param p2;
	p2.setValue("PATH:float",-1.0f,"PATH:float");
	p2.setValue("PATH:string","some string","PATH:string");
	p2.setValue("float",-2.0f,"float");
	p2.setValue("string","other string","string");

	p2.setValue("PATH:stringlist",std::vector<std::string>{"d","a","v","i","d"},"PATH:stringlist");
	p2.setValue("stringlist",std::vector<std::string>{"r","o","c","k","s"},"stringlist");
	p2.setValue("PATH:intlist2",ListUtils::create<Int>("14,9"),"PATH:intlist2");
	p2.setValue("intlist", ListUtils::create<Int>("16,9"),"intlist");
	p2.setValue("PATH:doublelist2",ListUtils::create<double>("6.66,6.16"),"PATH:doublelist2");
	p2.setValue("doublelist",ListUtils::create<double>("1.2,5.55"),"doublelist");

	TEST_EQUAL(p2.size(),10);

	p2.setDefaults(defaults);
	TEST_EQUAL(p2.size(),16);
	TEST_REAL_SIMILAR(float(p2.getValue("float")),-2.0);
	TEST_STRING_EQUAL(p2.getDescription("float"),"float");
	TEST_REAL_SIMILAR(float(p2.getValue("float2")),2.0);
	TEST_STRING_EQUAL(p2.getDescription("float2"),"float2");
	TEST_EQUAL(string(p2.getValue("string")),"other string");
	TEST_STRING_EQUAL(p2.getDescription("string"),"string");
	TEST_EQUAL(string(p2.getValue("string2")),"default string2");
	TEST_STRING_EQUAL(p2.getDescription("string2"),"string2");
	TEST_STRING_EQUAL(p2.getSectionDescription("PATH"),"PATHdesc");

	TEST_EQUAL(p2.getValue("stringlist") == ListUtils::create<std::string>("r,o,c,k,s"), true)
	TEST_EQUAL(p2.getValue("intlist") == ListUtils::create<Int>("16,9"), true)
	TEST_EQUAL(p2.getValue("doublelist") == ListUtils::create<double>("1.2,5.55"), true)
	TEST_EQUAL(p2.getValue("stringlist2") == ListUtils::create<std::string>("d,e,f"), true)
	TEST_EQUAL(p2.getValue("intlist2") == ListUtils::create<Int>("11,22,33"), true)
	TEST_EQUAL(p2.getValue("doublelist2") == ListUtils::create<double>("11.22,22.33"), true)



	p2.setDefaults(defaults,"PATH");

	TEST_EQUAL(p2.size(),22);
	TEST_REAL_SIMILAR(float(p2.getValue("PATH:float")),-1.0);
	TEST_STRING_EQUAL(p2.getDescription("PATH:float"),"PATH:float");
	TEST_REAL_SIMILAR(float(p2.getValue("PATH:float2")),2.0);
	TEST_STRING_EQUAL(p2.getDescription("PATH:float2"),"float2");
	TEST_EQUAL(string(p2.getValue("PATH:string")),"some string");
	TEST_STRING_EQUAL(p2.getDescription("PATH:string"),"PATH:string");
	TEST_EQUAL(string(p2.getValue("PATH:string2")),"default string2");
	TEST_STRING_EQUAL(p2.getDescription("PATH:string2"),"string2");
	TEST_STRING_EQUAL(p2.getSectionDescription("PATH"),"PATHdesc");
	TEST_STRING_EQUAL(p2.getSectionDescription("PATH:PATH"),"PATHdesc");

	TEST_EQUAL(p2.getValue("PATH:stringlist") == ListUtils::create<std::string>("d,a,v,i,d"), true)
	TEST_EQUAL(p2.getValue("PATH:intlist") == ListUtils::create<Int>("1,2,3"), true)
	TEST_EQUAL(p2.getValue("PATH:doublelist") == ListUtils::create<double>("1.2,2.3"), true)

END_SECTION

const char* a1 ="executable";
const char* a2 ="-a";
const char* a3 ="av";
const char* a4 ="-b";
const char* a5 ="bv";
const char* a6 ="-c";
const char* a7 ="cv";
const char* a8 ="rv1";
const char* a9 ="rv2";
const char* a10="-1.0";

const char* command_line[9]; // "executable -a av -b bv -c cv rv1 rv2"
command_line[0] = a1;
command_line[1] = a2;
command_line[2] = a3;
command_line[3] = a4;
command_line[4] = a5;
command_line[5] = a6;
command_line[6] = a7;
command_line[7] = a8;
command_line[8] = a9;

const char* command_line2[6]; // "executable -a av -b -c cv"
command_line2[0] = a1;
command_line2[1] = a2;
command_line2[2] = a3;
command_line2[3] = a4;
command_line2[4] = a6;
command_line2[5] = a7;

const char* command_line3[6]; // "executable -a -b -c cv rv1"
command_line3[0] = a1;
command_line3[1] = a2;
command_line3[2] = a4;
command_line3[3] = a6;
command_line3[4] = a7;
command_line3[5] = a8;

const char* command_line4[10]; // "executable -a -1.0 -b bv -c cv rv1 rv2 -1.0"
command_line4[0] = a1;
command_line4[1] = a2;
command_line4[2] = a10;
command_line4[3] = a4;
command_line4[4] = a5;
command_line4[5] = a6;
command_line4[6] = a7;
command_line4[7] = a8;
command_line4[8] = a9;
command_line4[9] = a10;

START_SECTION((void parseCommandLine(const int argc, const char **argv, const std::string& prefix="")))
	Param p2,p3;
	p2.parseCommandLine(9,command_line,"test4");
	p3.setValue("test4:-a","av");
	p3.setValue("test4:-b","bv");
	p3.setValue("test4:-c","cv");
	p3.setValue("test4:misc",std::vector<std::string>{"rv1","rv2"});
	TEST_EQUAL(p2==p3,true)

	Param p20,p30;
	p20.parseCommandLine(6,command_line2);
	p30.setValue("-a","av");
	p30.setValue("-b","");
	p30.setValue("-c","cv");
	TEST_EQUAL(p20==p30,true)

	Param p200,p300;
	p200.parseCommandLine(10,command_line4,"test4");
	p300.setValue("test4:-a","-1.0");
	p300.setValue("test4:-b","bv");
	p300.setValue("test4:-c","cv");
	p300.setValue("test4:misc",std::vector<std::string>{"rv1","rv2","-1.0"});
	TEST_EQUAL(p200==p300,true)

END_SECTION

const char* m1 ="mult";
const char* m2 ="-d";
const char* m3 ="1.333";
const char* m4 ="2.23";
const char* m5 ="3";
const char* m6 ="-e";
const char* m7 ="4";
const char* m8 ="-f";
const char* m9 ="-g";

const char* command_line_mult[9];	// "mult -d 1.333 2.23 3 -e 4 -f -g"
command_line_mult[0] = m1;
command_line_mult[1] = m2;
command_line_mult[2] = m3;
command_line_mult[3] = m4;
command_line_mult[4] = m5;
command_line_mult[5] = m6;
command_line_mult[6] = m7;
command_line_mult[7] = m8;
command_line_mult[8] = m9;

START_SECTION((void parseCommandLine(const int argc, const char **argv, const std::map< std::string, std::string > &options_with_one_argument, const std::map< std::string, std::string > &options_without_argument, const std::map< std::string, std::string > &options_with_multiple_argument, const std::string &misc="misc", const std::string &unknown="unknown")))

	std::map<std::string,std::string> with_one,without,with_multiple;
	with_one["-a"]="a";
	with_one["-b"]="b";
	with_one["-c"]="c";

	with_multiple["-d"] = "d";
	with_multiple["-e"] = "e";
	with_multiple["-f"] = "f";
	with_multiple["-g"] = "g";

	Param p2,p3;
	p2.parseCommandLine(10,command_line4,with_one,without,with_multiple,"misc_","unknown_");
	p3.setValue("a","-1.0");
	p3.setValue("b","bv");
	p3.setValue("c","cv");
	p3.setValue("misc_",std::vector<std::string>{"rv1","rv2","-1.0"});
	TEST_EQUAL(p2==p3,true)

	Param p4,p5;
	p4.parseCommandLine(9,command_line,with_one,without,with_multiple,"misc_","unknown_");
	p5.setValue("a","av");
	p5.setValue("b","bv");
	p5.setValue("c","cv");
	p5.setValue("misc_",std::vector<std::string>{"rv1","rv2"});
	TEST_EQUAL(p4==p5,true)

	with_one.clear();
	with_one["-a"]="a";
	without["-b"]="b";

	Param p40,p50;
	p40.parseCommandLine(9,command_line,with_one,without,with_multiple,"misc__","unknown__");
	p50.setValue("a","av");
	p50.setValue("b","true");
	p50.setValue("misc__",std::vector<std::string>{"bv","cv","rv1","rv2"});
	p50.setValue("unknown__",std::vector<std::string>{"-c"});
	TEST_EQUAL(p40==p50,true)
	TEST_EQUAL(p40,p50)
	//"executable -a av -b -c cv"
	Param p400,p500;
	p400.parseCommandLine(6,command_line2,with_one,without,with_multiple,"misc__","unknown__");
	p500.setValue("a","av");
	p500.setValue("b","true");
	p500.setValue("misc__",std::vector<std::string>{"cv"});
	p500.setValue("unknown__",std::vector<std::string>{"-c"});
	TEST_EQUAL(p400==p500,true)

	//"executable -a -b -c cv rv1"
	Param p4000,p5000;
	p4000.parseCommandLine(6,command_line3,with_one,without,with_multiple,"misc__","unknown__");
	p5000.setValue("a","");
	p5000.setValue("b","true");
	p5000.setValue("misc__",std::vector<std::string>{"cv","rv1"});
	p5000.setValue("unknown__",std::vector<std::string>{"-c"});
	TEST_EQUAL(p4000==p5000,true)

	// list options:
	Param p6,p7;
	p6.parseCommandLine(9,command_line_mult,with_one,without,with_multiple,"misc__","unkown__");
	p7.setValue("d",std::vector<std::string>{"1.333","2.23","3"});
	p7.setValue("e",std::vector<std::string>{"4"});
	p7.setValue("f",std::vector<std::string>());
	p7.setValue("g",std::vector<std::string>());
	TEST_EQUAL(p6,p7);

	Param p8,p9;
	p9.parseCommandLine(4,command_line_mult,with_one,without,with_multiple,"misc__","unkown__");
	p8.setValue("d", std::vector<std::string>{"1.333","2.23"});
	TEST_EQUAL(p9,p8);

END_SECTION

START_SECTION((void update(const Param& old_version, const bool add_unknown, Logger::LogStream& stream)))
	Param common;
	common.setValue("float",1.0f,"float");
	common.setValue("float2",2.0f,"float2");
	common.setValue("string","default string1","string");
	common.setValue("string2","default string2","string2");
	common.setValue("PATH:onlyfordescription",45.2);

	common.setValue("stringlist",std::vector<std::string>{"a","b","c"},"stringlist");
	common.setValue("stringlist2",std::vector<std::string>{"d","e","f"},"stringlist2");
	common.setValue("intlist",ListUtils::create<Int>("1,2,3"),"intlist");

  // copy and alter
  Param old = common;
  //old.setValue("recently_removed_float",1.1f,"float");  // should not make it into new param
  old.setValue("old_type","a string","string");
  old.setValue("some:version","1.2","old version");
  old.setValue("some:1:type","unlabeled","type");
  old.setValue("some:type","unlabeled","type");
	old.setValue("stringlist2",std::vector<std::string>{"d","e","f","altered"},"stringlist2"); // change some values, we expect them to show up after update()
	old.setValue("intlist",ListUtils::create<Int>("3"),"intlist");

  Param defaults = common;
  defaults.setValue("old_type",3,"old_type has evolved from string to int"); // as type has changed, this value should be kept
  defaults.setValue("some:version","1.9","new version"); // this value should be kept (due to its reserved name)
  defaults.setValue("some:1:type","information","type");   // this value should be kept (due to its reserved name at depth 2)
  defaults.setValue("some:type","information","type");   // this value should NOT be kept (wrong depth)
  defaults.setValue("new_value",3,"new param not present in old");

  Param expected = defaults;
	expected.setValue("stringlist2",std::vector<std::string>{"d","e","f","altered"},"stringlist2"); // change some values, we expect them to show up after update()
	expected.setValue("intlist",ListUtils::create<Int>("3"),"intlist");
  expected.setValue("some:type","unlabeled","type");

  // update()
  defaults.update(old);

  TEST_EQUAL(defaults,expected);
END_SECTION

START_SECTION((void merge(const Param& toMerge)))
{
  Param original;
  original.setValue("a", 2.0f, "a value");
  original.setMinFloat("a", 0.0f);
  original.setValue("b", "value", "b value");

  Param toMerge;
  toMerge.setValue("b", "value", "a value");
  toMerge.setValue("section:a", "a-value", "section:a");
  toMerge.setSectionDescription("section", "section description");
  toMerge.setValue("section:b", "b-value", "section:b");

  Param expected;
  expected.setValue("a", 2.0f, "a value");
  expected.setMinFloat("a", 0.0f);
  expected.setValue("b", "value", "b value");
  expected.setValue("section:a", "a-value", "section:a");
  expected.setValue("section:b", "b-value", "section:b");
  expected.setSectionDescription("section", "section description");

  original.merge(toMerge);
  TEST_EQUAL(original, expected)
  TEST_EQUAL(original.getSectionDescription("section"),expected.getSectionDescription("section"))

  Param p1;
  p1.setValue("in", "in-value", "in-description");
  p1.setValue("out", "out-value", "out-description");
	p1.setValue("reference:index", "reference:index value", "reference:index description");
  p1.setSectionDescription("reference", "reference description");
  p1.setValue("algorithm:sub_param", "algorithm:sub_param value", "algorithm:sub_param description");

  Param p2;
  p2.setValue("reference:index", "reference:index value", "reference:index description");
  p2.setSectionDescription("reference", "reference description");
  p2.setValue("algorithm:sub_param", "algorithm:sub_param value", "algorithm:sub_param description");
  p2.setValue("algorithm:superimposer:mz_pair_max_distance", "algorithm:superimposer:mz_pair_max_distance value", "algorithm:superimposer:mz_pair_max_distance description");
	p2.setSectionDescription("algorithm", "algorithm description");
	p2.setSectionDescription("algorithm:superimposer", "algorithm:superimposer description");

	Param expected_2;
  expected_2.setValue("in", "in-value", "in-description");
  expected_2.setValue("out", "out-value", "out-description");
  expected_2.setValue("algorithm:sub_param", "algorithm:sub_param value", "algorithm:sub_param description");
  expected_2.setValue("reference:index", "reference:index value", "reference:index description");
  expected_2.setSectionDescription("reference", "reference description");
  expected_2.setValue("algorithm:superimposer:mz_pair_max_distance", "algorithm:superimposer:mz_pair_max_distance value", "algorithm:superimposer:mz_pair_max_distance description");
	expected_2.setSectionDescription("algorithm", "algorithm description");
	expected_2.setSectionDescription("algorithm:superimposer", "algorithm:superimposer description");

  p1.merge(p2);
	TEST_EQUAL(p1, expected_2)
  TEST_EQUAL(p1.getSectionDescription("algorithm"),expected_2.getSectionDescription("algorithm"))
  TEST_EQUAL(p1.getSectionDescription("algorithm:superimposer"),expected_2.getSectionDescription("algorithm:superimposer"))
  TEST_EQUAL(p1.getSectionDescription("reference"),expected_2.getSectionDescription("reference"))
}
END_SECTION

START_SECTION((ParamIterator findFirst(const std::string &leaf) const ))
{
  Param p;
  p.setValue("a:b:leaf", "leaf_val1", "leaf 1");
  p.setValue("b:a:leaf", "leaf_val2", "leaf 2");
  p.setValue("a:c:leaf", "leaf_val3", "leaf 3");
  p.setValue("a:c:another-leaf", "leaf_val4", "leaf 3");

  Param::ParamIterator pI = p.findFirst("leaf");
  TEST_EQUAL(pI.getName(), "a:b:leaf")

  p.remove("a:b:leaf");
  pI = p.findFirst("leaf");
  TEST_EQUAL(pI.getName(), "a:c:leaf")

  p.remove("a:c:leaf");
  pI = p.findFirst("leaf");
  TEST_EQUAL(pI.getName(), "b:a:leaf")

  p.remove("b:a:leaf");
  pI = p.findFirst("leaf");
  TEST_EQUAL(pI == p.end(), true)
}
END_SECTION

START_SECTION((ParamIterator findNext(const std::string &leaf, const ParamIterator &start_leaf) const))
{
  Param p;
  p.setValue("a:b:leaf", "leaf_val1", "leaf 1");
  p.setValue("b:a:leaf", "leaf_val2", "leaf 2");
  p.setValue("a:c:leaf", "leaf_val3", "leaf 3");
  p.setValue("a:c:another-leaf", "leaf_val4", "leaf 3");

  Param::ParamIterator pI = p.findFirst("leaf");
  TEST_EQUAL(pI.getName(), "a:b:leaf")

  pI = p.findNext("leaf", pI);
  TEST_EQUAL(pI.getName(), "a:c:leaf")

  pI = p.findNext("leaf", pI);
  TEST_EQUAL(pI.getName(), "b:a:leaf")

  pI = p.findNext("leaf", pI);
  TEST_EQUAL(pI == p.end(), true)
}
END_SECTION

START_SECTION((ParamIterator begin() const))
        NOT_TESTABLE;
END_SECTION

START_SECTION((ParamIterator end() const))
	Param p;
	p.setValue("a",5);
	p.setValue("b:a",6);
	p.setValue("b:b",7);
	p.setValue("c",8);

	Param::ParamIterator it = p.begin();
	TEST_EQUAL(it->name, "a")
	TEST_EQUAL(it.getName(), "a")
	TEST_EQUAL((UInt)it->value, 5)

	++it;
	TEST_EQUAL(it->name, "c")
	TEST_EQUAL(it.getName(), "c")
	TEST_EQUAL((UInt)it->value, 8)

	++it;
	TEST_EQUAL(it->name, "a")
	TEST_EQUAL(it.getName(), "b:a")
	TEST_EQUAL((UInt)it->value, 6)

	++it;
	TEST_EQUAL(it->name, "b")
	TEST_EQUAL(it.getName(), "b:b")
	TEST_EQUAL((UInt)it->value, 7)

	++it;
	TEST_EQUAL(it==p.end(),true)
END_SECTION

START_SECTION((void setValidStrings(const std::string &key, const std::vector< std::string > &strings)))
  vector<std::string> strings;
  strings.push_back("bla");
  Param d;
  d.setValue("ok","string");
  d.setValue("dummy",5);

  d.setValidStrings("ok",strings);
  TEST_EQUAL(d.getEntry("ok").valid_strings==strings, true);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setValidStrings("dummy",strings))
  strings.push_back("sdf,sdfd");
  TEST_EXCEPTION(Exception::InvalidParameter, d.setValidStrings("ok",strings))
END_SECTION

START_SECTION((void setMinInt(const std::string &key, Int min)))
  Param d;
  d.setValue("ok",4);
  d.setValue("dummy",5.5);

  d.setMinInt("ok",4);
  TEST_EQUAL(d.getEntry("ok").min_int,4);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMinInt("dummy",4))
END_SECTION

START_SECTION((void setMaxInt(const std::string &key, Int max)))
  Param d;
  d.setValue("ok",4);
  d.setValue("dummy",5.5);

  d.setMaxInt("ok",4);
  TEST_EQUAL(d.getEntry("ok").max_int,4);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMaxInt("dummy",4))
END_SECTION

START_SECTION((void setMinFloat(const std::string &key, double min)))
  Param d;
  d.setValue("ok",4.5);
  d.setValue("dummy",4);

  d.setMinFloat("ok",4.0);
  TEST_REAL_SIMILAR(d.getEntry("ok").min_float,4.0);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMinFloat("dummy",4.5))
END_SECTION

START_SECTION((void setMaxFloat(const std::string &key, double max)))
  Param d;
  d.setValue("ok",4.5);
  d.setValue("dummy",4);

  d.setMaxFloat("ok",4.0);
  TEST_REAL_SIMILAR(d.getEntry("ok").max_float,4.0);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMaxFloat("dummy",4.5))
END_SECTION

// warnings for unknown parameters
// keep outside the scope of a single test to avoid destruction, leaving
// OpenMS_Log_warn in an undefined state
ostringstream os;
// checkDefaults sends its warnings to OPENMS_LOG_WARN so we register our own
// listener here to check the output
OpenMS_Log_warn.remove(cout);
OpenMS_Log_warn.insert(os);

START_SECTION((void checkDefaults(const std::string &name, const Param &defaults, const std::string& prefix="") const))
    Param p,d;
    p.setValue("string",std::string("bla"),"string");
    p.setValue("int",5,"int");
    p.setValue("double",47.11,"double");

    p.checkDefaults("Test",d,"");
    TEST_EQUAL(os.str().empty(),false)

    d.setValue("int",5,"int");
    d.setValue("double",47.11,"double");
    os.str("");
  os.clear();
    p.checkDefaults("Test",d,"");
    TEST_EQUAL(os.str().empty(),false)

    p.clear();
    p.setValue("pref:string",std::string("bla"),"pref:string");
    p.setValue("pref:int",5,"pref:int");
    p.setValue("pref:double",47.11,"pref:double");
    os.str("");
  os.clear();
    p.checkDefaults("Test",d,"pref");
    TEST_EQUAL(os.str().empty(),false)

    os.str("");
  os.clear();
    p.checkDefaults("Test2",d,"pref:");
    TEST_EQUAL(os.str().empty(),false)

    //check string restrictions
    vector<std::string> s_rest = {"a","b","c"};
    d.setValue("stringv","bla","desc");
    d.setValidStrings("stringv", s_rest);
    p.clear();
    p.setValue("stringv","a");
    p.checkDefaults("Param_test",d,"");
    p.setValue("stringv","d");
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))

    //check int restrictions
    d.setValue("intv",4,"desc");
    d.setMinInt("intv",-4);
    p.clear();
    p.setValue("intv",-4);
    p.checkDefaults("Param_test",d,"");
    p.setValue("intv",700);
    p.checkDefaults("Param_test",d,"");
    p.setValue("intv",-5);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))

    d.setValue("intv2",4,"desc");
    d.setMaxInt("intv2",4);
    p.clear();
    p.setValue("intv2",4);
    p.checkDefaults("Param_test",d,"");
    p.setValue("intv2",-700);
    p.checkDefaults("Param_test",d,"");
    p.setValue("intv2",5);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))

    //check double restrictions
    d.setValue("doublev",4.0,"desc");
    d.setMinFloat("doublev",-4.0);
    p.clear();
    p.setValue("doublev",-4.0);
    p.checkDefaults("Param_test",d,"");
    p.setValue("doublev",0.0);
    p.checkDefaults("Param_test",d,"");
    p.setValue("doublev",7.0);
    p.checkDefaults("Param_test",d,"");
    p.setValue("doublev",-4.1);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))

    d.setValue("doublev2",4.0,"desc");
    d.setMaxFloat("doublev2",4.0);
    p.clear();
    p.setValue("doublev2",4.0);
    p.checkDefaults("Param_test",d,"");
    p.setValue("doublev2",-700.0);
    p.checkDefaults("Param_test",d,"");
    p.setValue("doublev2",4.1);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))

    //check list restrictions
    vector<std::string> s_rest1 = {"a","b","c"};
    d.setValue("stringlist",std::vector<std::string>{"aaa","abc","cab"},"desc");
    d.setValidStrings("stringlist", s_rest);
    p.clear();
    p.setValue("stringlist",std::vector<std::string>{"a","c"});
    p.checkDefaults("Param_test",d,"");
    p.setValue("stringlist",std::vector<std::string>{"aa","dd","cc"});
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))


    //wrong type
    p.clear();
    p.setValue("doublev",4);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))
    p.clear();
    p.setValue("intv","bla");
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))
    p.clear();
    p.setValue("stringv",4.5);
    TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,""))
END_SECTION

START_SECTION((void update(const Param& old_version, const bool add_unknown = false)))
{
  NOT_TESTABLE // see full implementation below
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


