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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

// 
// Note: This is only the 1st part of Param_test
//       For the second part of Param_test see Param_test_2
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Param, "$Id$")

//////////////////// Param::ParamEntry /////////////////////////////
////////////////////////////////////////////////////////////////////

Param::ParamEntry* pe_ptr =0;
START_SECTION((Param::ParamEntry()))
	pe_ptr = new Param::ParamEntry();
	TEST_NOT_EQUAL(pe_ptr,0)
END_SECTION

START_SECTION((~Param::ParamEntry()))
	delete pe_ptr;
END_SECTION

START_SECTION((Param::ParamEntry(const String& n, const DataValue& v, const String& d, bool u)))
	Param::ParamEntry pe("n","v","d",StringList::create("advanced"));
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

START_SECTION((bool operator==(const Param::ParamEntry& rhs) const))
	Param::ParamEntry n1("n","d","v",StringList::create("advanced"));
	Param::ParamEntry n2("n","d","v",StringList::create("advanced"));
	
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

Param::ParamNode* pn_ptr =0;
START_SECTION((Param::ParamNode()))
	pn_ptr = new Param::ParamNode();
	TEST_NOT_EQUAL(pn_ptr,0)
END_SECTION

START_SECTION((~Param::ParamNode()))
	delete pn_ptr;
END_SECTION

START_SECTION((Param::ParamNode(const String& n, const String& d)))
	Param::ParamNode n("n","d");
	TEST_EQUAL(n.name,"n")
	TEST_EQUAL(n.description,"d")
	
	n = Param::ParamNode("n1","d1");
	TEST_EQUAL(n.name,"n1")
	TEST_EQUAL(n.description,"d1")
END_SECTION

START_SECTION((bool operator==(const Param::ParamNode& rhs) const))
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

START_SECTION((String suffix(const String& key)))
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


START_SECTION((UInt size() const))
	TEST_EQUAL(pn.size(),4)
	TEST_EQUAL(pn.nodes[0].size(),2)
	TEST_EQUAL(pn.nodes[1].size(),1)
END_SECTION

START_SECTION((EntryIterator findEntry(const String& name)))
	TEST_EQUAL(pn.findEntry("A")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("B")!=pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("C")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("D")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("E")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("F")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("G")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("H")==pn.entries.end(),true)
END_SECTION

START_SECTION((NodeIterator findNode(const String& name)))
	TEST_EQUAL(pn.findNode("A")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("B")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("C")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("D")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("E")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("F")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("G")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("H")==pn.nodes.end(),true)
END_SECTION

START_SECTION((Param::ParamNode* findParentOf(const String& name)))
	TEST_EQUAL(pn.findParentOf("A"),0)
	TEST_EQUAL(pn.findParentOf("B"),&pn)
	TEST_EQUAL(pn.findParentOf("C"),&pn)
	TEST_EQUAL(pn.findParentOf("C:D"),&(pn.nodes[0]))
	TEST_EQUAL(pn.findParentOf("C:E"),&(pn.nodes[0]))
	TEST_EQUAL(pn.findParentOf("F"),0)
	TEST_EQUAL(pn.findParentOf("B:G"),&(pn.nodes[1]))
	TEST_EQUAL(pn.findParentOf("X"),0)
	TEST_EQUAL(pn.findParentOf("H:X"),0)
	TEST_EQUAL(pn.findParentOf("H:C:X"),0)
	TEST_EQUAL(pn.findParentOf("H:C:"),0)
END_SECTION

START_SECTION((Param::ParamEntry* findEntryRecursive(const String& name)))
	TEST_EQUAL(pn.findEntryRecursive("A"),0)
	TEST_EQUAL(pn.findEntryRecursive("B"),&(pn.entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("C"),0)
	TEST_EQUAL(pn.findEntryRecursive("C:D"),&(pn.nodes[0].entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("C:E"),&(pn.nodes[0].entries[1]))
	TEST_EQUAL(pn.findEntryRecursive("F"),0)
	TEST_EQUAL(pn.findEntryRecursive("B:G"),&(pn.nodes[1].entries[0]))
	TEST_EQUAL(pn.findEntryRecursive("X"),0)
	TEST_EQUAL(pn.findEntryRecursive("H:X"),0)
	TEST_EQUAL(pn.findEntryRecursive("H:C:X"),0)
	TEST_EQUAL(pn.findEntryRecursive("H:C:"),0)
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

START_SECTION((void insert(const Param::ParamNode& node, const String& prefix = "")))
	Param::ParamNode node("","");
	node.entries.push_back(Param::ParamEntry("H",5,"",StringList::create("advanced")));
	pn.insert(node,"F");
	TEST_NOT_EQUAL(pn.findEntryRecursive("F:H"),0)

	pn.insert(node,"F:Z");
	TEST_NOT_EQUAL(pn.findEntryRecursive("F:Z:H"),0)

	pn.insert(node,"F:Z:");
	TEST_NOT_EQUAL(pn.findEntryRecursive("F:Z::H"),0)

	pn.insert(node,"FD:ZD:D");
	TEST_NOT_EQUAL(pn.findEntryRecursive("FD:ZD:D:H"),0)
	
	node.name = "W";
	pn.insert(node);
	TEST_NOT_EQUAL(pn.findEntryRecursive("W:H"),0)	

	pn.insert(node,"Q");
	TEST_NOT_EQUAL(pn.findEntryRecursive("QW:H"),0)	
END_SECTION

START_SECTION((void insert(const Param::ParamEntry& entry, const String& prefix = "")))
	Param::ParamEntry entry("H","",5,StringList::create("advanced"));

	pn.insert(entry);
	TEST_NOT_EQUAL(pn.findEntryRecursive("H"),0)
		
	pn.insert(entry,"F");
	TEST_NOT_EQUAL(pn.findEntryRecursive("FH"),0)

	pn.insert(entry,"G:");
	TEST_NOT_EQUAL(pn.findEntryRecursive("G:H"),0)

	pn.insert(entry,"FD:ZD:D");
	TEST_NOT_EQUAL(pn.findEntryRecursive("FD:ZD:DH"),0)
END_SECTION


////////////////// Param::ParamIterator ////////////////////////////
////////////////////////////////////////////////////////////////////


Param::ParamIterator* pi_ptr=0;
START_SECTION((ParamIterator()))
	pi_ptr = new Param::ParamIterator();
	TEST_NOT_EQUAL(pi_ptr,0)
END_SECTION

START_SECTION((~ParamIterator()))
	delete(pi_ptr);
END_SECTION

START_SECTION((ParamIterator(const Param::ParamNode& root)))
	Param::ParamNode node;
	pi_ptr = new Param::ParamIterator(node);
	TEST_NOT_EQUAL(pi_ptr,0)
END_SECTION

START_SECTION((const Param::ParamEntry& operator*()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",StringList::create("advanced")));
	Param::ParamIterator it(node);
	TEST_EQUAL((*it).name,"name")
	TEST_EQUAL((*it).value,"value");
	TEST_EQUAL((*it).description,"description")
	TEST_EQUAL((*it).tags.count("advanced")==1,true)
END_SECTION

START_SECTION((const Param::ParamEntry* operator->()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",StringList::create("advanced")));
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
//  | |
//	| S
//  | |-B=2
//  | |-C=3
//  |-T
//    |-D=4
Param::ParamNode root, r, s, t;
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
root.nodes.push_back(r);
root.nodes.push_back(t);

START_SECTION((ParamIterator& operator++()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",StringList::create("advanced")));
	node.entries.push_back(Param::ParamEntry("name2","value2","description2"));
	node.entries.push_back(Param::ParamEntry("name3","value3","description3",StringList::create("advanced")));

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

START_SECTION((ParamIterator operator++(Int)))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",StringList::create("advanced")));
	node.entries.push_back(Param::ParamEntry("name2","value2","description2"));
	node.entries.push_back(Param::ParamEntry("name3","value3","description3",StringList::create("advanced")));

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

START_SECTION((String getName() const))
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


START_SECTION((bool operator==(const ParamIterator& rhs) const))
	Param::ParamIterator begin(root), begin2(root), end;
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin==begin, true)
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(end==end, true)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	
	++begin2;
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	
	++begin2;
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	
	++begin2;
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)

	++begin;
	TEST_EQUAL(begin==begin2, false)
	TEST_EQUAL(begin==end, true)
	TEST_EQUAL(begin2==end, false)
	
	++begin2;
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==end, true)
	TEST_EQUAL(begin2==end, true)
END_SECTION

START_SECTION((bool operator!=(const ParamIterator& rhs) const))
	Param::ParamIterator begin(root), begin2(root), end;
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==begin, true)
	TEST_EQUAL(begin2==begin2, true)
	TEST_EQUAL(end==end, true)
END_SECTION


START_SECTION((const std::vector< TraceInfo>& getTrace() const))
	
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



Param* d10_ptr = 0;
START_SECTION((Param()))
	d10_ptr = new Param();
	TEST_NOT_EQUAL(d10_ptr, 0)
END_SECTION

START_SECTION((~Param()))
	delete d10_ptr;
END_SECTION

START_SECTION((bool exists(const String& key) const))
	Param p;
	TEST_EQUAL(p.exists(""), false)
	TEST_EQUAL(p.exists("key"), false)	
	TEST_EQUAL(p.exists("key:value"), false)
END_SECTION

START_SECTION((const DataValue& getValue(const String &key) const  ))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getValue("key:value"))
END_SECTION

START_SECTION((const String& getSectionDescription(const String& key) const))
	Param p;
	TEST_EQUAL(p.getSectionDescription(""),"")
	TEST_EQUAL(p.getSectionDescription("key"),"")
	TEST_EQUAL(p.getSectionDescription("key:value"),"")
END_SECTION

START_SECTION((const String& getDescription(const String &key) const  ))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key:value"))
END_SECTION

START_SECTION((const ParamEntry& getEntry(const String &key) const  ))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key:value"))
END_SECTION

START_SECTION((void setValue(const String& key, const String& value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key","value");
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), "value")
	TEST_EQUAL(p.getDescription("key"), "")
	TEST_EQUAL(p.hasTag("key","advanced"), false)

	p.setValue("key","value","description",StringList::create("advanced"));
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

START_SECTION((void setValue(const String& key, Int value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",-5,"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL((Int)p.getValue("key"),-5)
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, UInt value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",5u,"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL((Int)p.getValue("key"),5u)
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, Real value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",11.4f,"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_REAL_SIMILAR(p.getValue("key"), 11.4f)
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, DoubleReal value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",11.5,"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_REAL_SIMILAR(p.getValue("key"), 11.5)
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, StringList value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",StringList::create("a,b,c,d"),"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), StringList::create("a,b,c,d"))
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, IntList value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",IntList::create("1,2,3"),"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), IntList::create("1,2,3"))
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION((void setValue(const String& key, DoubleList value, const String& description="", bool advanced=false)))
	Param p;
	p.setValue("key",DoubleList::create("11.5,3.44"),"description",StringList::create("advanced"));
	TEST_EQUAL(p.exists("key"), true)
	TEST_EQUAL(p.getValue("key"), DoubleList::create("11.5,3.44"))
	TEST_EQUAL(p.getDescription("key"), "description")
	TEST_EQUAL(p.hasTag("key","advanced"), true)
END_SECTION

START_SECTION(StringList getTags(const String& key) const)
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getTags("key"))

	p.setValue("key","value");
	TEST_EQUAL(p.getTags("key").size(),0)
END_SECTION

START_SECTION(void addTag(const String& key, const String& tag))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTag("key","bla"))
	StringList error_list;
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

START_SECTION(bool hasTag(const String& key, const String& tag) const)
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

START_SECTION(void addTags(const String& key, const StringList& tags))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTags("key",StringList()))
	StringList error_list;
	error_list.push_back("a,b");
	TEST_EXCEPTION(Exception::ElementNotFound, p.addTags("key",error_list))

	p.setValue("key","value");
	TEST_EQUAL(p.hasTag("key","advanced"),false)
	TEST_EQUAL(p.hasTag("key","advanced2"),false)
	p.addTags("key",StringList::create("advanced,advanced2"));
	TEST_EQUAL(p.hasTag("key","advanced"),true)
	TEST_EQUAL(p.hasTag("key","advanced2"),true)
END_SECTION

START_SECTION(void clearTags(const String& key))
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

START_SECTION((UInt size() const))
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

START_SECTION((void setSectionDescription(const String &key, const String &description) ))
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
	p.setValue("key",17.4);
	stringstream ss;
	ss << p;
	TEST_EQUAL(ss.str(), "\"key\" -> \"17.4\"\n")
	
	ss.str("");
	p.setValue("key",17.4, "thiskey");
	ss<<p;
	TEST_EQUAL(ss.str(), "\"key\" -> \"17.4\" (thiskey)\n")

	ss.str("");
	p.clear();
	p.setValue("tree:key",17.5);
	ss<<p;
	TEST_EQUAL(ss.str(), "\"tree|key\" -> \"17.5\"\n")
END_SECTION

Param p;
p.setValue("test:float",17.4f,"floatdesc");
p.setValue("test:string","test,test,test","stringdesc");
p.setValue("test:int",17,"intdesc");
p.setValue("test2:float",17.5f);
p.setValue("test2:string","test2");
p.setValue("test2:int",18);
p.setSectionDescription("test","sectiondesc");

START_SECTION((void insert(String prefix, const Param &param)))
	Param p2;
	p2.insert("test3",p);
	
	TEST_REAL_SIMILAR(float(p2.getValue("test3test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test3test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test3test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test3test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test3test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test3test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test3test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test3test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test3test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test3test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test3test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test3test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test3test"),"sectiondesc")
		
	p2.insert("",p);
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)	
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")

	p2.insert("test3:",p);
	
	TEST_REAL_SIMILAR(float(p2.getValue("test3:test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test3:test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test3:test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test3:test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test3:test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test3:test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test3:test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test3:test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test3:test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test3:test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test3:test"),"sectiondesc")
		
	p2.insert("",p);
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)	
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



