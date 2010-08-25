// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

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
START_SECTION(([Param::ParamEntry] ParamEntry()))
	pe_ptr = new Param::ParamEntry();
	TEST_NOT_EQUAL(pe_ptr,0)
END_SECTION

START_SECTION(([Param::ParamEntry] ~ParamEntry()))
	delete pe_ptr;
END_SECTION

START_SECTION(([Param::ParamEntry] ParamEntry(const String& n, const DataValue& v, const String& d, bool u)))
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

START_SECTION(([Param::ParamEntry] bool isValid(String& message) const))

	Param p;
	String m;
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


	vector<String> strings;
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
START_SECTION(([Param::ParamNode] ParamNode()))
	pn_ptr = new Param::ParamNode();
	TEST_NOT_EQUAL(pn_ptr,0)
END_SECTION

START_SECTION(([Param::ParamNode] ~ParamNode()))
	delete pn_ptr;
END_SECTION

START_SECTION(([Param::ParamNode] ParamNode(const String& n, const String& d)))
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

START_SECTION(([Param::ParamNode] String suffix(const String& key)))
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

START_SECTION(([Param::ParamNode] EntryIterator findEntry(const String& name)))
	TEST_EQUAL(pn.findEntry("A")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("B")!=pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("C")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("D")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("E")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("F")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("G")==pn.entries.end(),true)
	TEST_EQUAL(pn.findEntry("H")==pn.entries.end(),true)
END_SECTION

START_SECTION(([Param::ParamNode] NodeIterator findNode(const String& name)))
	TEST_EQUAL(pn.findNode("A")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("B")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("C")!=pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("D")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("E")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("F")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("G")==pn.nodes.end(),true)
	TEST_EQUAL(pn.findNode("H")==pn.nodes.end(),true)
END_SECTION

START_SECTION(([Param::ParamNode] Param::ParamNode* findParentOf(const String& name)))
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

START_SECTION(([Param::ParamNode] ParamEntry* findEntryRecursive(const String& name)))
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

START_SECTION(([Param::ParamNode] void insert(const ParamNode& node, const String& prefix = "")))
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

START_SECTION(([Param::ParamNode] void insert(const ParamEntry& entry, const String& prefix = "")))
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
START_SECTION(([Param::ParamIterator] ParamIterator()))
	pi_ptr = new Param::ParamIterator();
	TEST_NOT_EQUAL(pi_ptr,0)
END_SECTION

START_SECTION(([Param::ParamIterator] ~ParamIterator()))
	delete(pi_ptr);
END_SECTION

START_SECTION(([Param::ParamIterator] ParamIterator(const Param::ParamNode& root)))
	Param::ParamNode node;
	pi_ptr = new Param::ParamIterator(node);
	TEST_NOT_EQUAL(pi_ptr,0)
END_SECTION

START_SECTION(([Param::ParamIterator] const Param::ParamEntry& operator*()))
	Param::ParamNode node;
	node.entries.push_back(Param::ParamEntry("name","value","description",StringList::create("advanced")));
	Param::ParamIterator it(node);
	TEST_EQUAL((*it).name,"name")
	TEST_EQUAL((*it).value,"value");
	TEST_EQUAL((*it).description,"description")
	TEST_EQUAL((*it).tags.count("advanced")==1,true)
END_SECTION

START_SECTION(([Param::ParamIterator] const Param::ParamEntry* operator->()))
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

START_SECTION(([Param::ParamIterator] ParamIterator& operator++()))
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

START_SECTION(([Param::ParamIterator] ParamIterator operator++(Int)))
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

START_SECTION(([Param::ParamIterator] String getName() const))
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

START_SECTION(([Param::ParamIterator] bool operator!=(const ParamIterator& rhs) const))
	Param::ParamIterator begin(root), begin2(root), end;
	TEST_EQUAL(begin==end, false)
	TEST_EQUAL(begin2==end, false)
	TEST_EQUAL(begin==begin2, true)
	TEST_EQUAL(begin==begin, true)
	TEST_EQUAL(begin2==begin2, true)
	TEST_EQUAL(end==end, true)
END_SECTION


START_SECTION(([Param::ParamIterator] const std::vector< TraceInfo>& getTrace() const))
	
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

START_SECTION((const DataValue& getValue(const String &key) const))
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

START_SECTION((const String& getDescription(const String &key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getDescription("key:value"))
END_SECTION

START_SECTION((const ParamEntry& getEntry(const String &key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry(""))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key"))
	TEST_EXCEPTION(Exception::ElementNotFound, p.getEntry("key:value"))
END_SECTION

START_SECTION((void setValue(const String &key, const DataValue& value, const String &description="", const StringList &tags=StringList())))
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

START_SECTION((StringList getTags(const String& key) const))
	Param p;
	TEST_EXCEPTION(Exception::ElementNotFound, p.getTags("key"))

	p.setValue("key","value");
	TEST_EQUAL(p.getTags("key").size(),0)
END_SECTION

START_SECTION((void addTag(const String& key, const String& tag)))
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

START_SECTION((bool hasTag(const String& key, const String& tag) const))
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

START_SECTION((void addTags(const String& key, const StringList& tags)))
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

START_SECTION((void clearTags(const String& key)))
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

START_SECTION((void setSectionDescription(const String &key, const String &description)))
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

START_SECTION((void insert(const String& prefix, const Param &param)))
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
	TEST_STRING_EQUAL(p2.getDescription("prefixn2:d"), String::EMPTY)
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
	TEST_STRING_EQUAL(p2.getDescription("n2:d"), String::EMPTY)
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
	TEST_STRING_EQUAL(p2.getDescription("n3:n2:d"), String::EMPTY)
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

Param p;
p.setValue("test:float",17.4f,"floatdesc");
p.setValue("test:string","test,test,test","stringdesc");
p.setValue("test:int",17,"intdesc");
p.setValue("test2:float",17.5f);
p.setValue("test2:string","test2");
p.setValue("test2:int",18);
p.setSectionDescription("test","sectiondesc");
p.addTags("test:float", StringList::create("a,b,c"));

START_SECTION((Param(const Param& rhs)))
	Param p2(p);
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
	TEST_EQUAL(p2.getTags("test:float").size(), 3)
	TEST_EQUAL(p2.getTags("test:float") == StringList::create("a,b,c"), true)
END_SECTION

START_SECTION((Param& operator = (const Param& rhs)))
	Param p2;
	p2=p;
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
	TEST_EQUAL(p2.getTags("test:float").size(), 3)
	TEST_EQUAL(p2.getTags("test:float") == StringList::create("a,b,c"), true)
END_SECTION

START_SECTION((Param copy(const String &prefix, bool remove_prefix=false) const))
	Param p2;

	p2 = p.copy("notthere:");
	TEST_EQUAL((p2==Param()),true)

	p2 = p.copy("test:");
	
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EXCEPTION(Exception::ElementNotFound, p2.getValue("test2:float"))

	p2 = p.copy("test:",true);
	TEST_REAL_SIMILAR(float(p2.getValue("float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("float"), "floatdesc")
	TEST_EQUAL(p2.getValue("string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("string"), "stringdesc")

	p2 = p.copy("test");
	TEST_REAL_SIMILAR(float(p2.getValue("test:float")), 17.4)
	TEST_STRING_EQUAL(p2.getDescription("test:float"), "floatdesc")
	TEST_EQUAL(p2.getValue("test:string"), "test,test,test")
	TEST_STRING_EQUAL(p2.getDescription("test:string"), "stringdesc")
	TEST_EQUAL(Int(p2.getValue("test:int")), 17)
	TEST_STRING_EQUAL(p2.getDescription("test:int"), "intdesc")
	TEST_REAL_SIMILAR(float(p2.getValue("test2:float")), 17.5)
	TEST_STRING_EQUAL(p2.getDescription("test2:float"), String::EMPTY)
	TEST_EQUAL(p2.getValue("test2:string"), "test2")
	TEST_STRING_EQUAL(p2.getDescription("test2:string"), String::EMPTY)
	TEST_EQUAL(Int(p2.getValue("test2:int")), 18)
	TEST_STRING_EQUAL(p2.getDescription("test2:int"), String::EMPTY)
	TEST_EQUAL(p2.getSectionDescription("test"),"sectiondesc")
END_SECTION

START_SECTION((void remove(const String& key)))
	Param p2(p);
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

END_SECTION

START_SECTION((void removeAll(const String& prefix)))
	Param p2(p);
	
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
	Param p2(p);
	TEST_EQUAL(p==p2, true)
	p2.setValue("test:float",17.5f);
	TEST_EQUAL(p==p2, false)
	p2 = p;
	p2.setValue("test:float3",17.4f);
	TEST_EQUAL(p==p2, false)
	p2 = p;
	p2.removeAll("test:float");
	TEST_EQUAL(p==p2, false)
	
	//it should be independent of entry order
	Param p3,p4;
	p3.setValue("1",1);
	p3.setValue("2",2);
	p4.setValue("2",2);
	p4.setValue("1",1);
	TEST_EQUAL(p3==p4, true)

	//it should be independent of node order
	Param p5,p6;
	p5.setValue("1:1",1);
	p5.setValue("2:1",1);
	p6.setValue("2:1",1);
	p6.setValue("1:1",1);
	TEST_EQUAL(p5==p6, true)

END_SECTION

START_SECTION((void load(const String& filename)))
	Param p2;
	TEST_EXCEPTION(Exception::FileNotFound, p2.load("FileDoesNotExist.xml"))	
END_SECTION

START_SECTION((void store(const String& filename) const))
	Param p2(p);
	p2.setValue("test:a:a1", 47.1,"a1desc\"<>\nnewline");
	p2.setValue("test:b:b1", 47.1);
	p2.setSectionDescription("test:b","bdesc\"<>\nnewline");
	p2.setValue("test2:a:a1", 47.1);
	p2.setValue("test2:b:b1", 47.1,"",StringList::create("advanced"));
	p2.setSectionDescription("test2:a","adesc");
	
	//exception
	Param p300;
	TEST_EXCEPTION(Exception::UnableToCreateFile, p300.store("/does/not/exist/FileDoesNotExist.xml"))	
	
	String filename;
	NEW_TMP_FILE(filename);
	p2.store(filename);
	Param p3;
	p3.load(filename);
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
	TEST_EQUAL(p2.isValid(filename),true)
	
	//advanced
	NEW_TMP_FILE(filename);
	Param p7;
	p7.setValue("true",5,"",StringList::create("advanced"));
	p7.setValue("false",5,"");
	
	p7.store(filename);
	TEST_EQUAL(p7.isValid(filename),true)
	Param p8;
	p8.load(filename);
	
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

	vector<String> strings;
	p5.setValue("string","bli");
	strings.push_back("bla");
	strings.push_back("bluff");
	p5.setValue("string_2","bla");
	p5.setValidStrings("string_2",strings);
	
		//list restrictions
	vector<String> strings2;
	strings2.push_back("xml");
	strings2.push_back("txt");
	p5.setValue("stringlist2",StringList::create("a.txt,b.xml,c.pdf"));
	p5.setValue("stringlist",StringList::create("aa.C,bb.h,c.doxygen"));
	p5.setValidStrings("stringlist2",strings2);
	
	p5.setValue("intlist",IntList::create("2,5,10"));
	p5.setValue("intlist2",IntList::create("2,5,10"));
	p5.setValue("intlist3",IntList::create("2,5,10"));
	p5.setValue("intlist4",IntList::create("2,5,10"));
	p5.setMinInt("intlist2",1);
	p5.setMaxInt("intlist3",11);
	p5.setMinInt("intlist4",0);
	p5.setMaxInt("intlist4",15);
	
	p5.setValue("doublelist",DoubleList::create("1.2,3.33,4.44"));
	p5.setValue("doublelist2",DoubleList::create("1.2,3.33,4.44"));
	p5.setValue("doublelist3",DoubleList::create("1.2,3.33,4.44"));
	p5.setValue("doublelist4",DoubleList::create("1.2,3.33,4.44"));			
	
	p5.setMinFloat("doublelist2",1.1);
	p5.setMaxFloat("doublelist3",4.45);
	p5.setMinFloat("doublelist4",0.1);
	p5.setMaxFloat("doublelist4",5.8);
	
	
	p5.store(filename);
	TEST_EQUAL(p5.isValid(filename),true)
	Param p6;
	p6.load(filename);
	
	
	TEST_EQUAL(p6.getEntry("int").min_int, -numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("int").max_int, numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("int_min").min_int, 4)
	TEST_EQUAL(p6.getEntry("int_min").max_int, numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("int_max").min_int, -numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("int_max").max_int, 6)
	TEST_EQUAL(p6.getEntry("int_min_max").min_int, 0)
	TEST_EQUAL(p6.getEntry("int_min_max").max_int, 10)

	TEST_REAL_SIMILAR(p6.getEntry("float").min_float, -numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("float").max_float, numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("float_min").min_float, 4.1)
	TEST_REAL_SIMILAR(p6.getEntry("float_min").max_float, numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("float_max").min_float, -numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("float_max").max_float, 6.1)
	TEST_REAL_SIMILAR(p6.getEntry("float_min_max").min_float, 0.1)
	TEST_REAL_SIMILAR(p6.getEntry("float_min_max").max_float, 10.1)

	TEST_EQUAL(p6.getEntry("string").valid_strings.size(),0)
	TEST_EQUAL(p6.getEntry("string_2").valid_strings.size(),2)
	TEST_EQUAL(p6.getEntry("string_2").valid_strings[0],"bla")
	TEST_EQUAL(p6.getEntry("string_2").valid_strings[1],"bluff")


	
	TEST_EQUAL(p6.getEntry("stringlist").valid_strings.size(),0)
	TEST_EQUAL(p6.getEntry("stringlist2").valid_strings.size(),2)
	TEST_EQUAL(p6.getEntry("stringlist2").valid_strings[0],"xml")
	TEST_EQUAL(p6.getEntry("stringlist2").valid_strings[1],"txt")
	
	TEST_EQUAL(p6.getEntry("intlist").min_int, -numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("intlist").max_int, numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("intlist2").min_int, 1)
	TEST_EQUAL(p6.getEntry("intlist2").max_int, numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("intlist3").min_int, -numeric_limits<Int>::max())
	TEST_EQUAL(p6.getEntry("intlist3").max_int, 11)
	TEST_EQUAL(p6.getEntry("intlist4").min_int, 0)
	TEST_EQUAL(p6.getEntry("intlist4").max_int, 15)
	
	TEST_REAL_SIMILAR(p6.getEntry("doublelist").min_float, -numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("doublelist").max_float, numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("doublelist2").min_float, 1.1)
	TEST_REAL_SIMILAR(p6.getEntry("doublelist2").max_float, numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("doublelist3").min_float, -numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(p6.getEntry("doublelist3").max_float, 4.45)
	TEST_REAL_SIMILAR(p6.getEntry("doublelist4").min_float, 0.1)
	TEST_REAL_SIMILAR(p6.getEntry("doublelist4").max_float, 5.8)
	//Test if an empty Param written to a file validates against the schema
	NEW_TMP_FILE(filename);
	Param p4;
	p4.store(filename);
	TEST_EQUAL(p4.isValid(filename),true)
END_SECTION


START_SECTION((void setDefaults(const Param& defaults, const String& prefix="", bool showMessage=false)))
	Param defaults;
	defaults.setValue("float",1.0f,"float");	
	defaults.setValue("float2",2.0f,"float2");
	defaults.setValue("string","default string1","string");
	defaults.setValue("string2","default string2","string2");
	defaults.setValue("PATH:onlyfordescription",45.2);
	
	defaults.setValue("stringlist",StringList::create("a,b,c"),"stringlist");
	defaults.setValue("stringlist2",StringList::create("d,e,f"),"stringlist2");
	defaults.setValue("intlist",IntList::create("1,2,3"),"intlist");
	defaults.setValue("intlist2",IntList::create("11,22,33"),"intlist2");
	defaults.setValue("doublelist",DoubleList::create("1.2,2.3"),"doublelist");
	defaults.setValue("doublelist2",DoubleList::create("11.22,22.33"),"doublelist2");
	defaults.setSectionDescription("PATH","PATHdesc");
	Param p2;
	p2.setValue("PATH:float",-1.0f,"PATH:float");
	p2.setValue("PATH:string","some string","PATH:string");
	p2.setValue("float",-2.0f,"float");
	p2.setValue("string","other string","string");
	
	p2.setValue("PATH:stringlist",StringList::create("d,a,v,i,d"),"PATH:stringlist");
	p2.setValue("stringlist",StringList::create("r,o,c,k,s"),"stringlist");
	p2.setValue("PATH:intlist2",IntList::create("14,9"),"PATH:intlist2");
	p2.setValue("intlist", IntList::create("16,9"),"intlist");
	p2.setValue("PATH:doublelist2",DoubleList::create("6.66,6.16"),"PATH:doublelist2");
	p2.setValue("doublelist",DoubleList::create("1.2,5.55"),"doublelist");
	
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
	
	TEST_EQUAL(p2.getValue("stringlist"), StringList::create("r,o,c,k,s"))
	TEST_EQUAL(p2.getValue("intlist"),IntList::create("16,9"))
	TEST_EQUAL(p2.getValue("doublelist"),DoubleList::create("1.2,5.55"))
	TEST_EQUAL(p2.getValue("stringlist2"),StringList::create("d,e,f"))
	TEST_EQUAL(p2.getValue("intlist2"),IntList::create("11,22,33"))
	TEST_EQUAL(p2.getValue("doublelist2"),DoubleList::create("11.22,22.33"))
	
	
	
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
	
	TEST_EQUAL(p2.getValue("PATH:stringlist"),StringList::create("d,a,v,i,d"))
	TEST_EQUAL(p2.getValue("PATH:intlist"), IntList::create("1,2,3"))
	TEST_EQUAL(p2.getValue("PATH:doublelist"),DoubleList::create("1.2,2.3"))
	
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

const char* command_line4[10]; // "executable -a -1.0 -b bv -c cv rv1 rv2"
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

START_SECTION((void parseCommandLine(const int argc, const char **argv, const String& prefix="")))
	Param p2,p3;
	p2.parseCommandLine(9,command_line,"test4");
	p3.setValue("test4:-a","av");
	p3.setValue("test4:-b","bv");
	p3.setValue("test4:-c","cv");
	p3.setValue("test4:misc",StringList::create("rv1,rv2"));
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
	p300.setValue("test4:misc",StringList::create("rv1,rv2,-1.0"));
	TEST_EQUAL(p200==p300,true)

END_SECTION


START_SECTION((void parseCommandLine(const int argc, const char **argv, const Map< String, String > &options_with_one_argument, const Map< String, String > &options_without_argument, const Map< String, String > &options_with_multiple_argument, const String &misc="misc", const String &unknown="unknown")))

	Map<String,String> with_one,without,with_multiple;
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
	p3.setValue("misc_",StringList::create("rv1,rv2,-1.0"));
	TEST_EQUAL(p2==p3,true)
	
	Param p4,p5;
	p4.parseCommandLine(9,command_line,with_one,without,with_multiple,"misc_","unknown_");
	p5.setValue("a","av");
	p5.setValue("b","bv");
	p5.setValue("c","cv");
	p5.setValue("misc_",StringList::create("rv1,rv2"));
	TEST_EQUAL(p4==p5,true)

	with_one.clear();
	with_one["-a"]="a";
	without["-b"]="b";
	
	Param p40,p50;
	p40.parseCommandLine(9,command_line,with_one,without,with_multiple,"misc__","unknown__");
	p50.setValue("a","av");
	p50.setValue("b","true");
	p50.setValue("misc__",StringList::create("bv,cv,rv1,rv2"));
	p50.setValue("unknown__",StringList::create("-c"));
	TEST_EQUAL(p40==p50,true)
	TEST_EQUAL(p40,p50)
	//"executable -a av -b -c cv"	
	Param p400,p500;
	p400.parseCommandLine(6,command_line2,with_one,without,with_multiple,"misc__","unknown__");
	p500.setValue("a","av");
	p500.setValue("b","true");
	p500.setValue("misc__",StringList::create("cv"));
	p500.setValue("unknown__",StringList::create("-c"));
	TEST_EQUAL(p400==p500,true)

	//"executable -a -b -c cv rv1"
	Param p4000,p5000;
	p4000.parseCommandLine(6,command_line3,with_one,without,with_multiple,"misc__","unknown__");
	p5000.setValue("a","");
	p5000.setValue("b","true");
	p5000.setValue("misc__",StringList::create("cv,rv1"));
	p5000.setValue("unknown__",StringList::create("-c"));
	TEST_EQUAL(p4000==p5000,true)
		
	//"mult -d 1.333 2.23 3 -e 4 -f -g"
	const char* a1 ="mult";
	const char* a2 ="-d";
	const char* a3 ="1.333";
	const char* a4 ="2.23";
	const char* a5 ="3";
	const char* a6 ="-e";
	const char* a7 ="4";
	const char* a8 ="-f";
	const char* a9 ="-g";

	const char* command_line_mult[9];
	command_line_mult[0] = a1;
	command_line_mult[1] = a2;
	command_line_mult[2] = a3;
	command_line_mult[3] = a4;
	command_line_mult[4] = a5;
	command_line_mult[5] = a6;
	command_line_mult[6] = a7;
	command_line_mult[7] = a8;
	command_line_mult[8] = a9;
	Param p6,p7;
	p6.parseCommandLine(9,command_line_mult,with_one,without,with_multiple,"misc__","unkown__");
	p7.setValue("d",StringList::create("1.333,2.23,3"));
	p7.setValue("e",StringList::create("4"));
	p7.setValue("f",StringList());
	p7.setValue("g",StringList());
	TEST_EQUAL(p6,p7);
	TEST_EQUAL(StringList(),StringList());
	const char* command_line_mult1[4];
	command_line_mult1[0] = a1;
	command_line_mult1[1] = a2;
	command_line_mult1[2] = a3;
	command_line_mult1[3] = a4;
	Param p8,p9;
	p9.parseCommandLine(4,command_line_mult,with_one,without,with_multiple,"misc__","unkown__");
	p8.setValue("d",StringList::create("1.333,2.23"));
	TEST_EQUAL(p9,p8);
	
END_SECTION

START_SECTION((void setValidStrings(const String &key, const std::vector< String > &strings)))
  vector<String> strings;
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

START_SECTION((void setMinInt(const String &key, Int min)))
  Param d;
  d.setValue("ok",4);
  d.setValue("dummy",5.5);
  
  d.setMinInt("ok",4);
  TEST_EQUAL(d.getEntry("ok").min_int,4);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMinInt("dummy",4))
END_SECTION

START_SECTION((void setMaxInt(const String &key, Int max)))
  Param d;
  d.setValue("ok",4);
  d.setValue("dummy",5.5);
  
  d.setMaxInt("ok",4);
  TEST_EQUAL(d.getEntry("ok").max_int,4);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMaxInt("dummy",4))
END_SECTION

START_SECTION((void setMinFloat(const String &key, DoubleReal min)))
  Param d;
  d.setValue("ok",4.5);
  d.setValue("dummy",4);
  
  d.setMinFloat("ok",4.0);
  TEST_REAL_SIMILAR(d.getEntry("ok").min_float,4.0);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMinFloat("dummy",4.5))
END_SECTION

START_SECTION((void setMaxFloat(const String &key, DoubleReal max)))
  Param d;
  d.setValue("ok",4.5);
  d.setValue("dummy",4);
  
  d.setMaxFloat("ok",4.0);
  TEST_REAL_SIMILAR(d.getEntry("ok").max_float,4.0);
  TEST_EXCEPTION(Exception::ElementNotFound, d.setMaxFloat("dummy",4.5))
END_SECTION


START_SECTION((void checkDefaults(const String &name, const Param &defaults, const String& prefix="", std::ostream &os=std::cout) const))
	//warnings for unknown parameters
	ostringstream os;
	Param p,d;
	p.setValue("string",String("bla"),"string");
	p.setValue("int",5,"int");
	p.setValue("double",47.11,"double");
		
	p.checkDefaults("Test",d,"",os);
	TEST_EQUAL(os.str()=="",false);
	
	d.setValue("int",5,"int");
	d.setValue("double",47.11,"double");
	os.str("");
	p.checkDefaults("Test",d,"",os);
	TEST_EQUAL(os.str()=="",false);
	
	p.clear();
	p.setValue("pref:string",String("bla"),"pref:string");
	p.setValue("pref:int",5,"pref:int");
	p.setValue("pref:double",47.11,"pref:double");
	os.str("");
	p.checkDefaults("Test",d,"pref",os);
	TEST_EQUAL(os.str()=="",false);

	os.str("");
	p.checkDefaults("Test",d,"pref:",os);
	TEST_EQUAL(os.str()=="",false);
	
	//check string restrictions
	vector<String> s_rest;
	s_rest.push_back("a");
	s_rest.push_back("b");
	s_rest.push_back("c");
	d.setValue("stringv","bla","desc");
	d.setValidStrings("stringv", s_rest);
	p.clear();
	p.setValue("stringv","a");
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("stringv","d");
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os))

	//check int restrictions
	d.setValue("intv",4,"desc");
	d.setMinInt("intv",-4);
	p.clear();
	p.setValue("intv",-4);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("intv",700);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("intv",-5);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));

	d.setValue("intv2",4,"desc");
	d.setMaxInt("intv2",4);
	p.clear();
	p.setValue("intv2",4);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("intv2",-700);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("intv2",5);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));

	//check double restrictions
	d.setValue("doublev",4.0,"desc");
	d.setMinFloat("doublev",-4.0);
	p.clear();
	p.setValue("doublev",-4.0);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("doublev",0.0);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("doublev",7.0);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("doublev",-4.1);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));
	
	d.setValue("doublev2",4.0,"desc");
	d.setMaxFloat("doublev2",4.0);
	p.clear();
	p.setValue("doublev2",4.0);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("doublev2",-700.0);
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("doublev2",4.1);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));
	
	//check list restrictions
	vector<String> s_rest1;
	s_rest1.push_back("a");
	s_rest1.push_back("b");
	s_rest1.push_back("c");
	d.setValue("stringlist",StringList::create("aaa,abc,cab"),"desc");
	d.setValidStrings("stringlist", s_rest);
	p.clear();
	p.setValue("stringlist",StringList::create("a,c"));
	p.checkDefaults("Param_test",d,"",os);
	p.setValue("stringlist",StringList::create("aa,dd,cc"));
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os))
	
	
	//wrong type
	p.clear();
	p.setValue("doublev",4);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));
	p.clear();
	p.setValue("intv","bla");
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));
	p.clear();
	p.setValue("stringv",4.5);
	TEST_EXCEPTION(Exception::InvalidParameter,p.checkDefaults("Param_test",d,"",os));
END_SECTION

START_SECTION((void update(const Param& old_version, const bool report_new_params=false)))
	Param common;
	common.setValue("float",1.0f,"float");	
	common.setValue("float2",2.0f,"float2");
	common.setValue("string","default string1","string");
	common.setValue("string2","default string2","string2");
	common.setValue("PATH:onlyfordescription",45.2);
	
	common.setValue("stringlist",StringList::create("a,b,c"),"stringlist");
	common.setValue("stringlist2",StringList::create("d,e,f"),"stringlist2");
	common.setValue("intlist",IntList::create("1,2,3"),"intlist");

  // copy and alter
  Param old = common;
  //old.setValue("recently_removed_float",1.1f,"float");  // should not make it into new param
  old.setValue("old_type","a string","string");
  old.setValue("some:version","1.2","old version");
  old.setValue("some:type","unlabeled","type");
	old.setValue("stringlist2",StringList::create("d,e,f,altered"),"stringlist2"); // change some values, we expect them to show up after update()
	old.setValue("intlist",IntList::create("3"),"intlist");
  
  Param defaults = common;
  defaults.setValue("old_type",3,"old_type has evolved from string to int"); // as type has changed, this value should be kept
  defaults.setValue("some:version","1.9","new version"); // this value should be kept (due to its reserved name)
  defaults.setValue("some:type","information","type");   // this value should be kept (due to its reserved name)
  defaults.setValue("new_value",3,"new param not present in old");
  
  Param expected = defaults;
	expected.setValue("stringlist2",StringList::create("d,e,f,altered"),"stringlist2"); // change some values, we expect them to show up after update()
	expected.setValue("intlist",IntList::create("3"),"intlist");
  
  // update()
  defaults.update(old);
  
  TEST_EQUAL(defaults,expected)
  
END_SECTION


START_SECTION((ParamIterator begin() const))
	NOT_TESTABLE
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

START_SECTION([EXTRA] loading and storing of lists)
	Param p;
	p.setValue("stringlist", StringList::create("a,bb,ccc"));
	p.setValue("intlist", IntList::create("1,22,333"));
	p.setValue("item", String("bla"));
	p.setValue("stringlist2", StringList::create(""));
	p.setValue("intlist2", IntList::create(""));
	p.setValue("item1", 7);
	p.setValue("intlist3", IntList::create("1"));	
	p.setValue("stringlist3", StringList::create("1"));
	p.setValue("item3", 7.6);
	p.setValue("doublelist",DoubleList::create("1.22,2.33,4.55"));
	p.setValue("doublelist2",DoubleList::create(""));
	p.setValue("doublelist3",DoubleList::create("1.4"));
	//store
	String filename;
	NEW_TMP_FILE(filename);
	p.store(filename);
	//load
	Param p2;
	p2.load(filename);
	
	TEST_EQUAL(p2.size(),12);
	
	TEST_EQUAL(p2.getValue("stringlist").valueType(), DataValue::STRING_LIST)
	StringList list = p2.getValue("stringlist");
	TEST_EQUAL(list.size(),3)
	TEST_EQUAL(list[0],"a")
	TEST_EQUAL(list[1],"bb")
	TEST_EQUAL(list[2],"ccc")
	
	TEST_EQUAL(p2.getValue("stringlist2").valueType(), DataValue::STRING_LIST)
	list = p2.getValue("stringlist2");
	TEST_EQUAL(list.size(),0)
	
	TEST_EQUAL(p2.getValue("stringlist").valueType(), DataValue::STRING_LIST)
	list = p2.getValue("stringlist3");
	TEST_EQUAL(list.size(),1)
	TEST_EQUAL(list[0],"1")
	
	TEST_EQUAL(p2.getValue("intlist").valueType(), DataValue::INT_LIST)
	IntList intlist = p2.getValue("intlist");
	TEST_EQUAL(intlist.size(),3);
	TEST_EQUAL(intlist[0], 1)
	TEST_EQUAL(intlist[1], 22)
	TEST_EQUAL(intlist[2], 333)
	
	TEST_EQUAL(p2.getValue("intlist2").valueType(),DataValue::INT_LIST)
	intlist = p2.getValue("intlist2");
	TEST_EQUAL(intlist.size(),0)
	
	TEST_EQUAL(p2.getValue("intlist3").valueType(),DataValue::INT_LIST)
	intlist = p2.getValue("intlist3");
	TEST_EQUAL(intlist.size(),1)
	TEST_EQUAL(intlist[0],1)
	
		TEST_EQUAL(p2.getValue("doublelist").valueType(), DataValue::DOUBLE_LIST)
	DoubleList doublelist = p2.getValue("doublelist");
	TEST_EQUAL(doublelist.size(),3);
	TEST_EQUAL(doublelist[0], 1.22)
	TEST_EQUAL(doublelist[1], 2.33)
	TEST_EQUAL(doublelist[2], 4.55)
	
	TEST_EQUAL(p2.getValue("doublelist2").valueType(),DataValue::DOUBLE_LIST)
	doublelist = p2.getValue("doublelist2");
	TEST_EQUAL(doublelist.size(),0)
	
	TEST_EQUAL(p2.getValue("doublelist3").valueType(),DataValue::DOUBLE_LIST)
	doublelist = p2.getValue("doublelist3");
	TEST_EQUAL(doublelist.size(),1)
	TEST_EQUAL(doublelist[0],1.4)
	
END_SECTION


START_SECTION(([EXTRA] Escapingi of characters))
	Param p;
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
	p.store(filename);

	Param p2;
	p2.load(filename);

	TEST_STRING_EQUAL(p2.getDescription("string"), "string")

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


