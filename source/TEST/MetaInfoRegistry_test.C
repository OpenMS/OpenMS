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

#include <OpenMS/METADATA/MetaInfoRegistry.h>

///////////////////////////

START_TEST(MetaInfoRegistry, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MetaInfoRegistry* test;
CHECK((MetaInfoRegistry()))
	test = new MetaInfoRegistry;
	TEST_NOT_EQUAL(test, 0)
RESULT


CHECK((~MetaInfoRegistry()))
	delete test;
RESULT

MetaInfoRegistry mir;

CHECK((UInt registerName(const String& name, const String& description, const String& unit="") const))
	UInt testname = mir.registerName("testname","this is just a test");
	TEST_EQUAL(1024,testname)
	UInt retention_time = mir.registerName("retention time","this is just another test","sec");
	TEST_EQUAL(1025,retention_time)
	UInt another_testname = mir.registerName("another testname", "i will be set later", "so am i");
	TEST_EQUAL(1026,another_testname);
RESULT

CHECK((void setDescription(UInt index, const String& description) ))
	mir.setDescription(1026, "foo");
	TEST_EQUAL(mir.getDescription(1026), "foo")
RESULT

CHECK((void setDescription(const String& name, const String& description) ))
	mir.setDescription("another testname", "bar");
	TEST_EQUAL(mir.getDescription(1026), "bar")
RESULT

CHECK((void setUnit(UInt index, const String& unit) ))
	mir.setUnit(1026, "foo");
	TEST_EQUAL(mir.getUnit(1026), "foo")
RESULT

CHECK((void setUnit(const String& name, const String& unit) ))
	mir.setUnit("another testname", "bar");
	TEST_EQUAL(mir.getUnit(1026), "bar")
RESULT

CHECK((UInt getIndex(const String& name) const))
	UInt tmp;
	tmp = mir.getIndex ("testname");
	TEST_EQUAL(1024,tmp)
	tmp = mir.getIndex ("retention time");
	TEST_EQUAL(1025,tmp)
	tmp = mir.getIndex ("isotopic_range");
	TEST_EQUAL(1,tmp)
	tmp = mir.getIndex ("cluster_id");
	TEST_EQUAL(2,tmp)
	tmp = mir.getIndex ("unregistered name");
	TEST_EQUAL(1027,tmp)
	tmp = mir.getIndex ("another unregistered name");
	TEST_EQUAL(1028,tmp)
	tmp = mir.getIndex ("unregistered name");
	TEST_EQUAL(1027,tmp)
RESULT

CHECK((String getName(UInt index) const ))
	string tmp;
	tmp = mir.getName (1);
	TEST_EQUAL("isotopic_range",tmp)
	tmp = mir.getName (2);
	TEST_EQUAL("cluster_id",tmp)
	tmp = mir.getName (3);
	TEST_EQUAL("label",tmp)
	tmp = mir.getName (4);
	TEST_EQUAL("icon",tmp)
	tmp = mir.getName (1024);
	TEST_EQUAL("testname",tmp)
	tmp = mir.getName (1025);
	TEST_EQUAL("retention time",tmp)
RESULT

CHECK((String getDescription(UInt index) const ))
	string tmp;
	tmp = mir.getDescription(1024);
	TEST_EQUAL(tmp,string("this is just a test"))
	tmp = mir.getDescription(1025);
	TEST_EQUAL(tmp,string("this is just another test"))
	tmp = mir.getDescription(1);
	TEST_EQUAL(tmp,string("consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak"))
	tmp = mir.getDescription(2);
	TEST_EQUAL(tmp,string("consecutive numbering of isotope clusters in a spectrum"))
RESULT

CHECK((String getDescription(const String& name) const ))
	string tmp;
	tmp = mir.getDescription("testname");
	TEST_EQUAL(tmp,string("this is just a test"))
	tmp = mir.getDescription("retention time");
	TEST_EQUAL(tmp,string("this is just another test"))
	tmp = mir.getDescription("isotopic_range");
	TEST_EQUAL(tmp,string("consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak"))
	tmp = mir.getDescription("cluster_id");
	TEST_EQUAL(tmp,string("consecutive numbering of isotope clusters in a spectrum"))	
RESULT

CHECK((String getUnit(UInt index) const ))
	string tmp;
	tmp = mir.getUnit(1024);
	TEST_EQUAL(tmp,string(""))
	tmp = mir.getUnit(1025);
	TEST_EQUAL(tmp,string("sec"))
	tmp = mir.getUnit(1);
	TEST_EQUAL(tmp,string(""))
	tmp = mir.getUnit(2);
	TEST_EQUAL(tmp,string(""))
RESULT

CHECK((String getUnit(const String& name) const ))
	string tmp;
	tmp = mir.getUnit("testname");
	TEST_EQUAL(tmp,string(""))
	tmp = mir.getUnit("retention time");
	TEST_EQUAL(tmp,string("sec"))
	tmp = mir.getUnit("isotopic_range");
	TEST_EQUAL(tmp,string(""))
	tmp = mir.getUnit("cluster_id");
	TEST_EQUAL(tmp,string(""))	
RESULT

CHECK((MetaInfoRegistry(const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2(mir);
	TEST_EQUAL(1024,mir2.getIndex ("testname"))
	TEST_EQUAL(1025,mir2.getIndex ("retention time"))
	TEST_EQUAL("isotopic_range",mir2.getName (1))
	TEST_EQUAL("testname",mir2.getName (1024))
	TEST_EQUAL("retention time",mir2.getName (1025))
	TEST_EQUAL(mir2.getDescription(1024),string("this is just a test"))
	TEST_EQUAL(mir2.getDescription(1025),string("this is just another test"))
	TEST_EQUAL(mir2.getDescription("testname"),string("this is just a test"))
	TEST_EQUAL(mir2.getDescription("retention time"),string("this is just another test"))	
	TEST_EQUAL(mir2.getUnit(1024),string(""))
	TEST_EQUAL(mir2.getUnit(1025),string("sec"))
	TEST_EQUAL(mir2.getUnit("testname"),string(""))
	TEST_EQUAL(mir2.getUnit("retention time"),string("sec"))	
RESULT

CHECK((MetaInfoRegistry& operator = (const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2;
	mir2 = mir;
	TEST_EQUAL(1024,mir2.getIndex ("testname"))
	TEST_EQUAL(1025,mir2.getIndex ("retention time"))
	TEST_EQUAL("isotopic_range",mir2.getName (1))
	TEST_EQUAL("testname",mir2.getName (1024))
	TEST_EQUAL("retention time",mir2.getName (1025))
	TEST_EQUAL(mir2.getDescription(1024),string("this is just a test"))
	TEST_EQUAL(mir2.getDescription(1025),string("this is just another test"))
	TEST_EQUAL(mir2.getDescription("testname"),string("this is just a test"))
	TEST_EQUAL(mir2.getDescription("retention time"),string("this is just another test"))	
	TEST_EQUAL(mir2.getUnit(1024),string(""))
	TEST_EQUAL(mir2.getUnit(1025),string("sec"))
	TEST_EQUAL(mir2.getUnit("testname"),string(""))
	TEST_EQUAL(mir2.getUnit("retention time"),string("sec"))	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
