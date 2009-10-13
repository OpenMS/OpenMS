// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Tagging.h>
#include <sstream>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

///////////////////////////

class Test: public SampleTreatment
{
  public:
    Test():
    SampleTreatment("Test")
    {

    };

    Test(const Test& source):
    SampleTreatment(source)
    {

    };

    virtual ~Test()
    {

    };

    Test& operator = (const Test& source)
    {
      if (&source != this)
      {
        SampleTreatment::operator=(source);
      }
      return *this;
    };


    virtual SampleTreatment* clone() const
    {
      return new Test(*this);
    };

    virtual bool operator== (const SampleTreatment& rhs) const
    {
      if (type_!=rhs.getType()) return false;

      const Test* tmp = dynamic_cast<const Test*>(&rhs);
      return
        SampleTreatment::operator==(*tmp)
        ;
    };
};


START_TEST(SampleTreatment, "$Id$")

/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.001)

Test* dv_ptr = 0;
START_SECTION((SampleTreatment(const String& type)))
	dv_ptr = new Test;
	TEST_NOT_EQUAL(dv_ptr, 0)
END_SECTION

START_SECTION((~SampleTreatment()))
	delete dv_ptr;
END_SECTION

START_SECTION(const String& getType() const)
	Test s;
	TEST_EQUAL(s.getType(),"Test")
END_SECTION

START_SECTION((const String& getComment() const))
	Test s;
	TEST_EQUAL(s.getComment(),"")
END_SECTION

START_SECTION(void setComment(const String& comment))
	Test s;
	s.setComment("blubb");
	TEST_EQUAL(s.getComment(),"blubb");
END_SECTION

START_SECTION([EXTRA] MetaInfo)
	Test s;
	//empty
	TEST_EQUAL(s.isMetaEmpty(),true)

	s.setMetaValue("origin",String("cow"));
	s.setMetaValue("size",1.0);
	TEST_EQUAL(s.isMetaEmpty(),false)
	TEST_EQUAL(String(s.getMetaValue("origin")),"cow")
	TEST_REAL_SIMILAR(double(s.getMetaValue("size")),1.0)
END_SECTION

START_SECTION((SampleTreatment(const SampleTreatment&)))
	Test s;
	//set
	s.setComment("TTEST");
	s.setMetaValue("origin",String("horse"));
	//copy
	Test s2(s);
	//get
	TEST_EQUAL(s2.getComment(),"TTEST")
	TEST_EQUAL(s.getMetaValue("origin"),"horse")
END_SECTION

START_SECTION((SampleTreatment& operator=(const SampleTreatment&)))
	Test s,s2;
	//set
	s.setComment("TTEST");
	s.setMetaValue("origin",String("horse"));
	//assign
	s2 = s;
	//get
	TEST_EQUAL(s2.getComment(),"TTEST")
	TEST_EQUAL(s.getMetaValue("origin"),"horse")
END_SECTION

START_SECTION((virtual SampleTreatment* clone() const=0))
	Test s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Test* dp;

	//set
	s.setComment("TTEST");
	s.setMetaValue("origin",String("horse"));

	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Test*>(st);

	//get
	TEST_EQUAL(dp->getComment(),"TTEST")
	TEST_EQUAL(dp->getMetaValue("origin"),"horse")
END_SECTION

START_SECTION((bool operator== (const SampleTreatment& rhs) const))
	Test edit,empty;

	edit.setComment("bla");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	Tagging t;
	TEST_EQUAL(t==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
