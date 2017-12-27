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
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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
    }

    Test(const Test& source):
    SampleTreatment(source)
    {
    }

    ~Test() override
    {
    }

    Test& operator = (const Test& source)
    {
      if (&source != this)
      {
        SampleTreatment::operator=(source);
      }
      return *this;
    }

    SampleTreatment* clone() const override
    {
      return new Test(*this);
    }

    bool operator== (const SampleTreatment& rhs) const override
    {
      if (type_ != rhs.getType()) 
      {
        return false;
      }

      const Test* tmp = dynamic_cast<const Test*>(&rhs);
      return SampleTreatment::operator==(*tmp);
    }
};


START_TEST(SampleTreatment, "$Id$")

/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.001)

Test* dv_ptr = nullptr;
Test* dv_nullPointer = nullptr;
START_SECTION((SampleTreatment(const String& type)))
	dv_ptr = new Test;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
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
