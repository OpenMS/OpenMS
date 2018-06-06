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

#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Digestion, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Digestion* dv_ptr = nullptr;
Digestion* dv_nullPointer = nullptr;
START_SECTION((Digestion()))
	dv_ptr = new Digestion;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Digestion()))
	delete dv_ptr;
END_SECTION

//basic accessors
START_SECTION((const String& getEnzyme() const))
	Digestion s;
	TEST_EQUAL(s.getEnzyme(),"")
END_SECTION

//basic accessors
START_SECTION((double getDigestionTime() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getDigestionTime(),0.0)
END_SECTION

//basic accessors
START_SECTION((double getTemperature() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getTemperature(),0.0)
END_SECTION

//basic accessors
START_SECTION((double getPh() const ))
	Digestion s;
	TEST_REAL_SIMILAR(s.getPh(),0.0)
END_SECTION

//basic accessors
START_SECTION((void setEnzyme(const String& enzyme)))
	Digestion s;
	s.setEnzyme("TTEST");
	TEST_EQUAL(s.getEnzyme(),"TTEST")
END_SECTION

//basic accessors
START_SECTION((void setDigestionTime(double digestion_time)))
	Digestion s;
	//set
	s.setDigestionTime(4711.2);
	TEST_REAL_SIMILAR(s.getDigestionTime(),4711.2)
END_SECTION

//basic accessors
START_SECTION((void setTemperature(double temperature)))
	Digestion s;
	s.setTemperature(4711.3);
	TEST_REAL_SIMILAR(s.getTemperature(),4711.3)
END_SECTION

//basic accessors
START_SECTION((void setPh(double ph)))
	Digestion s;
	s.setPh(4711.4);
	TEST_REAL_SIMILAR(s.getPh(),4711.4)
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Digestion s;
	TEST_EQUAL(s.getType(),"Digestion")
END_SECTION

//copy ctr
START_SECTION((Digestion(const Digestion&)))
	Digestion s;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));
	
	//copy
	Digestion s2(s);

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(s2.getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(s2.getTemperature(),4711.3)
	TEST_REAL_SIMILAR(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
END_SECTION

START_SECTION((Digestion& operator=(const Digestion&)))
	Digestion s,s2;
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s2.getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(s2.getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(s2.getTemperature(),4711.3)
	TEST_REAL_SIMILAR(s2.getPh(),4711.4)
	TEST_EQUAL(string(s.getMetaValue("color")),"red")
END_SECTION

START_SECTION((virtual SampleTreatment* clone() const ))
	Digestion s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Digestion* dp;
	
	//set
	s.setEnzyme("TTEST");
	s.setDigestionTime(4711.2);
	s.setTemperature(4711.3);
	s.setPh(4711.4);
	s.setMetaValue("color",String("red"));
	
	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Digestion*>(st);
	
	//get
	TEST_EQUAL(dp->getEnzyme(),"TTEST")
	TEST_REAL_SIMILAR(dp->getDigestionTime(),4711.2)
	TEST_REAL_SIMILAR(dp->getTemperature(),4711.3)
	TEST_REAL_SIMILAR(dp->getPh(),4711.4)
	TEST_EQUAL(string(dp->getMetaValue("color")),"red")
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Digestion empty,edit;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.setEnzyme("TTEST");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setDigestionTime(4711.2);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setTemperature(4711.3);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);			

	edit.setPh(4711.4);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);	
	
	Modification m;
	TEST_EQUAL(m==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
