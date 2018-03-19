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

#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/Modification.h>
#include <sstream>

///////////////////////////

START_TEST(Tagging, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Tagging* dv_ptr = nullptr;
Tagging* dv_nullPointer = nullptr;
START_SECTION((Tagging()))
	dv_ptr = new Tagging;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Tagging()))
	delete dv_ptr;
END_SECTION

START_SECTION((const IsotopeVariant& getVariant() const))
	Tagging s;
	TEST_EQUAL(s.getVariant(),Tagging::LIGHT)
END_SECTION

START_SECTION((double getMassShift() const ))
	Tagging s;
	TEST_REAL_SIMILAR(s.getMassShift(),0.0)
END_SECTION

START_SECTION((void setMassShift(double mass_shift)))
	Tagging s;
	s.setMassShift(4711.2);
	TEST_REAL_SIMILAR(s.getMassShift(),4711.2)
END_SECTION

START_SECTION((void setVariant(const IsotopeVariant& variant)))
	Tagging s;
	s.setVariant(Tagging::HEAVY);
	TEST_EQUAL(s.getVariant(),Tagging::HEAVY)
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Tagging s;
	TEST_EQUAL(s.getType(),"Tagging")
END_SECTION

//copy ctr
START_SECTION((Tagging(const Tagging&)))
	Tagging s;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//copy
	Tagging s2(s);

	//get
	TEST_REAL_SIMILAR(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(s2.getMass(),23.4)
END_SECTION

//assignment operator
START_SECTION((Tagging& operator=(const Tagging&)))
	Tagging s,s2;
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//assign
	s2 = s;

	//get
	TEST_REAL_SIMILAR(s2.getMassShift(),4711.2)
	TEST_EQUAL(s2.getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(s2.getMass(),23.4)
END_SECTION

//clone
START_SECTION((virtual SampleTreatment* clone() const ))
	Tagging s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Tagging* dp;
	
	//set
	s.setMassShift(4711.2);
	s.setVariant(Tagging::LIGHT);
	s.setMass(23.4);
	
	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Tagging*>(st);
	
	//get
	TEST_REAL_SIMILAR(dp->getMassShift(),4711.2)
	TEST_EQUAL(dp->getVariant(),Tagging::LIGHT)
	TEST_REAL_SIMILAR(dp->getMass(),23.4)
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Tagging empty,edit;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.setMassShift(4711.2);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setVariant(Tagging::HEAVY);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setMass(23.4);
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
