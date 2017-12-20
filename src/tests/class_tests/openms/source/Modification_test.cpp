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

#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <sstream>

///////////////////////////

START_TEST(Modification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Modification* dv_ptr = nullptr;
Modification* dv_nullPointer = nullptr;
START_SECTION((Modification()))
	dv_ptr = new Modification;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Modification()))
	delete dv_ptr;
END_SECTION

START_SECTION((const String& getReagentName() const))
	Modification s;
	TEST_EQUAL(s.getReagentName(),"")
END_SECTION

START_SECTION((double getMass() const ))
	Modification s;
	TEST_REAL_SIMILAR(s.getMass(),0.0)
END_SECTION

START_SECTION((const SpecificityType& getSpecificityType() const))
	Modification s;
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
END_SECTION

START_SECTION((const String& getAffectedAminoAcids() const))
	Modification s;
	TEST_EQUAL(s.getAffectedAminoAcids(),"")
END_SECTION

START_SECTION((void setReagentName(const String& reagent_name)))
	Modification s;
	s.setReagentName("TTEST");
	TEST_EQUAL(s.getReagentName(),"TTEST")
END_SECTION

START_SECTION((void setMass(double mass)))
	Modification s;
	s.setMass(11.9);
	TEST_REAL_SIMILAR(s.getMass(),11.9)
END_SECTION

START_SECTION((void setSpecificityType(const SpecificityType& specificity_type)))
	Modification s;
	s.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(s.getSpecificityType(),Modification::CTERM)
END_SECTION

START_SECTION((void setAffectedAminoAcids(const String& affected_amino_acids)))
	Modification s;
	s.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
END_SECTION

//getType
START_SECTION([EXTRA] getType)
	Modification s;
	TEST_EQUAL(s.getType(),"Modification")
END_SECTION

//copy ctr
START_SECTION((Modification(const Modification&)))
	Modification s;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",String("red"));

	//copy
	Modification s2(s);

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_SIMILAR(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
	
END_SECTION

//assignment operator
START_SECTION((Modification& operator=(const Modification&)))
	Modification s,s2;
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",String("red"));

	//assign
	s2 = s;

	//get
	TEST_EQUAL(s.getReagentName(),"TTEST")
	TEST_REAL_SIMILAR(s.getMass(),11.9)
	TEST_EQUAL(s.getSpecificityType(),Modification::AA)
	TEST_EQUAL(s.getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(s.getMetaValue("color")),"red")
END_SECTION

//clone
START_SECTION((virtual SampleTreatment* clone() const ))
	Modification s;
	SampleTreatment* st1;
	SampleTreatment* st;
	Modification* dp;
	
	//set
	s.setReagentName("TTEST");
	s.setMass(11.9);
	s.setSpecificityType(Modification::AA);
	s.setAffectedAminoAcids("ABCDE");
	s.setMetaValue("color",String("red"));

	//assign
	st1 = &s;
	st = st1->clone();
	dp = dynamic_cast<Modification*>(st);

	//get
	TEST_EQUAL(dp->getReagentName(),"TTEST")
	TEST_REAL_SIMILAR(dp->getMass(),11.9)
	TEST_EQUAL(dp->getSpecificityType(),Modification::AA)
	TEST_EQUAL(dp->getAffectedAminoAcids(),"ABCDE")
	TEST_EQUAL(String(dp->getMetaValue("color")),"red")
END_SECTION

START_SECTION((virtual bool operator==(const SampleTreatment &rhs) const ))
	Modification empty,edit;
	
	TEST_EQUAL(edit==empty, true);

	edit.setMass(11.9);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);

	edit.setSpecificityType(Modification::CTERM);
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);		

	edit.setAffectedAminoAcids("ABCDE");
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);			

	edit.setMetaValue("color",String("red"));
	TEST_EQUAL(edit==empty, false);
	edit = empty;
	TEST_EQUAL(edit==empty, true);	
	
	Tagging m;
	TEST_EQUAL(m==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
