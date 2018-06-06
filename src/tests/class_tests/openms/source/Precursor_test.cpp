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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/Precursor.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Precursor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Precursor* ptr = nullptr;
Precursor* nullPointer = nullptr;
START_SECTION((Precursor()))
	ptr = new Precursor();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~Precursor()))
	delete ptr;
END_SECTION

START_SECTION((double getActivationEnergy() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationEnergy(),0);
END_SECTION

START_SECTION((void setActivationEnergy(double activation_energy)))
  Precursor tmp;
  tmp.setActivationEnergy(47.11);
  TEST_REAL_SIMILAR(tmp.getActivationEnergy(),47.11);
END_SECTION

START_SECTION((double getDriftTime() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getDriftTime(),-1);
END_SECTION

START_SECTION((void setDriftTime(double dt)))
  Precursor tmp;
  tmp.setDriftTime(47.11);
  TEST_REAL_SIMILAR(tmp.getDriftTime(),47.11);
END_SECTION

START_SECTION((const set<ActivationMethod>& getActivationMethods() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationMethods().size(),0);
END_SECTION


START_SECTION((set<ActivationMethod>& getActivationMethods()))
  Precursor tmp;
	tmp.getActivationMethods().insert(Precursor::CID);
  TEST_EQUAL(tmp.getActivationMethods().size(),1);
END_SECTION

START_SECTION((void setActivationMethods(const set<ActivationMethod>& activation_methods)))
  Precursor tmp;
	set<Precursor::ActivationMethod> methods;
	methods.insert(Precursor::CID);
	tmp.setActivationMethods(methods);
  TEST_EQUAL(tmp.getActivationMethods().size(),1);
END_SECTION

START_SECTION((double getIsolationWindowUpperOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowUpperOffset(double bound)))
  Precursor tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 22.7);
END_SECTION

START_SECTION((double getIsolationWindowLowerOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowLowerOffset(double bound)))
  Precursor tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 22.8);
END_SECTION

START_SECTION((Int getCharge() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getCharge(), 0);
END_SECTION

START_SECTION((void setCharge(Int charge)))
  Precursor tmp;
  tmp.setCharge(2);
  TEST_EQUAL(tmp.getCharge(), 2);
END_SECTION

START_SECTION((const std::vector<Int>& getPossibleChargeStates() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getPossibleChargeStates().size(), 0);
END_SECTION

START_SECTION((std::vector<Int>& getPossibleChargeStates()))
  Precursor tmp;
  tmp.getPossibleChargeStates().resize(1);
  TEST_EQUAL(tmp.getPossibleChargeStates().size(), 1);
END_SECTION

START_SECTION((void setPossibleChargeStates(const std::vector<Int>& possible_charge_states)))
  Precursor tmp;
  vector<Int> states(1);
  tmp.setPossibleChargeStates(states);
  TEST_EQUAL(tmp.getPossibleChargeStates().size(), 1);
END_SECTION

START_SECTION((Precursor(const Precursor& source)))
	Precursor tmp;
	tmp.setActivationEnergy(47.11);
	tmp.setDriftTime(7.11);
	tmp.getActivationMethods().insert(Precursor::CID);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	Precursor tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getActivationMethods().size(),1);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),47.11);
	TEST_REAL_SIMILAR(tmp2.getDriftTime(),7.11);
END_SECTION

START_SECTION((Precursor& operator= (const Precursor& source)))
	Precursor tmp;
	tmp.setActivationEnergy(47.11);
	tmp.setDriftTime(7.11);
	tmp.getActivationMethods().insert(Precursor::CID);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	//normal assignment
	Precursor tmp2;
	tmp2 = tmp;
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getActivationMethods().size(),1);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),47.11);
	TEST_REAL_SIMILAR(tmp2.getDriftTime(),7.11);
		
	//assignment of empty object
	tmp2 = Precursor();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getActivationMethods().size(),0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),0.0);
	TEST_REAL_SIMILAR(tmp2.getDriftTime(),-1.0);
END_SECTION

START_SECTION((bool operator== (const Precursor& rhs) const))
	Precursor tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setActivationEnergy(47.11);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setDriftTime(5.0);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.getActivationMethods().insert(Precursor::CID);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setCharge(13);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.getPossibleChargeStates().resize(5);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION((bool operator!= (const Precursor& rhs) const))
	Precursor tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setActivationEnergy(47.11);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setDriftTime(5.0);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.getActivationMethods().insert(Precursor::CID);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
  tmp.setCharge(13);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
  tmp.getPossibleChargeStates().resize(5);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

START_SECTION(double getUnchargedMass() const)
  Precursor tmp;
  tmp.setMZ(123);
  tmp.setCharge(13);
  TEST_REAL_SIMILAR(tmp.getUnchargedMass(), 1585.90540593198);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



