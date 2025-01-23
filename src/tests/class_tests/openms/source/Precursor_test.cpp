// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

START_SECTION((double getDriftTimeUnit() const ))
  Precursor tmp;
  TEST_EQUAL(tmp.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
END_SECTION

START_SECTION((void setDriftTimeUnit(double dt)))
  Precursor tmp;
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(tmp.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
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

START_SECTION((double getDriftTimeWindowUpperOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getDriftTimeWindowUpperOffset(), 0);
END_SECTION

START_SECTION((void setDriftTimeWindowUpperOffset(double bound)))
  Precursor tmp;
  tmp.setDriftTimeWindowUpperOffset(22.7);
  TEST_REAL_SIMILAR(tmp.getDriftTimeWindowUpperOffset(), 22.7);
END_SECTION

START_SECTION((double getDriftTimeWindowLowerOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getDriftTimeWindowLowerOffset(), 0);
END_SECTION

START_SECTION((void setDriftTimeWindowLowerOffset(double bound)))
  Precursor tmp;
  tmp.setDriftTimeWindowLowerOffset(22.8);
  TEST_REAL_SIMILAR(tmp.getDriftTimeWindowLowerOffset(), 22.8);
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
{
  Precursor tmp;
  tmp.getActivationMethods().insert(Precursor::CID);
  tmp.setActivationEnergy(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
  tmp.setDriftTime(7.11);
  tmp.setDriftTimeWindowUpperOffset(12.8);
  tmp.setDriftTimeWindowLowerOffset(12.7);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  tmp.setCharge(2);
  tmp.getPossibleChargeStates().resize(2);
  tmp.setMetaValue("label",String("label"));
  
  Precursor tmp2(tmp);
  TEST_EQUAL(tmp2.getActivationMethods().size(),1);
  TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),47.11);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTime(),7.11);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowUpperOffset(), 12.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowLowerOffset(), 12.7);
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(tmp2.getCharge(),2);
  TEST_EQUAL(tmp2.getPossibleChargeStates().size(),2);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
}
END_SECTION

START_SECTION((Precursor(const Precursor&& source)))
{
  Precursor tmp;
  tmp.getActivationMethods().insert(Precursor::CID);
  tmp.getActivationMethods().insert(Precursor::BIRD);
  tmp.setActivationEnergy(40.11);
  tmp.setIsolationWindowUpperOffset(20.7);
  tmp.setIsolationWindowLowerOffset(20.8);
  tmp.setDriftTime(0.11);
  tmp.setDriftTimeWindowUpperOffset(10.8);
  tmp.setDriftTimeWindowLowerOffset(10.7);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  tmp.setCharge(8);
  tmp.getPossibleChargeStates().resize(4);
  tmp.setMetaValue("label",String("label2"));
  TEST_EQUAL(tmp.getActivationMethods().size(),2);

  //copy tmp so we can move one of them
  Precursor orig = tmp;

  Precursor tmp2(std::move(tmp));
  TEST_EQUAL(tmp2, orig);

  TEST_EQUAL(tmp2.getActivationMethods().size(),2);
  TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),40.11);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 20.7);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 20.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTime(),0.11);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowUpperOffset(), 10.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowLowerOffset(), 10.7);
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true);
  TEST_EQUAL(tmp2.getCharge(),8);
  TEST_EQUAL(tmp2.getPossibleChargeStates().size(),4);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label2");
}
END_SECTION

START_SECTION((Precursor& operator= (const Precursor& source)))
{
  Precursor tmp;
  tmp.getActivationMethods().insert(Precursor::CID);
  tmp.setActivationEnergy(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
  tmp.setDriftTime(7.11);
  tmp.setDriftTimeWindowUpperOffset(12.8);
  tmp.setDriftTimeWindowLowerOffset(12.7);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  tmp.setCharge(9);
  tmp.getPossibleChargeStates().resize(5);
  tmp.setMetaValue("label",String("label"));
  
  //normal assignment
  Precursor tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getActivationMethods().size(),1);
  TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),47.11);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTime(),7.11);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowUpperOffset(), 12.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowLowerOffset(), 12.7);
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(tmp2.getCharge(),9);
  TEST_EQUAL(tmp2.getPossibleChargeStates().size(),5);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
    
  //assignment of empty object
  tmp2 = Precursor();
  TEST_EQUAL(tmp2.getActivationMethods().size(),0);
  TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),0.0);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 0.0);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 0.0);
  TEST_REAL_SIMILAR(tmp2.getDriftTime(),-1.0);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowUpperOffset(), 0.0);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowLowerOffset(), 0.0);
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
  TEST_EQUAL(tmp2.getCharge(),0);
  TEST_EQUAL(tmp2.getPossibleChargeStates().size(),0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
}
END_SECTION

START_SECTION((Precursor& operator= (const Precursor&& source)))
{
  Precursor tmp;
  tmp.getActivationMethods().insert(Precursor::CID);
  tmp.getActivationMethods().insert(Precursor::BIRD);
  tmp.setActivationEnergy(40.11);
  tmp.setIsolationWindowUpperOffset(20.7);
  tmp.setIsolationWindowLowerOffset(20.8);
  tmp.setDriftTime(0.11);
  tmp.setDriftTimeWindowUpperOffset(10.8);
  tmp.setDriftTimeWindowLowerOffset(10.7);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  tmp.setCharge(8);
  tmp.getPossibleChargeStates().resize(4);
  tmp.setMetaValue("label",String("label2"));

  //copy tmp so we can move one of them
  Precursor orig = tmp;

  //move assignment
  Precursor tmp2;
  tmp2 = std::move(tmp);
  TEST_EQUAL(tmp2, orig);

  TEST_EQUAL(tmp2.getActivationMethods().size(),2);
  TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),40.11);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 20.7);
  TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 20.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTime(),0.11);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowUpperOffset(), 10.8);
  TEST_REAL_SIMILAR(tmp2.getDriftTimeWindowLowerOffset(), 10.7);
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true);
  TEST_EQUAL(tmp2.getCharge(),8);
  TEST_EQUAL(tmp2.getPossibleChargeStates().size(),4);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label2");
}
END_SECTION

START_SECTION((bool operator== (const Precursor& rhs) const))
	Precursor tmp,tmp2;
	
	TEST_TRUE(tmp == tmp2);
	
	tmp2.setActivationEnergy(47.11);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setDriftTime(5.0);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
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
  tmp.setDriftTimeWindowUpperOffset(22.7);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setDriftTimeWindowLowerOffset(22.8);
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
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setDriftTime(5.0);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
	tmp.getActivationMethods().insert(Precursor::CID);
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setDriftTimeWindowUpperOffset(22.7);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setDriftTimeWindowLowerOffset(22.8);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
  tmp.setCharge(13);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
  tmp.getPossibleChargeStates().resize(5);
	TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_FALSE(tmp == tmp2);
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



