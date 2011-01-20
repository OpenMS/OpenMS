// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/Precursor.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Precursor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Precursor* ptr = 0;
START_SECTION((Precursor()))
	ptr = new Precursor();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~Precursor()))
	delete ptr;
END_SECTION

START_SECTION((DoubleReal getActivationEnergy() const))
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationEnergy(),0);
END_SECTION

START_SECTION((void setActivationEnergy(DoubleReal activation_energy)))
  Precursor tmp;
  tmp.setActivationEnergy(47.11);
  TEST_REAL_SIMILAR(tmp.getActivationEnergy(),47.11);
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

START_SECTION((DoubleReal getIsolationWindowUpperOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowUpperOffset(DoubleReal bound)))
  Precursor tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 22.7);
END_SECTION

START_SECTION((DoubleReal getIsolationWindowLowerOffset() const))
  Precursor tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowLowerOffset(DoubleReal bound)))
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
END_SECTION

START_SECTION((Precursor& operator= (const Precursor& source)))
	Precursor tmp;
	tmp.setActivationEnergy(47.11);
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
		
	//assignment of empty object
	tmp2 = Precursor();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getActivationMethods().size(),0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getActivationEnergy(),0.0);
END_SECTION

START_SECTION((bool operator== (const Precursor& rhs) const))
	Precursor tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setActivationEnergy(47.11);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



