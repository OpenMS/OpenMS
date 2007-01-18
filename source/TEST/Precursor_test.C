// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/Precursor.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Precursor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Precursor* ptr = 0;
CHECK(Precursor())
	ptr = new Precursor();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Precursor())
	delete ptr;
RESULT

CHECK(float getActivationEnergy() const)
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationEnergy(),0);
RESULT

CHECK(void setActivationEnergy(float activation_energy))
  Precursor tmp;
  tmp.setActivationEnergy(47.11);
  TEST_REAL_EQUAL(tmp.getActivationEnergy(),47.11);
RESULT

CHECK(ActivationMethod getActivationMethod() const)
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationMethod(),Precursor::ACTMETHNULL);
RESULT

CHECK(void setActivationMethod(ActivationMethod activation_method))
  Precursor tmp;
  tmp.setActivationMethod(Precursor::CID);
  TEST_EQUAL(tmp.getActivationMethod(),Precursor::CID);
RESULT

CHECK(EnergyUnits getActivationEnergyUnit() const)
  Precursor tmp;
  TEST_EQUAL(tmp.getActivationEnergyUnit(),Precursor::UNITSNULL);
RESULT

CHECK(void setActivationEnergyUnit(EnergyUnits activation_energy_unit))
  Precursor tmp;
  tmp.setActivationEnergyUnit(Precursor::EV);
  TEST_EQUAL(tmp.getActivationEnergyUnit(),Precursor::EV);
RESULT

CHECK(float getWindowSize() const)
  Precursor tmp;
  TEST_REAL_EQUAL(tmp.getWindowSize(), 0);
RESULT

CHECK(void setWindowSize(float size))
  Precursor tmp;
  tmp.setWindowSize(22.7);
  TEST_REAL_EQUAL(tmp.getWindowSize(), 22.7);
RESULT

CHECK(Precursor(const Precursor& source))
	Precursor tmp;
	tmp.setActivationEnergy(47.11);
	tmp.setActivationMethod(Precursor::CID);
	tmp.setActivationEnergyUnit(Precursor::EV);
	tmp.setWindowSize(22.7);
	tmp.setMetaValue("label",String("label"));
	
	Precursor tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getActivationMethod(),Precursor::CID);
	TEST_EQUAL(tmp2.getActivationEnergyUnit(),Precursor::EV);
	TEST_REAL_EQUAL(tmp2.getWindowSize(), 22.7);
	TEST_REAL_EQUAL(tmp2.getActivationEnergy(),47.11);
RESULT

CHECK(Precursor& operator= (const Precursor& source))
	Precursor tmp;
	tmp.setActivationEnergy(47.11);
	tmp.setActivationMethod(Precursor::CID);
	tmp.setActivationEnergyUnit(Precursor::EV);
	tmp.setWindowSize(22.7);
	tmp.setMetaValue("label",String("label"));
	
	//normal assignment
	Precursor tmp2;
	tmp2 = tmp;
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getActivationMethod(),Precursor::CID);
	TEST_EQUAL(tmp2.getActivationEnergyUnit(),Precursor::EV);
	TEST_REAL_EQUAL(tmp2.getWindowSize(), 22.7);
	TEST_REAL_EQUAL(tmp2.getActivationEnergy(),47.11);
		
	//assignment of empty object
	tmp2 = Precursor();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getActivationMethod(),Precursor::ACTMETHNULL);
	TEST_EQUAL(tmp2.getActivationEnergyUnit(),Precursor::UNITSNULL);
	TEST_REAL_EQUAL(tmp2.getWindowSize(), 0.0);
	TEST_REAL_EQUAL(tmp2.getActivationEnergy(),0.0);
RESULT

CHECK(bool operator== (const Precursor& rhs) const)
	Precursor tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setActivationEnergy(47.11);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setActivationMethod(Precursor::CID);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setActivationEnergyUnit(Precursor::EV);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setWindowSize(22.7);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
RESULT

CHECK(bool operator!= (const Precursor& rhs) const)
	Precursor tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setActivationEnergy(47.11);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setActivationMethod(Precursor::CID);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setActivationEnergyUnit(Precursor::EV);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setWindowSize(22.7);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



