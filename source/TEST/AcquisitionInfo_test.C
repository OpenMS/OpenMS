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
#include <OpenMS/METADATA/AcquisitionInfo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AcquisitionInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AcquisitionInfo* ptr = 0;
START_SECTION(AcquisitionInfo())
	ptr = new AcquisitionInfo();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~AcquisitionInfo())
	delete ptr;
END_SECTION

START_SECTION(const String& getMethodOfCombination() const)
  AcquisitionInfo tmp;
  TEST_EQUAL(tmp.getMethodOfCombination(),"");
END_SECTION

START_SECTION(void setMethodOfCombination(const String& method_of_combination))
  AcquisitionInfo tmp;
  tmp.setMethodOfCombination("TEST");
  TEST_EQUAL(tmp.getMethodOfCombination(),"TEST");
END_SECTION

START_SECTION(AcquisitionInfo(const AcquisitionInfo& source))
	AcquisitionInfo tmp;
	Acquisition a;
	a.setIdentifier("4711");
	tmp.push_back(a);
	tmp.setMethodOfCombination("Combo");
	tmp.setMetaValue("bla",4.0);

	AcquisitionInfo tmp2(tmp);
	TEST_EQUAL(tmp2.size(), 1);
	TEST_EQUAL(tmp2[0].getIdentifier(), "4711");  
	TEST_EQUAL(tmp2.getMethodOfCombination(), "Combo");  
	TEST_REAL_SIMILAR((DoubleReal)(tmp2.getMetaValue("bla")), 4.0)
END_SECTION

START_SECTION(AcquisitionInfo& operator= (const AcquisitionInfo& source))
	AcquisitionInfo tmp;
	Acquisition a;
	a.setIdentifier("4711");
	tmp.push_back(a);
	tmp.setMethodOfCombination("Combo");
	tmp.setMetaValue("bla",4.0);

	//normal assignment
	AcquisitionInfo tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.size(), 1);
	TEST_EQUAL(tmp2[0].getIdentifier(), "4711");  
	TEST_EQUAL(tmp2.getMethodOfCombination(), "Combo");
	TEST_REAL_SIMILAR((DoubleReal)(tmp2.getMetaValue("bla")), 4.0)
		
	//assignment of a empty value
	tmp2 = AcquisitionInfo();
	TEST_EQUAL(tmp2.size(), 0);
	TEST_EQUAL(tmp2.getMethodOfCombination(), "");
	TEST_EQUAL(tmp2.metaValueExists("bla"), false)
END_SECTION

START_SECTION(bool operator== (const AcquisitionInfo& rhs) const)
  AcquisitionInfo empty,edit;
	TEST_EQUAL(empty==edit,true);
	
	Acquisition a;
	edit.push_back(a);
	TEST_EQUAL(empty==edit,false);

	edit.setMetaValue("bla",4.0);
	TEST_EQUAL(empty==edit,false);

	edit = empty;
	edit.setMethodOfCombination("Combo");
	TEST_EQUAL(empty==edit,false);
END_SECTION

START_SECTION(bool operator!= (const AcquisitionInfo& rhs) const)
  AcquisitionInfo empty,edit;
	TEST_EQUAL(empty!=edit,false);
	
	Acquisition a;
	edit.push_back(a);
	TEST_EQUAL(empty!=edit,true);

	edit.setMetaValue("bla",4.0);
	TEST_EQUAL(empty!=edit,true);

	edit = empty;
	edit.setMethodOfCombination("Combo");
	TEST_EQUAL(empty!=edit,true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



