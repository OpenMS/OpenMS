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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChargePair, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChargePair* ptr = 0;
START_SECTION(ChargePair())
{
	ptr = new ChargePair();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ChargePair())
{
	delete ptr;
}
END_SECTION

START_SECTION((ChargePair(const Size &index0, const Size &index1, const Int &charge0, const Int &charge1, const Size &compomer_id, const DoubleReal &mass_diff, const bool active)))
{
	ChargePair cp(34,45, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomerId(), 99);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, 99, 12.34, false);
	ChargePair cp (cp2);
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomerId(), 99);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 0);	
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((ChargePair& operator=(const ChargePair &rhs)))
{
	ChargePair cp2(34,45, 4,5, 99, 12.34, false);
	ChargePair cp = cp2;
	TEST_EQUAL(cp.getElementIndex(0), 34);
	TEST_EQUAL(cp.getElementIndex(1), 45);	
	TEST_EQUAL(cp.getCharge(0), 4);
	TEST_EQUAL(cp.getCharge(1), 5);	
	TEST_EQUAL(cp.getCompomerId(), 99);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 12.34);
	TEST_EQUAL(cp.getEdgeScore(), 0);
	TEST_EQUAL(cp.isActive(), false);	
}
END_SECTION

START_SECTION((virtual ~ChargePair()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((Int getCharge(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...	
}
END_SECTION

START_SECTION((void setCharge(UInt pairID, Int e)))
{
  ChargePair cp;
	cp.setCharge(0,123);
	cp.setCharge(1,321);
	TEST_EQUAL(cp.getCharge(0), 123)
	TEST_EQUAL(cp.getCharge(1), 321)
}
END_SECTION

START_SECTION((Size getElementIndex(UInt pairID) const ))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setElementIndex(UInt pairID, Size e)))
{
  ChargePair cp;
	cp.setElementIndex(0,123);
	cp.setElementIndex(1,321);
	TEST_EQUAL(cp.getElementIndex(0), 123)
	TEST_EQUAL(cp.getElementIndex(1), 321)
}
END_SECTION

START_SECTION((Size getCompomerId() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setCompomerId(Size compomer_id)))
{
  ChargePair cp;
	cp.setCompomerId(123);
	TEST_EQUAL(cp.getCompomerId(), 123)
}
END_SECTION

START_SECTION((DoubleReal getMassDiff() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setMassDiff(DoubleReal mass_diff)))
{
  ChargePair cp;
	cp.setMassDiff(123.432);
	TEST_REAL_SIMILAR(cp.getMassDiff(), 123.432)
}
END_SECTION

START_SECTION((Real getEdgeScore() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setEdgeScore(Real score)))
{
  ChargePair cp;
	cp.setEdgeScore(1123.432f);
	TEST_REAL_SIMILAR(cp.getEdgeScore(), 1123.432)
}
END_SECTION
		

START_SECTION((bool isActive() const))
{
	NOT_TESTABLE //well.. tested below...
}
END_SECTION

START_SECTION((void setActive(const bool active)))
{
  ChargePair cp;
	cp.setActive(true);
	TEST_EQUAL(cp.isActive(), true)
	cp.setActive(false);
	TEST_EQUAL(cp.isActive(), false)
}
END_SECTION

START_SECTION((virtual bool operator==(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, 99, 12.34, false);
	ChargePair cp2(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp1==cp2, false);
	ChargePair cp3(34,15, 4,5, 99, 12.34, true);
	ChargePair cp4(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp3==cp4, false);
	ChargePair cp5(34,15, 4,5, 99, 12.34, false);
	ChargePair cp6(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp5==cp6, true);
	
}
END_SECTION

START_SECTION((virtual bool operator!=(const ChargePair &i) const))
{
	ChargePair cp1(34,45, 4,5, 99, 12.34, false);
	ChargePair cp2(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp1!=cp2, true);
	ChargePair cp3(34,15, 4,5, 99, 12.34, true);
	ChargePair cp4(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp3!=cp4, true);
	ChargePair cp5(34,15, 4,5, 99, 12.34, false);
	ChargePair cp6(34,15, 4,5, 99, 12.34, false);
	TEST_EQUAL(cp5!=cp6, false);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



