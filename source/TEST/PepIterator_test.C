// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl,Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/PepIterator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PepIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepIterator* ptr = 0;
CHECK(PepIterator())
{
	NOT_TESTABLE
}
RESULT

CHECK((virtual ~PepIterator()))
{
  NOT_TESTABLE
}
RESULT

CHECK((PepIterator(const PepIterator &source)))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual FASTAEntry operator *()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual PepIterator& operator++()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual PepIterator* operator++(int)=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual void setFastaFile(const String &f)=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual String getFastaFile()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual void setSpectrum(const std::vector< float > &s)=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual const std::vector<float>& getSpectrum()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual void setTolerance(float t)=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual float getTolerance()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual bool begin()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((virtual bool isAtEnd()=0))
{
  NOT_TESTABLE
}
RESULT

CHECK((void registerChildren()))
{
	PepIterator* p1 = Factory<PepIterator>::create("EdwardsLippertIterator");
	TEST_STRING_EQUAL(p1->getName(), "EdwardsLippertIterator")
	p1 = Factory<PepIterator>::create("EdwardsLippertIteratorTryptic");
	TEST_STRING_EQUAL(p1->getName(), "EdwardsLippertIteratorTryptic")
	p1 = Factory<PepIterator>::create("TrypticIterator");
	TEST_STRING_EQUAL(p1->getName(), "TrypticIterator")
	p1 = Factory<PepIterator>::create("FastaIterator");
	TEST_STRING_EQUAL(p1->getName(), "FastaIterator")
	p1 = Factory<PepIterator>::create("FastaIteratorIntern");
	TEST_STRING_EQUAL(p1->getName(), "FastaIteratorIntern")
}
RESULT

CHECK(String getProductName())
	TEST_STRING_EQUAL(PepIterator::getProductName(), "PepIterator")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



