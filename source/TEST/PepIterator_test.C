// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/CHEMISTRY/PepIterator.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIterator.h>
#include <OpenMS/FORMAT/FastaIteratorIntern.h>
#include <OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>
#include <OpenMS/CHEMISTRY/TrypticIterator.h>

using namespace OpenMS;
using namespace std;

START_TEST(PepIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(PepIterator())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~PepIterator()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((PepIterator(const PepIterator &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual FASTAEntry operator *()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual PepIterator& operator++()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual PepIterator* operator++(int)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setFastaFile(const String &f)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual String getFastaFile()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setSpectrum(const std::vector< DoubleReal > &s)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual const std::vector<DoubleReal>& getSpectrum()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setTolerance(DoubleReal t)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual DoubleReal getTolerance()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool begin()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool isAtEnd()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void registerChildren()))
{
	PepIterator* p1 = Factory<PepIterator>::create("EdwardsLippertIterator");
	TEST_NOT_EQUAL(dynamic_cast<EdwardsLippertIterator*>(p1),0)
	p1 = Factory<PepIterator>::create("EdwardsLippertIteratorTryptic");
	TEST_NOT_EQUAL(dynamic_cast<EdwardsLippertIteratorTryptic*>(p1),0)
	p1 = Factory<PepIterator>::create("TrypticIterator");
	TEST_NOT_EQUAL(dynamic_cast<TrypticIterator*>(p1),0)
	p1 = Factory<PepIterator>::create("FastaIterator");
	TEST_NOT_EQUAL(dynamic_cast<FastaIterator*>(p1),0)
	p1 = Factory<PepIterator>::create("FastaIteratorIntern");
	TEST_NOT_EQUAL(dynamic_cast<FastaIteratorIntern*>(p1),0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



