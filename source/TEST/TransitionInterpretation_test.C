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
// $Maintainer: Andreas Bertsch $
// $Authors:Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MRM/TransitionInterpretation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TransitionInterpretation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionInterpretation* ptr = 0;
START_SECTION(TransitionInterpretation())
{
	ptr = new TransitionInterpretation();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~TransitionInterpretation())
{
	delete ptr;
}
END_SECTION

START_SECTION((TransitionInterpretation(const TransitionInterpretation &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setMZDelta(DoubleReal mz_delta)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getMZDelta() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrimary(bool primary)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getPrimary() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductAdjustment(const String &adjustment)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getProductAdjustment() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductOrdinal(Size ordinal)))
{
  // TODO
}
END_SECTION

START_SECTION((Size getProductOrdinal() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductSeries(const String &series)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getProductSeries() const ))
{
  // TODO
}
END_SECTION

START_SECTION((TransitionInterpretation& operator=(const TransitionInterpretation &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



