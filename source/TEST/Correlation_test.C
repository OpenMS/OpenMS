X-Powered-By: PHP/5.1.4
Content-type: text/html

// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TODO/Correlation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Correlation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Correlation* ptr = 0;
CHECK((Correlation()))
	ptr = new Correlation();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~Correlation()))
	delete ptr;
RESULT

CHECK((Correlation()))
  // ???
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)))
  // ???
RESULT

CHECK((double evaluate(const IndexSet& set, const BaseModel<2>& model)))
  // ???
RESULT

CHECK((static BaseQuality* create()))
  // ???
RESULT

CHECK((static const String getName()))
  // ???
RESULT

CHECK((~Correlation()))
  // ???
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



