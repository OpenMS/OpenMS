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
#include <OpenMS/METADATA/Acquisition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Acquisition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Acquisition* ptr = 0;
CHECK(Acquisition())
	ptr = new Acquisition();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Acquisition())
	delete ptr;
RESULT

CHECK(const SignedInt getNumber() const)
  Acquisition tmp;
  TEST_EQUAL(tmp.getNumber(), -1);
RESULT

CHECK(void setNumber(const SignedInt number))
	Acquisition tmp;
	tmp.setNumber(5);
  TEST_EQUAL(tmp.getNumber(), 5);
RESULT

CHECK(Acquisition(const Acquisition& source))
	Acquisition tmp;
	tmp.setNumber(5);
	tmp.setMetaValue("label",String("label"));
	Acquisition tmp2(tmp);
	TEST_EQUAL(tmp2.getNumber(), 5);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
RESULT

CHECK(Acquisition& operator= (const Acquisition& source))
	Acquisition tmp,tmp2,tmp3;
	// assignment of a modified object
	tmp2.setNumber(5);
	tmp2.setMetaValue("label",String("label"));
	tmp = tmp2;
	TEST_EQUAL(tmp.getNumber(), 5);
	TEST_EQUAL((String)(tmp.getMetaValue("label")), String("label"));
	
	// assignment of a default-constructed object
	tmp = tmp3;
	TEST_EQUAL(tmp.getNumber(), -1);
	TEST_EQUAL(tmp.isMetaEmpty(), true);	
RESULT

CHECK(bool operator== (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setNumber(5);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
RESULT

CHECK(bool operator!= (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setNumber(5);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



