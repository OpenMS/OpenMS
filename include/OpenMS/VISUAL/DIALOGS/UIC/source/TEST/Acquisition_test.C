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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/Acquisition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Acquisition, "$Id: Acquisition_test.C 6135 2009-10-19 16:05:59Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Acquisition* ptr = 0;
START_SECTION(Acquisition())
	ptr = new Acquisition();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~Acquisition())
	delete ptr;
END_SECTION

START_SECTION(const String& getIdentifier() const)
  Acquisition tmp;
  TEST_EQUAL(tmp.getIdentifier(), "");
END_SECTION

START_SECTION(void setIdentifier(const String& identifier))
	Acquisition tmp;
	tmp.setIdentifier("5");
  TEST_EQUAL(tmp.getIdentifier(), "5");
END_SECTION

START_SECTION(Acquisition(const Acquisition& source))
	Acquisition tmp;
	tmp.setIdentifier("5");
	tmp.setMetaValue("label",String("label"));
	Acquisition tmp2(tmp);
	TEST_EQUAL(tmp2.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
END_SECTION

START_SECTION(Acquisition& operator= (const Acquisition& source))
	Acquisition tmp,tmp2,tmp3;
	// assignment of a modified object
	tmp2.setIdentifier("5");
	tmp2.setMetaValue("label",String("label"));
	tmp = tmp2;
	TEST_EQUAL(tmp.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp.getMetaValue("label")), String("label"));
	
	// assignment of a default-constructed object
	tmp = tmp3;
	TEST_EQUAL(tmp.getIdentifier(), "");
	TEST_EQUAL(tmp.isMetaEmpty(), true);	
END_SECTION

START_SECTION(bool operator== (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setIdentifier("5");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION(bool operator!= (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setIdentifier("5");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



