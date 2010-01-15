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
#include <OpenMS/METADATA/Product.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Product, "$Id: Product_test.C 6135 2009-10-19 16:05:59Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Product* ptr = 0;
START_SECTION((Product()))
	ptr = new Product();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~Product()))
	delete ptr;
END_SECTION

START_SECTION((DoubleReal getMZ() const ))
  Product tmp;
  TEST_EQUAL(tmp.getMZ(),0);
END_SECTION

START_SECTION((void setMZ(DoubleReal mz)))
  Product tmp;
  tmp.setMZ(47.11);
  TEST_REAL_SIMILAR(tmp.getMZ(),47.11);
END_SECTION

START_SECTION((DoubleReal getIsolationWindowUpperOffset() const ))
  Product tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowUpperOffset(DoubleReal bound)))
  Product tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 22.7);
END_SECTION

START_SECTION((DoubleReal getIsolationWindowLowerOffset() const ))
  Product tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowLowerOffset(DoubleReal bound)))
  Product tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 22.8);
END_SECTION

START_SECTION((Product(const Product& source)))
	Product tmp;
	tmp.setMZ(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	Product tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getMZ(),47.11);
END_SECTION

START_SECTION((Product& operator= (const Product& source)))
	Product tmp;
	tmp.setMZ(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	//normal assignment
	Product tmp2;
	tmp2 = tmp;
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getMZ(),47.11);
		
	//assignment of empty object
	tmp2 = Product();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getMZ(),0.0);
END_SECTION

START_SECTION((bool operator== (const Product& rhs) const))
	Product tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setMZ(47.11);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION((bool operator!= (const Product& rhs) const))
	Product tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setMZ(47.11);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



