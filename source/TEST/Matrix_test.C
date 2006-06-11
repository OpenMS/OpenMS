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
// $Id: Matrix_test.C,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

// Includes go here....
#include <OpenMS/DATASTRUCTURES/Matrix.h>

///////////////////////////

START_TEST(Matrix, "$Id: Matrix_test.C,v 1.4 2006/03/28 12:53:13 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

STATUS("\n\nDISCLAIMER: This ist not really a full test, but it shows that the thing compiles and you may inspect the STATUS() output manually ...\n\n\n")

Matrix<int>* ptr = 0;
CHECK(Matrix<int>())
	ptr = new Matrix<int>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Matrix<int>())
	delete ptr;
RESULT

CHECK(Matrix())
  Matrix<int> mi1;
  STATUS("mi1:\n"<< mi1)
RESULT

Matrix<int> mi;

CHECK((void resize(size_type i, size_type j, value_type value = value_type())))
	mi.resize(2,2,3);
  STATUS("mi1:\n"<< mi)
	mi.resize(2,3,7);
  STATUS("mi1:\n"<< mi)
RESULT

CHECK(Matrix(const Matrix & source))
  Matrix<int> mi2(mi);
  STATUS("mi2:\n"<< mi2)
RESULT

CHECK(Matrix& operator = (const Matrix & rhs))
	Matrix<int> mi3;
  STATUS("mi3:\n"<<mi3)
	mi3=mi;
  STATUS("mi3:\n"<<mi3)
RESULT

CHECK((const_reference getValue(size_type const i, size_type const j) const))
 	mi(1,1)=17;
  STATUS("mi:\n"<<mi)
RESULT

CHECK((const_reference operator() (size_type const i, size_type const j) const))
	Matrix<int> const & micr = mi;
  STATUS(micr(1,2));
	Matrix<int> & mir = mi;
  mir(1,2) = 22;
  STATUS(mir(1,2));
RESULT

CHECK((reference getValue(size_type const i, size_type const j)))
		STATUS(mi.getValue(1,2));
RESULT

CHECK((reference operator() (size_type const i, size_type const j)))
  // ???
RESULT

CHECK(void clear())
	Matrix<int> mi4(mi);
  STATUS("mi4:\n"<<mi4)
	mi4.clear();
  TEST_EQUAL(mi4.empty(),true);
  STATUS("mi4:\n"<<mi4)
	mi4 = mi;
  STATUS("mi4:\n"<<mi4)
RESULT

CHECK((void setValue(size_type const i, size_type const j, value_type value)))
	mi.setValue(1,1,18);
STATUS("mi:\n"<<mi)
RESULT

Matrix<int> mi5(4,5,6);

CHECK((Matrix(SizeType rows, SizeType cols, ValueType value = ValueType())))
  STATUS("mi5:\n"<<mi5)
	TEST_EQUAL(mi5.size(),20)
RESULT

CHECK(SizeType cols() const throw())
	TEST_EQUAL(mi5.rows(),4)
RESULT

CHECK(SizeType rows() const throw())
	TEST_EQUAL(mi5.cols(),5)
RESULT



Matrix<float> mf(6,7,8);

CHECK(SizeType colIndex(SizeType index) const)
	TEST_EQUAL(mf.colIndex(30),2)
RESULT

CHECK((SizeType const index(SizeType row, SizeType col) const))
	TEST_EQUAL(mf.index(5,5),40)
RESULT

CHECK(SizeType rowIndex(SizeType index) const)
  TEST_EQUAL(mf.rowIndex(30),4)
RESULT

CHECK((std::pair<Size,Size> const indexPair(Size index) const))
	std::pair<Size,Size> result = mf.indexPair(30);
  TEST_EQUAL(result.first,4)
	TEST_EQUAL(result.second,2)
RESULT

CHECK(~Matrix())
  // yeah, see above
RESULT

STATUS("\n\nDISCLAIMER: This ist not really a full test, but it shows that the thing compiles and you may inspect the STATUS() output manually ...\n\n\n")


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


