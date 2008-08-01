// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------
//


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

// Includes go here....
#include <OpenMS/DATASTRUCTURES/Matrix.h>

///////////////////////////

START_TEST(Matrix, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Matrix<int>* ptr = 0;
CHECK((Matrix()))
{
	ptr = new Matrix<int>;
	TEST_NOT_EQUAL(ptr, 0);

  Matrix<int> mi1;
	TEST_EQUAL(mi1.size(),0);
	TEST_EQUAL(mi1.cols(),0);
	TEST_EQUAL(mi1.rows(),0);
	TEST_EQUAL(mi1.empty(),true);
  STATUS("mi1:\n"<< mi1);
}
RESULT;

CHECK((~Matrix()))
{
	delete ptr;
}
RESULT;

Matrix<int> mi;

CHECK((void resize(size_type i, size_type j, value_type value = value_type())))
{
	mi.resize(2,2,3);
  STATUS("mi1:\n"<< mi);
	mi.resize(2,3,7);
  STATUS("mi1:\n"<< mi);
	TEST_EQUAL(mi(0,0),3);
	TEST_EQUAL(mi(0,1),3);
	TEST_EQUAL(mi(0,2),3);
	TEST_EQUAL(mi(1,0),3);
	TEST_EQUAL(mi(1,1),7);
	TEST_EQUAL(mi(1,2),7);
}
RESULT

CHECK((void resize(std::pair<UInt,UInt> const & size_pair, value_type value = value_type())))
{
	std::pair<UInt,UInt> const sizepair(2,2);
	mi.resize(sizepair,3);
  STATUS("mi1:\n"<< mi);
	mi.resize(2,3,7);
  STATUS("mi1:\n"<< mi);
	TEST_EQUAL(mi(0,0),3);
	TEST_EQUAL(mi(0,1),3);
	TEST_EQUAL(mi(0,2),3);
	TEST_EQUAL(mi(1,0),3);
	TEST_EQUAL(mi(1,1),7);
	TEST_EQUAL(mi(1,2),7);
}
RESULT

CHECK((Matrix(const Matrix & source)))
{
  Matrix<int> mi2(mi);
  STATUS("mi2:\n"<< mi2);
	TEST_EQUAL(mi2.cols(),3);
	TEST_EQUAL(mi2.rows(),2);
	TEST_EQUAL(mi2(0,0),3);
	TEST_EQUAL(mi2(0,1),3);
	TEST_EQUAL(mi2(0,2),3);
	TEST_EQUAL(mi2(1,0),3);
	TEST_EQUAL(mi2(1,1),7);
	TEST_EQUAL(mi2(1,2),7);
}
RESULT

CHECK((Matrix& operator = (const Matrix & rhs)))
{
	Matrix<int> mi3;
  STATUS("mi3:\n"<<mi3);
	mi3=mi;
  STATUS("mi3:\n"<<mi3);
	TEST_EQUAL(mi3.cols(),3);
	TEST_EQUAL(mi3.rows(),2);
	TEST_EQUAL(mi3(0,0),3);
	TEST_EQUAL(mi3(0,1),3);
	TEST_EQUAL(mi3(0,2),3);
	TEST_EQUAL(mi3(1,0),3);
	TEST_EQUAL(mi3(1,1),7);
	TEST_EQUAL(mi3(1,2),7);
}
RESULT

mi(1,1)=17;

CHECK((const_reference getValue(size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr.getValue(1,1),17);
}
RESULT

CHECK((const_reference operator() (size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr(1,1),17);
}
RESULT

CHECK((reference getValue(size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,2));
	mi.getValue(1,2)=33;
	STATUS(mi.getValue(1,2));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,2),33);
}
RESULT

CHECK((reference operator() (size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,0));
	mi.getValue(1,0)=44;
	STATUS(mi.getValue(1,0));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,0),44);
}
RESULT

CHECK((void clear()))
{
	Matrix<int> mi4(mi);
  STATUS("mi4:\n"<<mi4);
	mi4.clear();
  STATUS("mi4:\n"<<mi4);
  TEST_EQUAL(mi4.empty(),true);
}
RESULT

CHECK((void setValue(size_type const i, size_type const j, value_type value)))
{
	mi.setValue(1,1,18);
	STATUS("mi:\n"<<mi);
	TEST_EQUAL(mi(1,1),18);
}
RESULT;

Matrix<int> mi5(4,5,6);

CHECK((Matrix(SizeType rows, SizeType cols, ValueType value = ValueType())))
{
  STATUS("mi5:\n"<<mi5);
	TEST_EQUAL(mi5.size(),20);
}
RESULT;

CHECK((SizeType cols() const))
{
	TEST_EQUAL(mi5.rows(),4);
}
RESULT;

CHECK((SizeType rows() const))
{
	TEST_EQUAL(mi5.cols(),5);
}
RESULT;

Matrix<float> mf(6,7,8);

CHECK((SizeType colIndex(SizeType index) const))
{
	TEST_EQUAL(mf.colIndex(30),2);
}
RESULT;

CHECK((SizeType const index(SizeType row, SizeType col) const))
{
	TEST_EQUAL(mf.index(5,5),40);
}
RESULT;

CHECK((SizeType rowIndex(SizeType index) const))
{
  TEST_EQUAL(mf.rowIndex(30),4);
}
RESULT;

CHECK((std::pair<UInt,UInt> const indexPair(UInt index) const))
{
	std::pair<UInt,UInt> result = mf.indexPair(30);
  TEST_EQUAL(result.first,4);
	TEST_EQUAL(result.second,2);
}
RESULT

CHECK((std::pair<UInt,UInt> sizePair() const))
{
	Matrix<float> const mf(6,7,8);
	TEST_EQUAL(mf.sizePair().first,6);
	TEST_EQUAL(mf.sizePair().second,7);
}
RESULT

CHECK((bool operator == ( Matrix const & rhs ) const throw (Exception::Precondition)))
{
	Matrix<int> mi1(4,5,6);
	mi1(2,3)=17;
	Matrix<int> mi2(4,5,6);
	TEST_NOT_EQUAL(mi1,mi2);
	mi1(2,3)=6;
	TEST_EQUAL(mi1,mi2);

	Matrix<int> mi3(5,4,6);
	Matrix<int> mi4(4,4,6);
	Matrix<int> mi5(5,5,6);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi3);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi4);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi5);
}
RESULT

CHECK((bool operator < ( Matrix const & rhs ) const throw (Exception::Precondition)))
{
	Matrix<int> mi1(4,5,6);
	TEST_EQUAL(mi1<mi1,false);
	mi1(2,3)=17;
	TEST_EQUAL(mi1<mi1,false);
	Matrix<int> mi2(4,5,6);
	TEST_EQUAL(mi1<mi2,false);
	TEST_EQUAL(mi2<mi1,true);
	mi2(2,3)=18;
	TEST_EQUAL(mi1<mi2,true);

	Matrix<int> mi3(5,4,6);
	Matrix<int> mi4(4,4,6);
	Matrix<int> mi5(5,5,6);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi3);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi4);
	TEST_EXCEPTION(Exception::Precondition,mi1==mi5);
}
RESULT

CHECK((template <typename Value> std::ostream & operator<<(std::ostream &os, const Matrix< Value > &matrix)))
{
	Matrix<int> mi(2,3,6);
	mi(1,2)=112;
	mi(0,0)=100;
	mi(1,1)=111;
	mi(0,2)=103;
	std::ostringstream os;
	os << mi;
	// Uh, finally I got the whitespace right
	char matrix_dump[] =
	"   100      6    103 \n"
	"     6    111    112 \n";
	TEST_EQUAL(os.str(),matrix_dump);
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


