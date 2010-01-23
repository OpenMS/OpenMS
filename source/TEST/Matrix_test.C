// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// $Authors: $
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
START_SECTION((Matrix()))
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
END_SECTION;

START_SECTION((~Matrix()))
{
	delete ptr;
}
END_SECTION;

Matrix<int> mi;

START_SECTION((void resize(size_type i, size_type j, value_type value = value_type())))
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
END_SECTION

START_SECTION((void resize(std::pair<Size, Size> const & size_pair, value_type value = value_type())))
{
	std::pair<Size, Size> const sizepair(2,2);
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
END_SECTION

START_SECTION((Matrix(const Matrix & source)))
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
END_SECTION

START_SECTION((Matrix& operator = (const Matrix & rhs)))
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
END_SECTION

mi(1,1)=17;

START_SECTION((const_reference getValue(size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr.getValue(1,1),17);
}
END_SECTION

START_SECTION((const_reference operator() (size_type const i, size_type const j) const))
{
	Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
	TEST_EQUAL(micr(1,1),17);
}
END_SECTION

START_SECTION((reference getValue(size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,2));
	mi.getValue(1,2)=33;
	STATUS(mi.getValue(1,2));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,2),33);
}
END_SECTION

START_SECTION((reference operator() (size_type const i, size_type const j)))
{
	STATUS(mi.getValue(1,0));
	mi.getValue(1,0)=44;
	STATUS(mi.getValue(1,0));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr.getValue(1,0),44);
}
END_SECTION

START_SECTION((void clear()))
{
	Matrix<int> mi4(mi);
  STATUS("mi4:\n"<<mi4);
	mi4.clear();
  STATUS("mi4:\n"<<mi4);
  TEST_EQUAL(mi4.empty(),true);
}
END_SECTION

START_SECTION((void setValue(size_type const i, size_type const j, value_type value)))
{
	mi.setValue(1,1,18);
	STATUS("mi:\n"<<mi);
	TEST_EQUAL(mi(1,1),18);
}
END_SECTION;

Matrix<int> mi5(4,5,6);

START_SECTION((Matrix(const SizeType rows, const SizeType cols, ValueType value = ValueType())))
{
  STATUS("mi5:\n"<<mi5);
	TEST_EQUAL(mi5.size(),20);
}
END_SECTION;

START_SECTION((SizeType cols() const))
{
	TEST_EQUAL(mi5.rows(),4);
}
END_SECTION;

START_SECTION((SizeType rows() const))
{
	TEST_EQUAL(mi5.cols(),5);
}
END_SECTION;

Matrix<float> mf(6,7,8);

START_SECTION((SizeType colIndex(SizeType index) const))
{
	TEST_EQUAL(mf.colIndex(30),2);
}
END_SECTION;

START_SECTION((SizeType const index(SizeType row, SizeType col) const))
{
	TEST_EQUAL(mf.index(5,5),40);
}
END_SECTION;

START_SECTION((SizeType rowIndex(SizeType index) const))
{
  TEST_EQUAL(mf.rowIndex(30),4);
}
END_SECTION;

START_SECTION((std::pair<Size,Size> const indexPair(Size index) const))
{
	std::pair<Size,Size> result = mf.indexPair(30);
  TEST_EQUAL(result.first,4);
	TEST_EQUAL(result.second,2);
}
END_SECTION

START_SECTION((std::pair<Size,Size> sizePair() const))
{
	Matrix<float> const mf(6,7,8);
	TEST_EQUAL(mf.sizePair().first,6);
	TEST_EQUAL(mf.sizePair().second,7);
}
END_SECTION

START_SECTION((bool operator==(Matrix const &rhs) const))
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
	TEST_PRECONDITION_VIOLATED(mi1==mi3);
	TEST_PRECONDITION_VIOLATED(mi1==mi4);
	TEST_PRECONDITION_VIOLATED(mi1==mi5);
}
END_SECTION

START_SECTION((bool operator<(Matrix const &rhs) const))
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
	TEST_PRECONDITION_VIOLATED(mi1==mi3);
	TEST_PRECONDITION_VIOLATED(mi1==mi4);
	TEST_PRECONDITION_VIOLATED(mi1==mi5);
}
END_SECTION

START_SECTION((template <int ROWS, int COLS> void setMatrix(const ValueType matrix[ROWS][COLS])))
{
	double test_matrix[4][4] = {
		{0, 2.5,   3, 0.1},
		{0,   1, 5.9, 0.2},
		{0,   2, 5.6, 0.1},
		{0,   2,   3, 0.1}
	};

	Matrix<double> myMatrix;
	myMatrix.setMatrix<4,4>(test_matrix);
	for (size_t i=0;i<4;++i)
	{
		for (size_t j=0;j<4;++j)
		{
			TEST_EQUAL( myMatrix(i,j), test_matrix[i][j] )
		}
	}

}
END_SECTION

START_SECTION((gsl_matrix * toGslMatrix()))
{
	Matrix<double> mi(2,3,6);
	mi(1,2)=112;
	mi(0,0)=100;
	mi(1,1)=111;
	mi(0,2)=103;
	gsl_matrix* gsl_m = mi.toGslMatrix();
	for (size_t i=0;i<2;++i)
	{
		for (size_t j=0;j<3;++j)
		{
			TEST_EQUAL(mi(i,j),gsl_matrix_get (gsl_m, i, j))
		}
	}
	gsl_matrix_free (gsl_m);
}
END_SECTION

START_SECTION((template <typename  Value > std::ostream & operator<<(std::ostream &os, const Matrix< Value > &matrix)))
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
END_SECTION

#if 0
// actually seems to *generate* a warning for me! - Clemens

START_SECTION((OPENMS_DLLAPI gsl_matrix * toGslMatrix()))
{
	NOT_TESTABLE // avoid warning
}
END_SECTION

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


