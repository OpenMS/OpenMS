// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

// Includes go here....
#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <sstream>

///////////////////////////

START_TEST(Matrix, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Matrix<int>* ptr = nullptr;
Matrix<int>* nullPointer = nullptr;
START_SECTION((Matrix()))
{
  ptr = new Matrix<int>;
  TEST_NOT_EQUAL(ptr, nullPointer);

  Matrix<int> mi1;
  TEST_EQUAL(mi1.size(), 0);
  TEST_EQUAL(mi1.cols(), 0);
  TEST_EQUAL(mi1.rows(), 0);

  for (auto & i : mi1.getEigenMatrix().reshaped())
  {
	TEST_EQUAL(i, i - 1); // this should not be executed on empty matrix
  }

  for (const auto & i : mi1.getEigenMatrix().reshaped())
  {
	TEST_EQUAL(i, i - 1); // this should not be executed on empty matrix
  }  
  STATUS("mi1:\n"<< mi1);
}
END_SECTION;

START_SECTION((~Matrix()))
{
	delete ptr;
}
END_SECTION;

Matrix<int> mi;

START_SECTION((void getEigenMatrix().resize(size_type i, size_type j)))
{
  mi.getEigenMatrix().resize(2,2);
  mi.getEigenMatrix().fill(3);
  mi.getEigenMatrix().resize(2,3);
  mi.getEigenMatrix().fill(7);
  STATUS("mi1:\n"<< mi);
  TEST_EQUAL(mi(0,0),7);
  TEST_EQUAL(mi(0,1),7);
  TEST_EQUAL(mi(0,2),7);
  TEST_EQUAL(mi(1,0),7);
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
  TEST_EQUAL(mi2(0,0),7);
  TEST_EQUAL(mi2(0,1),7);
  TEST_EQUAL(mi2(0,2),7);
  TEST_EQUAL(mi2(1,0),7);
  TEST_EQUAL(mi2(1,1),7);
  TEST_EQUAL(mi2(1,2),7);

  // test iterators and confirm column first order
  size_t row{}, col{};
  for (auto & i : mi2.getEigenMatrix().reshaped())
  {
	TEST_EQUAL(i, mi.getValue(row, col));
	++col;
	if (col == mi2.cols()) { col = 0; ++row;}
  }

  row = 0; col = 0;
  for (const auto & i : mi2.getEigenMatrix().reshaped())
  {
	TEST_EQUAL(i, mi.getValue(row, col));
	++col;
	if (col == mi2.cols()) { col = 0; ++row;}
  }  
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
	TEST_EQUAL(mi3(0,0),7);
	TEST_EQUAL(mi3(0,1),7);
	TEST_EQUAL(mi3(0,2),7);
	TEST_EQUAL(mi3(1,0),7);
	TEST_EQUAL(mi3(1,1),7);
	TEST_EQUAL(mi3(1,2),7);
}
END_SECTION

mi(1,1)=17;

START_SECTION((const_reference operator()(size_type const i, size_type const j) const))
{
  Matrix<int> const & micr = mi;
  STATUS("micr:\n"<<micr);
  TEST_EQUAL(micr(1,1),17);
}
END_SECTION

START_SECTION((reference operator()(size_type const i, size_type const j)))
{
	STATUS(mi(1,2));
	mi(1,2)=33;
	STATUS(mi(1,2));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr(1,2), 33);
}
END_SECTION

START_SECTION((reference operator() (size_type const i, size_type const j)))
{
	STATUS(mi(1,0));
	mi(1,0) = 44;
	STATUS(mi(1,0));
	Matrix<int> const & micr = mi;
	TEST_EQUAL(micr(1,0), 44);
}
END_SECTION

START_SECTION((void operator()(size_type const i, size_type const j) = value_type value))
{
  mi(1,1) = 18;
  STATUS("mi:\n" << mi);
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
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi3);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi4);(void) comparison);
	TEST_PRECONDITION_VIOLATED(bool comparison = (mi1==mi5);(void) comparison);
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
	myMatrix.setMatrix<double, 4 , 4>(test_matrix);
	for (size_t i=0; i<4; ++i)
	{
		for (size_t j=0; j<4; ++j)
		{
			TEST_EQUAL( myMatrix(i,j), test_matrix[i][j] )
		}
	}

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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


