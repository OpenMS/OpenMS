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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

// Includes go here....
#include <OpenMS/DATASTRUCTURES/Matrix.h>

///////////////////////////

START_TEST(Matrix, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Matrix<int>* ptr = 0;
CHECK(Matrix<int>())
	ptr = new Matrix<int>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Matrix<int>())
	delete ptr;
RESULT

CHECK(Matrix())
{
  Matrix<int> mi1;
	TEST_EQUAL(mi1.size(),0);
	TEST_EQUAL(mi1.cols(),0);
	TEST_EQUAL(mi1.rows(),0);
	TEST_EQUAL(mi1.empty(),true);
  STATUS("mi1:\n"<< mi1);
}
RESULT

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

CHECK(Matrix(const Matrix & source))
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

CHECK(Matrix& operator = (const Matrix & rhs))
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

CHECK(void clear())
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
{
	std::pair<Size,Size> result = mf.indexPair(30);
  TEST_EQUAL(result.first,4);
	TEST_EQUAL(result.second,2);
}
RESULT
	
CHECK(writePGM())
{
	int feep[] =
		{ // The "feep" image is from http://netpbm.sourceforge.net/doc/pgm.html
			0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //
			0,  3,  3,  3,  3,  0,  0,  7,  7,  7,  7,  0,  0, 11, 11, 11, 11,  0,  0, 15, 15, 15, 15,  0, //
			0,  3,  0,  0,  0,  0,  0,  7,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0, 15,  0,  0, 15,  0, //
			0,  3,  3,  3,  0,  0,  0,  7,  7,  7,  0,  0,  0, 11, 11, 11,  0,  0,  0, 15, 15, 15, 15,  0, //
			0,  3,  0,  0,  0,  0,  0,  7,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0, 15,  0,  0,  0,  0, //
			0,  3,  0,  0,  0,  0,  0,  7,  7,  7,  7,  0,  0, 11, 11, 11, 11,  0,  0, 15,  0,  0,  0,  0, //
			0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //
		};

	Matrix<int> matrix;
	matrix.resize(7,24);
	matrix.assign( feep, feep + 24*7 );
	Matrix<int> const & matrixcr(matrix);
	
	{ // using defaults (black background, no scaling)
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output);
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_d.pgm");
	}
	{ // black background, with comment
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,15,0,1.,false,"One comment line\nAnother comment line.");
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_c.pgm");
	}
	{ // white background, reduced gray scale
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,4,0,1.,true);
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_r.pgm");
	}
	{ // using 1 bit (get "  ep" only)
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,1,1./15,1.,false,"binary feep (is \"  ep\" only)");
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_1.pgm");
	}
	{ // using 1 bit, reverse video (get "  ep" only)
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,1,1./15,1.,true,"binary feep (is \"  ep\" only)");
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_1r.pgm");
	}
	{ // no scaling, but reverse video
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,0,1,1.,true);
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_nr.pgm");
	}
	{ // using 16 bit
		std::string feep_pgm;
		NEW_TMP_FILE(feep_pgm);
		std::ofstream output (feep_pgm.c_str());
		matrixcr.writePGM(output,(1<<16)-1);
		output.close();
		TEST_FILE(feep_pgm.c_str(),"data/Matrix_test_16.pgm");
	}
#if 0
#endif
}
RESULT
	

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


