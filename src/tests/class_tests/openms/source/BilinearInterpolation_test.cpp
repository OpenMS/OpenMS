// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

///////////////////////////

// This one is going to be tested.
#include <OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <vector>
#include <cstdlib>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/MATH/MathFunctions.h>

///////////////////////////

namespace OpenMS
{
	double rand01 ()
	{
		return double(rand()) / RAND_MAX;
	}
}

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( BilinearInterpolation, "$Id$" )

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef BilinearInterpolation < float, double > BIFD;

BIFD * bifd_ptr = nullptr;
BIFD * bifd_nullPointer = nullptr;

START_SECTION(BilinearInterpolation())
{
	BIFD bifd;
	bifd_ptr = new BIFD;
  TEST_NOT_EQUAL(bifd_ptr,bifd_nullPointer);
}
END_SECTION

START_SECTION(~BilinearInterpolation())
{
	delete bifd_ptr;
}
END_SECTION

START_SECTION(BilinearInterpolation& operator= ( BilinearInterpolation const & arg ))
{
	BIFD::ContainerType v;
	v.getEigenMatrix().resize(2,3);
	v(0,0) = 17;
	v(0,1) = 18.9;
	v(0,2) = 20.333;
	v(1,0) = -.1;
	v(1,1) = -.13;
	v(1,2) = -.001;

	BIFD bifd;
	bifd.setData(v);
	bifd.setMapping_0( 13, 230, 14, 250 );
	bifd.setMapping_1( 15, 2100, 17, 2900 );

	BIFD bifd2;

	// do it
	bifd2 = bifd;

	TEST_REAL_SIMILAR(bifd2.getScale_0(), 20);
	TEST_REAL_SIMILAR(bifd2.getScale_1(), 400);
	TEST_REAL_SIMILAR(bifd2.getOffset_0(), -30);
	TEST_REAL_SIMILAR(bifd2.getOffset_1(), -3900);
	TEST_REAL_SIMILAR(bifd2.getInsideReferencePoint_0(), 13);
	TEST_REAL_SIMILAR(bifd2.getOutsideReferencePoint_0(), 230);
	TEST_REAL_SIMILAR(bifd2.getInsideReferencePoint_1(), 15);
	TEST_REAL_SIMILAR(bifd2.getOutsideReferencePoint_1(), 2100);
	for ( UInt i = 0; i < v.rows(); ++i )
		for ( UInt j = 0; j < v.cols(); ++j )
			TEST_EQUAL(bifd2.getData()(i,j),v(i,j));
}
END_SECTION

START_SECTION(BilinearInterpolation( BilinearInterpolation const & arg ))
{
	BIFD::ContainerType v;
	v.getEigenMatrix().resize(2,3);
	v(0,0) = 17;
	v(0,1) = 18.9;
	v(0,2) = 20.333;
	v(1,0) = -.1;
	v(1,1) = -.13;
	v(1,2) = -.001;

	BIFD bifd;
	bifd.setData(v);
	bifd.setMapping_0( 13, 230, 14, 250 );
	bifd.setMapping_1( 15, 2100, 17, 2900 );

	// do it
	BIFD bifd2 = bifd;

	TEST_REAL_SIMILAR(bifd2.getScale_0(), 20);
	TEST_REAL_SIMILAR(bifd2.getScale_1(), 400);
	TEST_REAL_SIMILAR(bifd2.getOffset_0(), -30);
	TEST_REAL_SIMILAR(bifd2.getOffset_1(), -3900);
	TEST_REAL_SIMILAR(bifd2.getInsideReferencePoint_0(), 13);
	TEST_REAL_SIMILAR(bifd2.getOutsideReferencePoint_0(), 230);
	TEST_REAL_SIMILAR(bifd2.getInsideReferencePoint_1(), 15);
	TEST_REAL_SIMILAR(bifd2.getOutsideReferencePoint_1(), 2100);
	for ( UInt i = 0; i < v.rows(); ++i )
		for ( UInt j = 0; j < v.cols(); ++j )
			TEST_EQUAL(bifd2.getData()(i,j),v(i,j));
}
END_SECTION

START_SECTION(ContainerType& getData())
{
  BIFD bifd;
	bifd.getData().getEigenMatrix().resize(2,3);
	bifd.getData()(1,2) = 10012;
	bifd.getData()(0,0) = 10000;
	bifd.getData()(1,0) = 10010;

	BIFD const & bifd_cr(bifd);
	TEST_REAL_SIMILAR(bifd_cr.getData()(1,2),10012);
	TEST_REAL_SIMILAR(bifd_cr.getData()(0,0),10000);
	TEST_REAL_SIMILAR(bifd_cr.getData()(1,0),10010);
}
END_SECTION

START_SECTION(ContainerType const& getData() const)
{
  // see above,  ContainerType& getData()
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(template< typename SourceContainer > void setData( SourceContainer const & data ))
{
  BIFD bifd;
	bifd.getData().getEigenMatrix().resize(2,3);
	bifd.getData()(1,2) = 10012;
	bifd.getData()(0,0) = 10000;
	bifd.getData()(1,0) = 10010;

	BIFD const & bifd_cr(bifd);

	BIFD bifd2;
	bifd2.setData(bifd_cr.getData());

	TEST_EQUAL(bifd.getData(),bifd2.getData());

	BIFD bifd3;
	bifd3.getData().getEigenMatrix().resize(2,3);
	TEST_NOT_EQUAL(bifd.getData(),bifd3.getData());
}
END_SECTION


// Lots of methods with _0 and _1.
// I'll deal with the _0 case first,
// then copy/rename this for _1

START_SECTION((void setMapping_0( KeyType const & inside_low, KeyType const & outside_low, KeyType const & inside_high, KeyType const & outside_high )))
{
	BIFD bifd;
	bifd.setMapping_0(1,2,3,8);
	TEST_REAL_SIMILAR(bifd.getScale_0(),3);
	TEST_REAL_SIMILAR(bifd.getOffset_0(),-1);
	TEST_REAL_SIMILAR(bifd.getScale_1(),1);
	TEST_REAL_SIMILAR(bifd.getOffset_1(),0);
}
END_SECTION

START_SECTION((void setMapping_0( KeyType const & scale, KeyType const & inside_low, KeyType const & outside_low )))
{
	BIFD bifd;
	bifd.setMapping_0(3,1,2);
	TEST_REAL_SIMILAR(bifd.getScale_0(),3);
	TEST_REAL_SIMILAR(bifd.getOffset_0(),-1);
	TEST_REAL_SIMILAR(bifd.getScale_1(),1);
	TEST_REAL_SIMILAR(bifd.getOffset_1(),0);
}
END_SECTION

START_SECTION(void setOffset_0( KeyType const & offset ))
{
	BIFD bifd;
	bifd.setOffset_0(987);
	BIFD const& bifd_cr(bifd);
	TEST_REAL_SIMILAR(bifd_cr.getOffset_0(),987);
}
END_SECTION

START_SECTION(KeyType const& getOffset_0() const)
{
  // see above,  void setOffset_0( KeyType const & offset )
	NOT_TESTABLE;
}
END_SECTION
START_SECTION(void setScale_0( KeyType const & scale ))
{
	BIFD bifd;
	bifd.setScale_0(987);
	BIFD const& bifd_cr(bifd);
	TEST_REAL_SIMILAR(bifd_cr.getScale_0(),987);
}
END_SECTION

START_SECTION(KeyType const& getScale_0() const)
{
  // see above,  void setScale_0( KeyType const & scale )
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(KeyType const& getInsideReferencePoint_0() const)
{
 	BIFD bifd;
	bifd.setMapping_0(1,4,3,8);
	TEST_REAL_SIMILAR(bifd.getInsideReferencePoint_0(),1);
	TEST_REAL_SIMILAR(bifd.getOutsideReferencePoint_0(),4);
	TEST_REAL_SIMILAR(bifd.getInsideReferencePoint_1(),0);
	TEST_REAL_SIMILAR(bifd.getOutsideReferencePoint_1(),0);
}
END_SECTION

START_SECTION(KeyType const& getOutsideReferencePoint_0() const)
{
  // see above,  getInsideReferencePoint_0()
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(KeyType index2key_0( KeyType pos ) const)
{
 	BIFD bifd;
	bifd.setMapping_0(3,1,2);
	TEST_REAL_SIMILAR(bifd.index2key_0(0),-1);
	TEST_REAL_SIMILAR(bifd.index2key_1(0),0);
}
END_SECTION

START_SECTION(KeyType key2index_0( KeyType pos ) const)
{
 	BIFD bifd;
	bifd.setMapping_0(3,1,2);
	TEST_REAL_SIMILAR(bifd.key2index_0(-1),0);
	TEST_REAL_SIMILAR(bifd.key2index_1(0),0);
}
END_SECTION

START_SECTION(KeyType supportMax_0() const)
{
 	BIFD bifd;

	bifd.setMapping_0(3,1,2);
	bifd.setMapping_1(5,3,4);

	bifd.getData().getEigenMatrix().resize(2,3);

	TEST_REAL_SIMILAR(bifd.index2key_0(0),-1);
	TEST_REAL_SIMILAR(bifd.index2key_0(1),2);
	TEST_REAL_SIMILAR(bifd.supportMin_0(),-4);
	TEST_REAL_SIMILAR(bifd.supportMax_0(),5);

	TEST_REAL_SIMILAR(bifd.index2key_1(0),-11);
	TEST_REAL_SIMILAR(bifd.index2key_1(2),-1);
	TEST_REAL_SIMILAR(bifd.supportMin_1(),-16);
	TEST_REAL_SIMILAR(bifd.supportMax_1(),4);
}
END_SECTION

START_SECTION(KeyType supportMin_0() const)
{
  // see above,  supportMax_0
	NOT_TESTABLE;
}
END_SECTION



// here is the same stuff with _1 and _0 exchanged


START_SECTION((void setMapping_1( KeyType const & inside_low, KeyType const & outside_low, KeyType const & inside_high, KeyType const & outside_high )))
{
	BIFD bifd;
	bifd.setMapping_1(1,2,3,8);
	TEST_REAL_SIMILAR(bifd.getScale_1(),3);
	TEST_REAL_SIMILAR(bifd.getOffset_1(),-1);
	TEST_REAL_SIMILAR(bifd.getScale_0(),1);
	TEST_REAL_SIMILAR(bifd.getOffset_0(),0);
}
END_SECTION

START_SECTION((void setMapping_1( KeyType const & scale, KeyType const & inside_low, KeyType const & outside_low )))
{
	BIFD bifd;
	bifd.setMapping_1(3,1,2);
	TEST_REAL_SIMILAR(bifd.getScale_1(),3);
	TEST_REAL_SIMILAR(bifd.getOffset_1(),-1);
	TEST_REAL_SIMILAR(bifd.getScale_0(),1);
	TEST_REAL_SIMILAR(bifd.getOffset_0(),0);
}
END_SECTION

START_SECTION(void setOffset_1( KeyType const & offset ))
{
	BIFD bifd;
	bifd.setOffset_1(987);
	BIFD const& bifd_cr(bifd);
	TEST_REAL_SIMILAR(bifd_cr.getOffset_1(),987);
}
END_SECTION

START_SECTION(KeyType const& getOffset_1() const)
{
  // see above,  void setOffset_1( KeyType const & offset )
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(void setScale_1( KeyType const & scale ))
{
	BIFD bifd;
	bifd.setScale_1(987);
	BIFD const& bifd_cr(bifd);
	TEST_REAL_SIMILAR(bifd_cr.getScale_1(),987);
}
END_SECTION

START_SECTION(KeyType const& getScale_1() const)
{
  // see above,  void setScale_1( KeyType const & scale )
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(KeyType const& getInsideReferencePoint_1() const)
{
 	BIFD bifd;
	bifd.setMapping_1(1,4,3,8);
	TEST_REAL_SIMILAR(bifd.getInsideReferencePoint_1(),1);
	TEST_REAL_SIMILAR(bifd.getOutsideReferencePoint_1(),4);
	TEST_REAL_SIMILAR(bifd.getInsideReferencePoint_0(),0);
	TEST_REAL_SIMILAR(bifd.getOutsideReferencePoint_0(),0);
}
END_SECTION

START_SECTION(KeyType const& getOutsideReferencePoint_1() const)
{
  // see above,  getInsideReferencePoint_1()
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(KeyType index2key_1( KeyType pos ) const)
{
 	BIFD bifd;
	bifd.setMapping_1(3,1,2);
	TEST_REAL_SIMILAR(bifd.index2key_1(0),-1);
	TEST_REAL_SIMILAR(bifd.index2key_0(0),0);
}
END_SECTION

START_SECTION(KeyType key2index_1( KeyType pos ) const)
{
 	BIFD bifd;
	bifd.setMapping_1(3,1,2);
	TEST_REAL_SIMILAR(bifd.key2index_1(-1),0);
	TEST_REAL_SIMILAR(bifd.key2index_0(0),0);
}
END_SECTION

START_SECTION(KeyType supportMax_1() const)
{
 	BIFD bifd;

	bifd.setMapping_1(3,1,2);
	bifd.setMapping_0(5,3,4);

	bifd.getData().getEigenMatrix().resize(3,2);

	TEST_REAL_SIMILAR(bifd.index2key_1(0),-1);
	TEST_REAL_SIMILAR(bifd.index2key_1(1),2);
	TEST_REAL_SIMILAR(bifd.supportMin_1(),-4);
	TEST_REAL_SIMILAR(bifd.supportMax_1(),5);

	TEST_REAL_SIMILAR(bifd.index2key_0(0),-11);
	TEST_REAL_SIMILAR(bifd.index2key_0(2),-1);
	TEST_REAL_SIMILAR(bifd.supportMin_0(),-16);
	TEST_REAL_SIMILAR(bifd.supportMax_0(),4);
}
END_SECTION

START_SECTION(KeyType supportMin_1() const)
{
  // see above,  supportMax_1
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(bool empty() const)
{
	BIFD bifd;
	TEST_EQUAL(bifd.empty(),true);
	bifd.getData().getEigenMatrix().resize(1,2);
	TEST_EQUAL(bifd.empty(),false);
	bifd.getData().getEigenMatrix().resize(0,0);
	TEST_EQUAL(bifd.empty(),true);
	bifd.getData().getEigenMatrix().resize(1,2);
	TEST_EQUAL(bifd.empty(),false);
	bifd.getData().getEigenMatrix().resize(1,0);
	TEST_EQUAL(bifd.empty(),true);
	bifd.getData().getEigenMatrix().resize(1,2);
	TEST_EQUAL(bifd.empty(),false);
	bifd.getData().getEigenMatrix().resize(0,0);
	TEST_EQUAL(bifd.empty(),true);
	bifd.getData().getEigenMatrix().resize(2,2);
	TEST_EQUAL(bifd.empty(),false);
}
END_SECTION






START_SECTION((void addValue( KeyType arg_pos_0, KeyType arg_pos_1, ValueType arg_value )))
{

#define verbose(a)
	// #define verbose(a) a

	for ( int i = -50; i <= 100; ++i )
	{
		float p = i / 10.f;
		verbose(STATUS(i));

		for ( int j = -50; j <= 100; ++j )
		{
			float q = j / 10.f;
			verbose(STATUS("i: " << i));
			verbose(STATUS("j: " << j));

			BIFD bifd_small;
			bifd_small.getData().getEigenMatrix().resize(5,5);
			bifd_small.getData().getEigenMatrix().setZero();
			bifd_small.setMapping_0( 0, 0, 5, 5 );
			bifd_small.setMapping_1( 0, 0, 5, 5 );
			bifd_small.addValue( p, q, 100 );
			bifd_small.getData().getEigenMatrix() = 
				bifd_small.getData().getEigenMatrix().unaryExpr([](double val) {
					return Math::round(val);
				});
			verbose(STATUS("          " << bifd_small.getData()));

			BIFD bifd_big;
			bifd_big.getData().getEigenMatrix().resize(15,15);
			bifd_big.getData().getEigenMatrix().setZero();
			bifd_big.setMapping_0( 5, 0, 10, 5 );
			bifd_big.setMapping_1( 5, 0, 10, 5 );
			bifd_big.addValue( p, q, 100 );
			// Round the entries of the matrix in place
			bifd_big.getData().getEigenMatrix() = 
				bifd_big.getData().getEigenMatrix().unaryExpr([](double val) {
					return Math::round(val);
				});
			verbose(STATUS(bifd_big.getData()));

			BIFD::ContainerType big_submatrix;
			big_submatrix.getEigenMatrix().resize(5,5);
			for ( int m = 0; m < 5; ++m )
				for ( int n = 0; n < 5; ++n )
					big_submatrix(m,n)=bifd_big.getData()(m+5,n+5);

			TEST_EQUAL(bifd_small.getData(),big_submatrix);
		}
	}
#undef verbose

}
END_SECTION



START_SECTION((ValueType value( KeyType arg_pos_0, KeyType arg_pos_1 ) const))
{
#define verbose(a)
	// #define verbose(a) a

	// initialize random number generator (not platform independent, but at
	// least reproducible)
	srand(2007);

	TOLERANCE_ABSOLUTE(0.01);

	BIFD bifd_small;
	BIFD bifd_big;

	bifd_small.getData().getEigenMatrix().resize(5,5);
	bifd_small.getData().getEigenMatrix().setZero();
	bifd_big.getData().getEigenMatrix().resize(15,15);
	bifd_big.getData().getEigenMatrix().setZero();
	for ( int i = 0; i < 5; ++i )
	{
		for ( int j = 0; j < 5; ++j )
		{
			int num = int( floor(rand01() * 100.) );
			bifd_small.getData()(i,j) = num;
			bifd_big.getData()(i+5,j+5) = num;
		}
	}

	bifd_small.setMapping_0( 0, 0, 5, 5 );
	bifd_small.setMapping_1( 0, 0, 5, 5 );
	bifd_big.setMapping_0( 5, 0, 10, 5 );
	bifd_big.setMapping_1( 5, 0, 10, 5 );

	Matrix<int> interpolated(151,151);

	// you can view this as an image in PGM format
	verbose(std::cout << "P2\n151 151\n100\n" << std::endl);

	for ( int i = -50; i <= 100; ++i )
	{
		float p = i / 10.f;
		verbose(STATUS(i));

		for ( int j = -50; j <= 100; ++j )
		{
			float q = j / 10.f;
			verbose(STATUS("i: " << i));
			verbose(STATUS("j: " << j));

			TEST_REAL_SIMILAR(bifd_small.value(p,q),bifd_big.value(p,q));

			interpolated(i+50,j+50) = int(bifd_small.value(p,q));
		}
	}

	verbose(std::cout << interpolated << std::endl);


#undef verbose

}
END_SECTION


//-----------------------------------------------------------
//-----------------------------------------------------------

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
