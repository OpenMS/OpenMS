// -*- mode: C++; tab-width: 2; -*-
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

///////////////////////////

// This one is going to be tested.
#include <OpenMS/MATH/MISC/BilinearInterpolation.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <vector>

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

namespace OpenMS
{
	template < typename T >
	std::ostream & operator << ( std::ostream & os, std::vector < T > const & cont )
	{
		for ( typename std::vector<T>::const_iterator iter = cont.begin(); iter != cont.end(); ++iter )
		{
			os << ' ' << *iter;
		}
		return os;
	}
}

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( BilinearInterpolation, "$Id$" )

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
CHECK( typedefs )
{
	typedef BilinearInterpolation < float, double > BIFD;
	BIFD::ValueType     * value;
	BIFD::KeyType       * key;
	BIFD::ContainerType * container;
	BIFD::ContainerType::value_type * containerValue;
	value = 0;
	key = 0;
	container = 0;
	containerValue = 0;
}
RESULT

typedef BilinearInterpolation < float, double > BIFD;

//-----------------------------------------------------------
CHECK( default constructor )
{
	BIFD bifd;
}
RESULT

#if 0

//-----------------------------------------------------------
CHECK( value() )
{
	typedef BilinearInterpolation < float, double > BIFD;
	BIFD bifd;
	int const rows = 5, cols = 5;
	bifd.getData().resize(rows, cols);
	bifd.addValue( 0, 0, 0 );
	bifd.addValue( 0.5, 0, 10 );
	STATUS("dump of matrix bifd follows");
	std::ofstream original_pgm("original.pgm");
	bifd.getData().writePGM(std::cout,100,0);
	bifd.getData().writePGM(original_pgm,100,0);
	STATUS("dump of matrix bifd finished");

	for ( float x = -2; x <= rows + 1; x += .25 )
	{
		std::cout << "row:" << x << ":\t";
		for ( float y = -2; y <= cols + 1; y += .25 )
		{
			std::cout << bifd.value(x,y) << ' ' ;
		}
		std::cout << '\n';
	}

	BIFD resampled;
	resampled.getData().resize((rows+3)*4+1,(cols+3)*4+1);
	bifd.setMapping_0 ( 4, 0, 8 );
	bifd.setMapping_1 ( 4, 0, 8 );
	for ( int p = 0; p != resampled.getData().size(); ++p )
	{
		const int i = resampled.getData().rowIndex(p);
		const int j = resampled.getData().colIndex(p);
		// std::cout << i <<' ' << j << ' ';
		resampled.getData()(i,j) = bifd.value(i,j);
	}
	STATUS("dump of matrix resampled follows");
	std::ofstream resampled_pgm("resampled.pgm");
	resampled.getData().writePGM(std::cout,1000,100);
	resampled.getData().writePGM(resampled_pgm,1000,100);
	STATUS("dump of matrix resampled finished");

}
RESULT
//-----------------------------------------------------------
CHECK( simple painting application ... \n )
{
	typedef BilinearInterpolation < float, double > BIFD;
	BIFD bifd;
	BIFD resampled;
	int const rows = 5, cols = 5;
	bifd.getData().resize(rows, cols);
	float x = 0, y = 0, it = 0;
	for (;;)
	{
		std::cout << "Enter x y intensity ...\n";
		std::cin >> x >> y >> it;
		if ( ! std::cin ) break;
		std::cout << "Read " << x <<' '<< y <<' '<< it <<'\n';
		bifd.addValue( x, y, it );
		std::ofstream original_pgm("original.pgm");
		bifd.getData().writePGM(original_pgm,-100,0);
		bifd.getData().writePGM(std::cout,100,0);
		original_pgm.close();
#if 0
		resampled.getData().clear();
		resampled.getData().resize((rows+3)*4+1,(cols+3)*4+1);
		bifd.setMapping_0 ( 4, 0, 8 );
		bifd.setMapping_1 ( 4, 0, 8 );
		for ( int p = 0; p != resampled.getData().size(); ++p )
		{
			const int i = resampled.getData().rowIndex(p);
			const int j = resampled.getData().colIndex(p);
			resampled.getData()(i,j) = bifd.value(i,j);
		}
		std::ofstream resampled_pgm("resampled.pgm");
		resampled.getData().writePGM(resampled_pgm,-100,0);
		resampled_pgm.close();
#endif
	}
	STATUS("Stop.");

}
RESULT

#endif

// Please do not remove the {} inside the CHECK...RESULT blocks.
// Emacs will completely mess up the indentation otherwise.

CHECK(BilinearInterpolation& operator= ( BilinearInterpolation const & arg ))
{
	BIFD::ContainerType v;
	v.resize(2,3);
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

	TEST_REAL_EQUAL(bifd2.getScale_0(), 20);
	TEST_REAL_EQUAL(bifd2.getScale_1(), 400);
	TEST_REAL_EQUAL(bifd2.getOffset_0(), -30);
	TEST_REAL_EQUAL(bifd2.getOffset_1(), -3900);
	TEST_REAL_EQUAL(bifd2.getInsideReferencePoint_0(), 13);
	TEST_REAL_EQUAL(bifd2.getOutsideReferencePoint_0(), 230);
	TEST_REAL_EQUAL(bifd2.getInsideReferencePoint_1(), 15);
	TEST_REAL_EQUAL(bifd2.getOutsideReferencePoint_1(), 2100);
	for ( UInt i = 0; i < v.rows(); ++i )
		for ( UInt j = 0; j < v.cols(); ++j )
			TEST_EQUAL(bifd2.getData()(i,j),v(i,j));
}
RESULT

CHECK(BilinearInterpolation( BilinearInterpolation const & arg ))
{
	BIFD::ContainerType v;
	v.resize(2,3);
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

	TEST_REAL_EQUAL(bifd2.getScale_0(), 20);
	TEST_REAL_EQUAL(bifd2.getScale_1(), 400);
	TEST_REAL_EQUAL(bifd2.getOffset_0(), -30);
	TEST_REAL_EQUAL(bifd2.getOffset_1(), -3900);
	TEST_REAL_EQUAL(bifd2.getInsideReferencePoint_0(), 13);
	TEST_REAL_EQUAL(bifd2.getOutsideReferencePoint_0(), 230);
	TEST_REAL_EQUAL(bifd2.getInsideReferencePoint_1(), 15);
	TEST_REAL_EQUAL(bifd2.getOutsideReferencePoint_1(), 2100);
	for ( UInt i = 0; i < v.rows(); ++i )
		for ( UInt j = 0; j < v.cols(); ++j )
			TEST_EQUAL(bifd2.getData()(i,j),v(i,j));
}
RESULT

CHECK((BilinearInterpolation( KeyType scale_0 = 1., KeyType offset_0 = 0., KeyType scale_1 = 1., KeyType offset_1 = 0. )))
{
  // ???
}
RESULT

CHECK(ContainerType const& getData() const)
{
  // ???
}
RESULT

CHECK(ContainerType& getData())
{
  // ???
}
RESULT

// Lots of methods with _0 and _1.
// I'll deal with the _0 case first,
// then copy/rename this for _1

CHECK((void setMapping_0( KeyType const & inside_low, KeyType const & outside_low, KeyType const & inside_high, KeyType const & outside_high )))
{
  // ???
}
RESULT

CHECK((void setMapping_0( KeyType const & scale, KeyType const & inside, KeyType const & outside )))
{
  // ???
}
RESULT

CHECK(void setOffset_0( KeyType const & offset ))
{
  // ???
}
RESULT

CHECK(void setScale_0( KeyType const & scale ))
{
  // ???
}
RESULT

CHECK(KeyType const& getInsideReferencePoint_0() const)
{
  // ???
}
RESULT

CHECK(KeyType const& getOffset_0() const)
{
  // ???
}
RESULT
CHECK(KeyType const& getOutsideReferencePoint_0() const)
{
  // ???
}
RESULT

CHECK(KeyType const& getScale_0() const)
{
  // ???
}
RESULT

CHECK(KeyType index2key_0( KeyType pos ) const)
{
  // ???
}
RESULT

CHECK(KeyType key2index_0( KeyType pos ) const)
{
  // ???
}
RESULT

CHECK(KeyType supportMax_0() const)
{
  // ???
}
RESULT

CHECK(KeyType supportMin_0() const)
{
  // ???
}
RESULT

CHECK((ValueType value( KeyType arg_pos_0, KeyType arg_pos_1 ) const))
{
  // ???
}
RESULT

CHECK(bool empty() const)
{
  // ???
}
RESULT

CHECK(template< typename SourceContainer > void setData( SourceContainer const & data ))
{
  // ???
}
RESULT

CHECK((void addValue( KeyType arg_pos_0, KeyType arg_pos_1, ValueType arg_value )))
{

#define verbose(a)
	// #define verbose(a) a

	for ( int i = -50; i <= 100; ++i )
	{
		float p = i / 10.;
		verbose(STATUS(i));

		for ( int j = -50; j <= 100; ++j )
		{
			float q = j / 10.;
			verbose(STATUS("i: " << i));
			verbose(STATUS("j: " << j));

			BIFD bifd_small;
			bifd_small.getData().resize(5,5,0);
			bifd_small.setMapping_0( 0, 0, 5, 5 );
			bifd_small.setMapping_1( 0, 0, 5, 5 );
			bifd_small.addValue( p, q, 100 );
			for ( BIFD::ContainerType::iterator iter = bifd_small.getData().begin();
						iter != bifd_small.getData().end();
						++iter
					) *iter = round(*iter);
			verbose(STATUS("          " << bifd_small.getData()));

			BIFD bifd_big;
			bifd_big.getData().resize(15,15,0);
			bifd_big.setMapping_0( 5, 0, 10, 5 );
			bifd_big.setMapping_1( 5, 0, 10, 5 );
			bifd_big.addValue( p, q, 100 );
			for ( BIFD::ContainerType::iterator iter = bifd_big.getData().begin();
						iter != bifd_big.getData().end();
						++iter
					) *iter = round(*iter);
			verbose(STATUS(bifd_big.getData()));

			BIFD::ContainerType big_submatrix;
			big_submatrix.resize(5,5);
			for ( int m = 0; m < 5; ++m )
				for ( int n = 0; n < 5; ++n )
					big_submatrix(m,n)=bifd_big.getData()(m+5,n+5);

			TEST_EQUAL(bifd_small.getData(),big_submatrix);
		}
	}
}
RESULT

CHECK(~BilinearInterpolation())
{
  // ???
}
RESULT



//-----------------------------------------------------------
//-----------------------------------------------------------

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
