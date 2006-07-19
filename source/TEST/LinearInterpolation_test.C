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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

///////////////////////////

// This one is going to be tested.
#include <OpenMS/MATH/MISC/LinearInterpolation.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <vector>

#include <OpenMS/CONCEPT/ClassTest.h>


///////////////////////////

namespace OpenMS
{
	// no extra stuff required for this test
}

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( LinearInterpolation, "$Id$" )

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
CHECK( typedefs )
{
	typedef LinearInterpolation < float, double > SMDFD;
	SMDFD::ValueType     * value;
	SMDFD::KeyType       * key;
	SMDFD::ContainerType * container;
	SMDFD::ContainerType::value_type * containerValue;
	value = 0;
	key = 0;
	container = 0;
	containerValue = 0;
}
RESULT
//-----------------------------------------------------------
CHECK( constructors and accessors )
{
	typedef LinearInterpolation < float, double > SMDFD;
	SMDFD smdfd0;
	SMDFD smdfd1 ( 1.125 );
	SMDFD smdfd2 ( 1.125, 3.5 );

	TEST_EQUAL ( smdfd0.getScale(), 1 )
	TEST_EQUAL ( smdfd0.getOffset(), 0 )
		smdfd0.setScale ( 2.5 );
		smdfd0.setOffset ( -1 );
	TEST_EQUAL ( smdfd0.getScale(), 2.5 )
	TEST_EQUAL ( smdfd0.getOffset(), -1 )

	TEST_EQUAL ( smdfd1.getScale(), 1.125 )
	TEST_EQUAL ( smdfd1.getOffset(), 0 )

	TEST_EQUAL ( smdfd2.getScale(), 1.125 )
	TEST_EQUAL ( smdfd2.getOffset(), 3.5 )

	TEST_EQUAL ( smdfd0.getData().size(), 0 )
	TEST_EQUAL ( smdfd1.getData().size(), 0 )
	TEST_EQUAL ( smdfd2.getData().size(), 0 )

}
RESULT
//-----------------------------------------------------------
CHECK( supportMin() and supportMax() )
{
	typedef LinearInterpolation < float, double > SMDFD;

	SMDFD smdfd2 ( 1.125, 3.5 );

	TEST_REAL_EQUAL ( smdfd2.supportMin(), 3.5 )
	TEST_REAL_EQUAL ( smdfd2.supportMax(), 3.5 )

	smdfd2.getData().push_back(1);

	TEST_REAL_EQUAL ( smdfd2.supportMin(), 3.5-1.125 )
	TEST_REAL_EQUAL ( smdfd2.supportMax(), 3.5+1.125 )

	SMDFD smdfd3;

	TEST_EQUAL ( smdfd3.empty(), true )

	smdfd3.setScale ( smdfd2.getScale() );
	smdfd3.setOffset ( smdfd2.getOffset() );
	smdfd3.setData( smdfd2.getData() );

	TEST_EQUAL ( smdfd3.empty(), false )

	TEST_REAL_EQUAL ( smdfd3.supportMin(), 3.5-1.125 )
	TEST_REAL_EQUAL ( smdfd3.supportMax(), 3.5+1.125 )

}
RESULT
//-----------------------------------------------------------
CHECK( copy constructor )
{
	typedef LinearInterpolation < float, double > SMDFD;

	SMDFD smdfd2 ( 1.125, 3.5 );

	TEST_REAL_EQUAL ( smdfd2.supportMin(), 3.5 )
	TEST_REAL_EQUAL ( smdfd2.supportMax(), 3.5 )

	smdfd2.getData().push_back(1);

	TEST_REAL_EQUAL ( smdfd2.supportMin(), 3.5-1.125 )
	TEST_REAL_EQUAL ( smdfd2.supportMax(), 3.5+1.125 )

	SMDFD smdfd3 ( smdfd2 );

	TEST_EQUAL ( smdfd3.empty(), false )

	TEST_REAL_EQUAL ( smdfd3.supportMin(), 3.5-1.125 )
	TEST_REAL_EQUAL ( smdfd3.supportMax(), 3.5+1.125 )

}
RESULT
//-----------------------------------------------------------
CHECK( value() and key2index() and index2key() )
{
	typedef LinearInterpolation < float, double > SMDFD;

	SMDFD smdfd0;

	double values[] = { 1, 2, 0, 1 };
	int const num_values = sizeof (values) / sizeof (*values);
	std::copy ( values, values + num_values, std::back_inserter ( smdfd0.getData() ) );

	TEST_EQUAL ( (int) smdfd0.getData().size(), num_values );

	for ( int i = 0; i < num_values; ++i )
	{
		TEST_EQUAL ( smdfd0.value ( i ), values[i] );
	}

	double inter_values[] =
		{ 0, 0.00, 0.00, 0.00,
			0, 0.25, 0.50, 0.75,
			1, 1.25, 1.50, 1.75,
			2, 1.50, 1.00, 0.50,
			0, 0.25, 0.50, 0.75,
			1, 0.75, 0.50, 0.25,
			0, 0.00, 0.00, 0.00,
			0
		};
	
	for ( int i = 0; i < num_values+4; ++i )
	{
		TEST_REAL_EQUAL ( smdfd0.value ( i-2 ), inter_values[4*i] );
	}

	int const num_inter_values = sizeof (inter_values) / sizeof (*inter_values);
	for ( int i = 0; i < num_inter_values; ++i )
	{
		TEST_REAL_EQUAL ( smdfd0.value ( (i-8.)/4. ), inter_values[i] );
	}

	SMDFD smdfd1 (smdfd0);

	double const scale = 1;
	double const offset = 100;
	smdfd1.setScale( scale );
	smdfd1.setOffset( offset );

	for ( int i = -8; i < num_inter_values-8; ++i )
	{
		double pos = i/4.;
		TEST_REAL_EQUAL ( smdfd1.key2index(smdfd1.index2key(pos)), pos );
	}

	for ( int i = -8; i < num_inter_values-8; ++i )
	{
		double pos = i/4.; 
		TEST_REAL_EQUAL ( smdfd1.value ( pos*scale+offset ), smdfd0.value ( pos ) ) ;
	}

}
RESULT
//-----------------------------------------------------------
CHECK( derivative()  )
{
	typedef LinearInterpolation < float, double > SMDFD;

	SMDFD smdfd0;

	double values[] = { 1, 2, 0, 1 };
	int const num_values = sizeof (values) / sizeof (*values);
	std::copy ( values, values + num_values, std::back_inserter ( smdfd0.getData() ) );

	TEST_EQUAL ( (int) smdfd0.getData().size(), num_values );

	for ( int i = 0; i < num_values; ++i )
	{
		TEST_EQUAL ( smdfd0.value ( i ), values[i] );
	}

	double inter_values[] =         // left .. (derivative) .. right
		{ +0.00,  0.00,  0.00,  0.25, // 0 .. (0) .. 0
			+0.50,  0.75,  1.00,  1.00, // 0 .. (1) .. 1
			+1.00,  1.00,  1.00,  0.25, // 1 .. (1) .. 2
			-0.50, -1.25, -2.00, -1.25, // 2 .. (-2) .. 0
			-0.50,  0.25,  1.00,  0.50, // 0 .. (1) .. 1
			+0.00, -0.50, -1.00, -0.75, // 1 .. (-1) .. 0
			-0.50, -0.25,  0.00,  0.00, // 0 .. (0) .. 0
			0
		};
	
	int const num_inter_values = sizeof (inter_values) / sizeof (*inter_values);
	for ( int i = -8; i < num_inter_values-8; ++i )
	{
		double key = i/4.;
		int index = i+8;
		STATUS( "key:" << key << "  index:" << index << '\n')
		TEST_REAL_EQUAL ( smdfd0.derivative ( key ), inter_values[index] );
	}

}
RESULT
//-----------------------------------------------------------
CHECK( setMapping() and getInsideReferencePoint() and getOutsideReferencePoint() )
{
	LinearInterpolation < float, double > lininterpol;

	lininterpol.setMapping( 1, 23, 53 );
	TEST_REAL_EQUAL(lininterpol.getScale(), 1)
	TEST_REAL_EQUAL(lininterpol.getInsideReferencePoint(), 23)
	TEST_REAL_EQUAL(lininterpol.getOutsideReferencePoint(), 53)
		
	lininterpol.setMapping( 1, 0, 53 );
	TEST_REAL_EQUAL(lininterpol.supportMin(), 53)
	TEST_REAL_EQUAL(lininterpol.supportMax(), 53)

	lininterpol.setMapping( 1, 500, 53 );
	TEST_REAL_EQUAL(lininterpol.supportMin(), 53-500)
	TEST_REAL_EQUAL(lininterpol.supportMax(), 53-500)

	lininterpol.getData().resize(300);

	lininterpol.setMapping( 10, 0, 1000 );
	TEST_REAL_EQUAL(lininterpol.supportMin(), 990)
	TEST_REAL_EQUAL(lininterpol.supportMax(), 4000)

	lininterpol.setMapping( 10, 200, 1000 );
	TEST_REAL_EQUAL(lininterpol.getScale(), 10)
	TEST_REAL_EQUAL(lininterpol.getInsideReferencePoint(), 200)
	TEST_REAL_EQUAL(lininterpol.getOutsideReferencePoint(), 1000)
	TEST_REAL_EQUAL(lininterpol.getOffset(), -1000)
	TEST_REAL_EQUAL(lininterpol.supportMin(), -1010)
	TEST_REAL_EQUAL(lininterpol.supportMax(), 2000)

}
RESULT
//-----------------------------------------------------------

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
