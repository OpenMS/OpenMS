// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

///////////////////////////

// This one is going to be tested.
#include <OpenMS/ML/INTERPOLATION/LinearInterpolation.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <vector>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/test_config.h>
#include <OpenMS/MATH/MathFunctions.h>

///////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( LinearInterpolation, "$Id$" )

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION([EXTRA] typedefs )
{
	typedef LinearInterpolation < float, double > LIFD;
	LIFD::ValueType     * value = new LIFD::ValueType();
	LIFD::KeyType       * key = new LIFD::KeyType();
	LIFD::ContainerType * container = new LIFD::ContainerType();
	LIFD::ContainerType::value_type * containerValue = new LIFD::ContainerType::value_type();

	LIFD::ValueType     * nullValue = nullptr;
	LIFD::KeyType       * nullKey = nullptr;
	LIFD::ContainerType * nullContainer = nullptr;
	LIFD::ContainerType::value_type * nullContainerValue = nullptr;

  TEST_NOT_EQUAL(value, nullValue)
  TEST_NOT_EQUAL(key, nullKey)
  TEST_NOT_EQUAL(container, nullContainer)
  TEST_NOT_EQUAL(containerValue, nullContainerValue)
  delete value;
  delete key;
  delete container;
  delete containerValue;
  delete nullValue;
  delete nullKey;
  delete nullContainer;
  delete nullContainerValue;
}
END_SECTION

typedef LinearInterpolation < float, double > LIFD;

//-----------------------------------------------------------
// Without these extra parens, check_test will not recognize this test...
START_SECTION((LinearInterpolation(KeyType scale=1., KeyType offset=0.)))
{
	LIFD lifd0;
	LIFD lifd1 ( 1.125 );
	LIFD lifd2 ( 1.125, 3.5 );

	TEST_EQUAL ( lifd0.getScale(), 1 );
	TEST_EQUAL ( lifd0.getOffset(), 0 );

	TEST_EQUAL ( lifd1.getScale(), 1.125 );
	TEST_EQUAL ( lifd1.getOffset(), 0 );

	TEST_EQUAL ( lifd2.getScale(), 1.125 );
	TEST_EQUAL ( lifd2.getOffset(), 3.5 );
}
END_SECTION

LIFD * lifd_nullPointer = nullptr;

START_SECTION(~LinearInterpolation())
{
	LIFD * lifd_ptr = nullptr;
	lifd_ptr = new LIFD;
  TEST_NOT_EQUAL(lifd_ptr,lifd_nullPointer);
	delete lifd_ptr;
}
END_SECTION

START_SECTION(ContainerType const& getData() const )
{
	LIFD lifd;
	std::vector < LIFD::ValueType > v;
	v.push_back(17);
	v.push_back(18.9);
	v.push_back(20.333);
	v.push_back(-.1);
	lifd.setData(v);
	TEST_EQUAL(lifd.getData().size(),v.size());
	for ( UInt i = 0; i < v.size(); ++i )
		TEST_EQUAL(lifd.getData()[i],v[i]);
}
END_SECTION

START_SECTION(template< typename SourceContainer > void setData( SourceContainer const & data ) )
{
	// see above, getData()
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(ContainerType& getData() )
{
	LIFD lifd;
	std::vector < LIFD::ValueType > v;
	v.push_back(17);
	v.push_back(18.9);
	v.push_back(20.333);
	v.push_back(-.1);
	lifd.setData(v);
	LIFD const & lifd_cr = lifd;
	TEST_EQUAL(lifd_cr.getData().size(),v.size());
	for ( UInt i = 0; i < v.size(); ++i )
		TEST_EQUAL(lifd_cr.getData()[i],v[i]);
}
END_SECTION

START_SECTION(bool empty() const )
{
	LIFD lifd;
	TEST_EQUAL(lifd.getData().empty(),true);
	lifd.getData().push_back(3);
	TEST_EQUAL(lifd.getData().empty(),false);
	lifd.getData().push_back(3);
	TEST_EQUAL(lifd.getData().empty(),false);
	lifd.getData().push_back(3);
	TEST_EQUAL(lifd.getData().empty(),false);
	lifd.getData().clear();
	TEST_EQUAL(lifd.getData().empty(),true);
}
END_SECTION

START_SECTION((void setMapping( KeyType const & scale, KeyType const & inside, KeyType const & outside )))
{
	LIFD lifd;
	lifd.setMapping( 13, 23, 53 );
	TEST_REAL_SIMILAR(lifd.getScale(), 13)
	TEST_REAL_SIMILAR(lifd.getInsideReferencePoint(), 23)
	TEST_REAL_SIMILAR(lifd.getOutsideReferencePoint(), 53)
}
END_SECTION

START_SECTION(KeyType const& getScale() const )
{
	LIFD lifd;
	lifd.setMapping( 13, 23, 53 );
	TEST_REAL_SIMILAR(lifd.getScale(), 13.F);
}
END_SECTION

START_SECTION(KeyType const& getInsideReferencePoint() const )
{
	LIFD lifd;
	lifd.setMapping( 13, 23, 53 );
	TEST_REAL_SIMILAR(lifd.getInsideReferencePoint(), 23.F)
}
END_SECTION

START_SECTION(KeyType const& getOutsideReferencePoint() const )
{
	LIFD lifd;
	lifd.setMapping( 13, 23, 53 );
	TEST_REAL_SIMILAR(lifd.getOutsideReferencePoint(), 53.F)
}
END_SECTION

START_SECTION(void setScale( KeyType const & scale ) )
{
	LIFD lifd;
	lifd.setMapping( 13, 23, 53 );
	TEST_REAL_SIMILAR(lifd.getScale(), 13.f);
	lifd.setScale(88.88f);
	TEST_REAL_SIMILAR(lifd.getScale(), 88.88f);
}
END_SECTION

START_SECTION(KeyType const& getOffset() const )
{
	LIFD lifd ( 1.125, 3.5 );
	TEST_EQUAL ( lifd.getOffset(), 3.5f );
}
END_SECTION

START_SECTION(void setOffset( KeyType const & offset ) )
{
	LIFD lifd ( 1.125, 3.5 );
	TEST_EQUAL ( lifd.getOffset(), 3.5f );
	lifd.setOffset(88.88f);
	TEST_EQUAL ( lifd.getOffset(), 88.88f );
}
END_SECTION

START_SECTION((void setMapping( KeyType const & inside_low, KeyType const & outside_low, KeyType const & inside_high, KeyType const & outside_high )))
{
	LIFD lifd;

	lifd.setMapping( 13, 130, 14, 140 );
	TEST_REAL_SIMILAR(lifd.getScale(), 10)
	TEST_REAL_SIMILAR(lifd.getInsideReferencePoint(), 13)
	TEST_REAL_SIMILAR(lifd.getOutsideReferencePoint(), 130)
}
END_SECTION


START_SECTION(LinearInterpolation( LinearInterpolation const & arg ))
{
	LIFD lifd;
	lifd.setMapping( 13, 130, 14, 140 );
	LIFD::ContainerType v;
	v.push_back(17);
	v.push_back(18.9);
	v.push_back(20.333);
	v.push_back(-.1);
	lifd.setData(v);

	LIFD lifd2 = lifd;
	TEST_REAL_SIMILAR(lifd2.getScale(), 10);
	TEST_REAL_SIMILAR(lifd2.getInsideReferencePoint(), 13);
	TEST_REAL_SIMILAR(lifd2.getOutsideReferencePoint(), 130);
	for ( UInt i = 0; i < v.size(); ++i )
		TEST_EQUAL(lifd2.getData()[i],v[i]);
}
END_SECTION

START_SECTION(LinearInterpolation& operator= ( LinearInterpolation const & arg ))
{
	LIFD lifd;
	lifd.setMapping( 13, 130, 14, 140 );
	LIFD::ContainerType v;
	v.push_back(17);
	v.push_back(18.9);
	v.push_back(20.333);
	v.push_back(-.1);
	lifd.setData(v);

	LIFD lifd2;
	lifd2 = lifd;
	TEST_REAL_SIMILAR(lifd2.getScale(), 10);
	TEST_REAL_SIMILAR(lifd2.getInsideReferencePoint(), 13);
	TEST_REAL_SIMILAR(lifd2.getOutsideReferencePoint(), 130);
	for ( UInt i = 0; i < v.size(); ++i )
		TEST_EQUAL(lifd2.getData()[i],v[i]);

}
END_SECTION

START_SECTION(KeyType index2key( KeyType pos ) const )
{
	LIFD lifd(100,3456);

	TEST_REAL_SIMILAR(lifd.index2key(0),3456);
	TEST_REAL_SIMILAR(lifd.index2key(-1),3356);
	TEST_REAL_SIMILAR(lifd.index2key(1),3556);
}
END_SECTION

START_SECTION(KeyType key2index( KeyType pos ) const )
{
	LIFD lifd(100,3456);

	TEST_REAL_SIMILAR(lifd.key2index(3456),0);
	TEST_REAL_SIMILAR(lifd.key2index(3356),-1);
	TEST_REAL_SIMILAR(lifd.key2index(3556),1);

	lifd.setScale(0);
	TEST_REAL_SIMILAR(lifd.key2index(3456),0);
	TEST_REAL_SIMILAR(lifd.key2index(3356),0);
	TEST_REAL_SIMILAR(lifd.key2index(3556),0);
}
END_SECTION

START_SECTION(KeyType supportMin() const )
{
	LIFD lifd ( 1.125, 3.5 );
	TEST_REAL_SIMILAR ( lifd.supportMin(), 3.5 );
	lifd.getData().push_back(11111);
	TEST_REAL_SIMILAR ( lifd.supportMin(), 3.5-1.125 );
	lifd.getData().push_back(99999);
	TEST_REAL_SIMILAR ( lifd.supportMin(), 3.5-1.125 );
	lifd.getData().clear();
	TEST_REAL_SIMILAR ( lifd.supportMin(), 3.5 );
}
END_SECTION

START_SECTION(KeyType supportMax() const )
{
	LIFD lifd ( 1.125, 3.5 );
	TEST_REAL_SIMILAR ( lifd.supportMax(), 3.5 );
	lifd.getData().push_back(11111);
	TEST_REAL_SIMILAR ( lifd.supportMax(), 3.5+1.125 );
	lifd.getData().push_back(99999);
	TEST_REAL_SIMILAR ( lifd.supportMax(), 3.5+2*1.125 );
	lifd.getData().clear();
	TEST_REAL_SIMILAR ( lifd.supportMax(), 3.5 );
}
END_SECTION

//-----------------------------------------------------------
START_SECTION(ValueType value( KeyType arg_pos ) const )
{
	LIFD lifd0;

	double values[] = { 1, 2, 0, 1 };
	int const num_values = sizeof (values) / sizeof (*values);
	std::copy ( values, values + num_values, std::back_inserter ( lifd0.getData() ) );

	TEST_EQUAL ( (int) lifd0.getData().size(), num_values );

	for ( int i = 0; i < num_values; ++i )
	{
		TEST_EQUAL ( lifd0.value ( LIFD::KeyType(i) ), values[i] );
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
		TEST_REAL_SIMILAR ( lifd0.value ( i-2.f ), inter_values[4*i] );
	}

	int const num_inter_values = sizeof (inter_values) / sizeof (*inter_values);
	for ( int i = 0; i < num_inter_values; ++i )
	{
		TEST_REAL_SIMILAR ( lifd0.value ( (i-8.f)/4.f ), inter_values[i] );
	}

	LIFD lifd1 (lifd0);

	float const scale = 1;
	float const offset = 100;
	lifd1.setScale( scale );
	lifd1.setOffset( offset );

	for ( int i = -8; i < num_inter_values-8; ++i )
	{
		float pos = i/4.f;
		TEST_REAL_SIMILAR ( lifd1.key2index(lifd1.index2key(pos)), pos );
	}

	for ( int i = -8; i < num_inter_values-8; ++i )
	{
		float pos = i/4.f;
		TEST_REAL_SIMILAR ( lifd1.value ( pos*scale+offset ), lifd0.value ( pos ) ) ;
	}

	{

		TOLERANCE_ABSOLUTE(0.001);

		LIFD lifd_small;
		lifd_small.getData().resize(5,0);
		lifd_small.setMapping( 0, 0, 5, 5 );

		LIFD lifd_big;
		lifd_big.getData().resize(15,0);
		lifd_big.setMapping( 5, 0, 10, 5 );

		for ( int i = 0; i < 5; ++i )
		{
			lifd_small.getData()[i] = lifd_big.getData()[i+5] = i * 25 + 100;
		}
		STATUS("          " << lifd_small.getData());
		STATUS(lifd_big.getData());

		for ( int i = -50; i <= 100; ++i )
		{
			float pos = i / 10.f;
			STATUS(i);
			TEST_REAL_SIMILAR(lifd_small.value(pos),lifd_big.value(pos));
		}

	}



}
END_SECTION
//-----------------------------------------------------------
START_SECTION(ValueType derivative( KeyType arg_pos ) const )
{
	typedef LinearInterpolation < float, double > LIFD;

	LIFD lifd0;

	double values[] = { 1, 2, 0, 1 };
	int const num_values = sizeof (values) / sizeof (*values);
	std::copy ( values, values + num_values, std::back_inserter ( lifd0.getData() ) );

	TEST_EQUAL ( (int) lifd0.getData().size(), num_values );

	for ( int i = 0; i < num_values; ++i )
	{
		TEST_EQUAL ( lifd0.value ( float(i) ), values[i] );
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
		float key = i/4.f;
		int index = i+8;
		STATUS( "key:" << key << "  index:" << index << '\n')
		TEST_REAL_SIMILAR ( lifd0.derivative ( key ), inter_values[index] );
	}

}
END_SECTION
//-----------------------------------------------------------
START_SECTION((void addValue( KeyType arg_pos, ValueType arg_value ) ))
{

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(2.3,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[2],7);
		TEST_REAL_SIMILAR(lininterpol.getData()[3],3);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(0.3,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[0],7);
		TEST_REAL_SIMILAR(lininterpol.getData()[1],3);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(-0.7,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[0],3);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(-1.7,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[0],0);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(3.3,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i])
				}
		TEST_REAL_SIMILAR(lininterpol.getData()[3],7);
		TEST_REAL_SIMILAR(lininterpol.getData()[4],3);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(4.3,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[4],7);
	}

	{
		LinearInterpolation < double, double > lininterpol;

		lininterpol.getData().resize(5);
		lininterpol.addValue(5.3,10);
		for ( UInt i = 0; i != lininterpol.getData().size(); ++i )
		{
			STATUS(i << ": " << lininterpol.getData()[i]);
		}
		TEST_REAL_SIMILAR(lininterpol.getData()[4],0);
	}

	{

		for ( int i = -50; i <= 100; ++i )
		{
			float pos = i / 10.f;
			STATUS(i);

			LIFD lifd_small;
			lifd_small.getData().resize(5);
			lifd_small.setMapping( 0, 0, 5, 5 );
			lifd_small.addValue( pos, 10 );

      for ( LIFD::ContainerType::iterator iter = lifd_small.getData().begin(); iter != lifd_small.getData().end(); ++ iter ) *iter = Math::round(*iter);
			STATUS("          " << lifd_small.getData());

			LIFD lifd_big;
			lifd_big.getData().resize(15);
			lifd_big.setMapping( 5, 0, 10, 5 );
      lifd_big.addValue( pos, 10 );

      for ( LIFD::ContainerType::iterator iter = lifd_big.getData().begin(); iter != lifd_big.getData().end(); ++ iter ) *iter = Math::round(*iter);
			STATUS(lifd_big.getData());

			std::vector < LIFD::ContainerType::value_type > big_infix ( lifd_big.getData().begin()+5, lifd_big.getData().begin()+10 );
			TEST_EQUAL(lifd_small.getData().size(), big_infix.size())
			ABORT_IF(lifd_small.getData().size() != big_infix.size())

			// test in loop to avoid clang++ compiler error
			LIFD::ContainerType::const_iterator lifd_it = lifd_small.getData().begin();
			std::vector < LIFD::ContainerType::value_type >::const_iterator big_infix_it = big_infix.begin();
      for( ; lifd_it != lifd_small.getData().end() && big_infix_it != big_infix.end() ; ++lifd_it , ++big_infix_it)
			{
				TEST_EQUAL(*lifd_it,*big_infix_it);
			}
		}

	}

}
END_SECTION
//-----------------------------------------------------------

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
