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
// $Maintainer: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DistanceMatrix, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DistanceMatrix<double>* ptr = 0;
CHECK(DistanceMatrix())
{
	ptr = new DistanceMatrix<double>();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~DistanceMatrix())
{
	delete ptr;
}
RESULT

CHECK((DistanceMatrix(Value se=1)))
{
	DistanceMatrix<double> dm3(1.0);
	TEST_EQUAL(dm3.dimensionsize(), 1)
}
RESULT

DistanceMatrix<double> dm(8,1.0,1.0);

CHECK((DistanceMatrix(SizeType dimensionsize, Value value=Value(), Value se=1)))
{
	TEST_EQUAL(dm.dimensionsize(), 8)
	TEST_EQUAL(dm(6,7),1)
}
RESULT

DistanceMatrix<double> dm2(dm);

CHECK((DistanceMatrix(const DistanceMatrix &source)))
{
	TEST_EQUAL(dm2.dimensionsize(), 8)
	TEST_EQUAL(dm2(2,3),1)
}
RESULT

CHECK((void resize(size_type i) throw (Exception::OutOfRange)))
{
	dm2.resize(5);
	TEST_EQUAL(dm2.dimensionsize(),5)
	TEST_EQUAL(dm2(2,3),1)
}
RESULT

CHECK((DistanceMatrix& operator=(const DistanceMatrix &rhs)))
{
	dm=dm2;
	TEST_EQUAL(dm.dimensionsize(),5)
	TEST_EQUAL(dm(2,3),1)
}
RESULT

CHECK((SizeType dimensionsize() const ))
{
	TEST_EQUAL(dm.dimensionsize(),5)
}
RESULT

CHECK((void setValue(size_type const i, size_type const j, value_type value)))
{
	dm.setValue(0,1,10);
	dm.setValue(0,2,9);
	dm.setValue(0,3,8);
	dm.setValue(0,4,7);
	dm.setValue(1,2,6);
	dm.setValue(1,3,5);
	dm.setValue(1,4,4);
	dm.setValue(2,3,3);
	dm.setValue(2,4,2);
	dm.setValue(3,4,1);
	TEST_EQUAL(dm.getValue(0,1),10)
	//more tested below
}
RESULT

CHECK((const value_type getValue(size_type const i, size_type const j) const ))
{
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(0,2),9)
	TEST_EQUAL(dm.getValue(0,3),8)
	TEST_EQUAL(dm.getValue(0,4),7)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),5)
	TEST_EQUAL(dm.getValue(1,4),4)
	TEST_EQUAL(dm.getValue(2,3),3)
	TEST_EQUAL(dm.getValue(2,4),2)
	TEST_EQUAL(dm.getValue(3,4),1)
}
RESULT

CHECK((value_type getValue(size_type const i, size_type const j)))
{
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(0,2),9)
	TEST_EQUAL(dm.getValue(0,3),8)
	TEST_EQUAL(dm.getValue(0,4),7)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),5)
	TEST_EQUAL(dm.getValue(1,4),4)
	TEST_EQUAL(dm.getValue(2,3),3)
	TEST_EQUAL(dm.getValue(2,4),2)
	TEST_EQUAL(dm.getValue(3,4),1)
}
RESULT

CHECK((void clear()))
{
	dm2.clear();
	TEST_EQUAL(dm2.dimensionsize(),1)
}
RESULT

CHECK((void setValueQuick(size_type const i, size_type const j, value_type value)))
{
	dm.setValueQuick(0,1,1);
	dm.setValueQuick(0,2,2);
	dm.setValueQuick(0,3,3);
	dm.setValueQuick(0,4,4);
	dm.setValueQuick(1,2,5);
	dm.setValueQuick(1,3,6);
	dm.setValueQuick(1,4,7);
	dm.setValueQuick(2,3,8);
	dm.setValueQuick(2,4,9);
	dm.setValueQuick(3,4,10);
	TEST_EQUAL(dm.getValue(0,1),1)
	TEST_EQUAL(dm.getValue(0,2),2)
	TEST_EQUAL(dm.getValue(0,3),3)
	TEST_EQUAL(dm.getValue(0,4),4)
	TEST_EQUAL(dm.getValue(1,2),5)
	TEST_EQUAL(dm.getValue(1,3),6)
	TEST_EQUAL(dm.getValue(1,4),7)
	TEST_EQUAL(dm.getValue(2,3),8)
	TEST_EQUAL(dm.getValue(2,4),9)
	TEST_EQUAL(dm.getValue(3,4),10)
}
RESULT

CHECK((const value_type operator()(size_type const i, size_type const j) const ))
{
	TEST_EQUAL(dm.getValue(0,1),dm(0,1))
	TEST_EQUAL(dm.getValue(0,2),dm(0,2))
	TEST_EQUAL(dm.getValue(0,3),dm(0,3))
	TEST_EQUAL(dm.getValue(0,4),dm(0,4))
	TEST_EQUAL(dm.getValue(1,2),dm(1,2))
	TEST_EQUAL(dm.getValue(1,3),dm(1,3))
	TEST_EQUAL(dm.getValue(1,4),dm(1,4))
	TEST_EQUAL(dm.getValue(2,3),dm(2,3))
	TEST_EQUAL(dm.getValue(2,4),dm(2,4))
	TEST_EQUAL(dm.getValue(3,4),dm(3,4))
}
RESULT

CHECK((value_type operator()(size_type const i, size_type const j)))
{
	TEST_EQUAL(dm.getValue(0,1),dm(0,1))
	TEST_EQUAL(dm.getValue(0,2),dm(0,2))
	TEST_EQUAL(dm.getValue(0,3),dm(0,3))
	TEST_EQUAL(dm.getValue(0,4),dm(0,4))
	TEST_EQUAL(dm.getValue(1,2),dm(1,2))
	TEST_EQUAL(dm.getValue(1,3),dm(1,3))
	TEST_EQUAL(dm.getValue(1,4),dm(1,4))
	TEST_EQUAL(dm.getValue(2,3),dm(2,3))
	TEST_EQUAL(dm.getValue(2,4),dm(2,4))
	TEST_EQUAL(dm.getValue(3,4),dm(3,4))
}
RESULT

CHECK((void reduce(size_type j)))
{
	dm.reduce(2);

	TEST_EQUAL(dm.getValue(0,1),1)
	TEST_EQUAL(dm.getValue(0,2),3)
	TEST_EQUAL(dm.getValue(0,3),4)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),7)
	TEST_EQUAL(dm.getValue(2,3),10)
	TEST_EQUAL(dm.dimensionsize(),4)
}
RESULT

CHECK((std::pair<UInt,UInt> getMinElementCoordinates() const  throw (Exception::IndexOverflow)))
{
	dm.updateMinElement();
	pair<UInt,UInt> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,0)
	TEST_EQUAL(min.second,1)
}
RESULT

CHECK((void updateMinElement() throw (Exception::OutOfRange)))
{
	dm.setValueQuick(2,3,0.5);
	dm.updateMinElement();
	std::pair<UInt,UInt> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,2)
	TEST_EQUAL(min.second,3)	
}
RESULT

CHECK((UInt index(UInt row, UInt col) const ))
{
  NOT_TESTABLE
}
RESULT

CHECK((std::pair<UInt,UInt> const indexPair(UInt index) const ))
{
	NOT_TESTABLE
}
RESULT

CHECK((bool operator==(DistanceMatrix< Value > const &rhs) const  throw (Exception::Precondition)))
{
	dm2=dm;
	TEST_EQUAL((dm==dm2),true)
}
RESULT


CHECK((bool operator<(DistanceMatrix< Value > const &rhs) const  throw (Exception::Precondition)))
{
	dm2.setValue(2,3,10);
	TEST_EQUAL((dm<dm2),true)
}
RESULT

CHECK((template <typename Value> std::ostream & operator<<(std::ostream &os, const DistanceMatrix< Value > &matrix)))
{
  NOT_TESTABLE
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



