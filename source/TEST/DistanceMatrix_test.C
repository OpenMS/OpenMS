// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
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
START_SECTION(DistanceMatrix())
{
	ptr = new DistanceMatrix<double>();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DistanceMatrix())
{
	delete ptr;
}
END_SECTION

DistanceMatrix<double> dm(8,1.0);

START_SECTION((DistanceMatrix(SizeType dimensionsize, Value value=Value())))
{
	TEST_EQUAL(dm.dimensionsize(), 8)
	TEST_EQUAL(dm(6,7),1)
}
END_SECTION

DistanceMatrix<double> dm2(dm);

START_SECTION((DistanceMatrix(const DistanceMatrix &source)))
{
	TEST_EQUAL(dm2.dimensionsize(), 8)
	TEST_EQUAL(dm2(2,3),1)
}
END_SECTION

START_SECTION(void resize(SizeType dimensionsize, Value value=Value()))
{
	dm2.resize(5);
	TEST_EQUAL(dm2.dimensionsize(),5)
}
END_SECTION

START_SECTION((SizeType dimensionsize() const))
{
	TEST_EQUAL(dm2.dimensionsize(),5)
}
END_SECTION

START_SECTION((void setValue(SizeType i, SizeType j, ValueType value)))
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
	dm.setValue(3,4,0.5);
	TEST_EQUAL(dm.getValue(0,1),10)
	TEST_EQUAL(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),0.5)
	dm.setValue(3,4,1);
	TEST_EQUAL(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),1.0)
	//more tested below
}
END_SECTION

START_SECTION((const ValueType getValue(SizeType i, SizeType j) const))
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
END_SECTION

START_SECTION((ValueType getValue(SizeType i, SizeType j)))
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
END_SECTION

START_SECTION((void clear()))
{
	dm2.clear();
	TEST_EQUAL(dm2.dimensionsize(),0)
}
END_SECTION

START_SECTION((void setValueQuick(SizeType i, SizeType j, ValueType value)))
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
END_SECTION

START_SECTION((const ValueType operator()(SizeType i, SizeType j) const))
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
END_SECTION

START_SECTION((ValueType operator()(SizeType i, SizeType j)))
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
END_SECTION

START_SECTION((void reduce(SizeType j)))
{
	dm.reduce(2);
	TEST_EQUAL(dm.getValue(0,1),1)
	TEST_EQUAL(dm.getValue(0,2),3)
	TEST_EQUAL(dm.getValue(0,3),4)
	TEST_EQUAL(dm.getValue(1,2),6)
	TEST_EQUAL(dm.getValue(1,3),7)
	TEST_EQUAL(dm.getValue(2,3),10)
	TEST_EQUAL(dm.dimensionsize(),7)
}
END_SECTION

START_SECTION((std::pair<SizeType,SizeType> getMinElementCoordinates() const))
{
	dm.updateMinElement();
	pair<Size,Size> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,1)
	TEST_EQUAL(min.second,0)
}
END_SECTION

START_SECTION((void updateMinElement()))
{
	dm.setValueQuick(2,3,0.5);
	dm.updateMinElement();
	std::pair<Size,Size> min = dm.getMinElementCoordinates();
	TEST_EQUAL(min.first,3)
	TEST_EQUAL(min.second,2)
}
END_SECTION

DistanceMatrix<double> dm3(dm);

START_SECTION(bool operator==(DistanceMatrix< ValueType > const &rhs) const)
{
	TEST_EQUAL((dm==dm3),true)
}
END_SECTION


START_SECTION((template <typename Value> std::ostream & operator<<(std::ostream &os, const DistanceMatrix< Value > &matrix)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



