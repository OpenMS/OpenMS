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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <iterator>

/////////////////////////////////////////////////////////////

using namespace OpenMS;

START_TEST(DPosition<D>, "$Id$")

/////////////////////////////////////////////////////////////

std::cout.precision(writtenDigits<>(DoubleReal()));
std::cerr.precision(writtenDigits<>(DoubleReal()));

DPosition<10>* d10_ptr = 0;
DPosition<10>* d10_nullPointer = 0;

START_SECTION((DPosition()))
	d10_ptr = new DPosition<10>;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~DPosition()))
	delete d10_ptr;
END_SECTION

START_SECTION((CoordinateType operator[](Size index) const))
  const DPosition<3> i;
  TEST_EQUAL(i[0], 0.0)
  TEST_EQUAL(i[1], 0.0)
  TEST_EQUAL(i[2], 0.0)
  TEST_PRECONDITION_VIOLATED(i[3])
END_SECTION

START_SECTION((CoordinateType& operator[](Size index)))
  DPosition<3> i;
  const DPosition<3>& c_i(i);
  i[0] = 1.0;
  TEST_REAL_SIMILAR(c_i[0], 1.0)
  TEST_EQUAL(c_i[1], 0.0)
  TEST_EQUAL(c_i[2], 0.0)
  i[1] = 2.0;
  TEST_REAL_SIMILAR(c_i[0], 1.0)
  TEST_REAL_SIMILAR(c_i[1], 2.0)
  TEST_EQUAL(c_i[2], 0.0)
  i[2] = 3.0;
  TEST_REAL_SIMILAR(c_i[0], 1.0)
  TEST_REAL_SIMILAR(c_i[1], 2.0)
  TEST_REAL_SIMILAR(c_i[2], 3.0)
  TEST_PRECONDITION_VIOLATED(i[3] = 4.0)
END_SECTION

START_SECTION((DPosition(const DPosition& pos)))
  DPosition<3> p;
  p[0] = 12.3;
  p[1] = 23.4;
  p[2] = 34.5;
  DPosition<3> copy_of_p(p);
  TEST_EQUAL(copy_of_p[0], p[0])
  TEST_EQUAL(copy_of_p[1], p[1])
  TEST_EQUAL(copy_of_p[2], p[2])
	TEST_EQUAL(copy_of_p.size(), p.size())
END_SECTION

START_SECTION((DPosition& operator=(const DPosition &source)))
  DPosition<3> p;
  p[0] = 12.3;
  p[1] = 23.4;
  p[2] = 34.5;
  DPosition<3> copy_of_p;
  copy_of_p = p;
  TEST_EQUAL(copy_of_p[0], p[0])
  TEST_EQUAL(copy_of_p[1], p[1])
  TEST_EQUAL(copy_of_p[2], p[2])
	TEST_EQUAL(copy_of_p.size(), p.size())
END_SECTION

START_SECTION((DPosition(CoordinateType x)))
  DPosition<3> p(12.34);
  TEST_REAL_SIMILAR(p[0], 12.34)
  TEST_REAL_SIMILAR(p[1], 12.34)
  TEST_REAL_SIMILAR(p[2], 12.34)
END_SECTION

START_SECTION((CoordinateType operator *(const DPosition &point) const))
	DPosition<3> i;
	i[0] = 2.0;
	i[1] = 3.0;
	i[2] = 4.0;
	DPosition<3> j;
	j[0] = 3.0;
	j[1] = 4.0;
	j[2] = 5.0;
	TEST_REAL_SIMILAR(i * j, 6.0 + 12.0 + 20.0)
END_SECTION

DPosition<10> i;
i[0] = 1.0;
i[1] = 2.0;
i[2] = 3.0;
i[3] = 4.0;
i[4] = 5.0;
i[5] = 6.0;
i[6] = 7.0;
i[7] = 8.0;
i[8] = 9.0;
i[9] = 10.0;
START_SECTION((ConstIterator begin() const))
  const DPosition<10>& c_i(i);
  TEST_EQUAL(*c_i.begin(), 1.0)
END_SECTION

START_SECTION((ConstIterator end() const))
  const DPosition<10>& c_i(i);
  TEST_NOT_EQUAL(c_i.end(), c_i.begin())
  std::vector<double> v;
  std::copy(c_i.begin(), c_i.end(), std::back_inserter<std::vector<double> >(v));
  TEST_EQUAL(v.size(), 10)
  ABORT_IF(v.size() != 10)
  TEST_REAL_SIMILAR(v[0], 1.0)
  TEST_REAL_SIMILAR(v[1], 2.0)
  TEST_REAL_SIMILAR(v[2], 3.0)
  TEST_REAL_SIMILAR(v[3], 4.0)
  TEST_REAL_SIMILAR(v[4], 5.0)
  TEST_REAL_SIMILAR(v[5], 6.0)
  TEST_REAL_SIMILAR(v[6], 7.0)
  TEST_REAL_SIMILAR(v[7], 8.0)
  TEST_REAL_SIMILAR(v[8], 9.0)
  TEST_REAL_SIMILAR(v[9], 10.0)
END_SECTION

START_SECTION((Iterator begin()))
  TEST_EQUAL(*i.begin(), 1.0)
  TEST_EQUAL(i.begin(), &(i[0]))
  *i.begin() = 11.0;
  TEST_EQUAL(i[0], 11.0)
END_SECTION

START_SECTION((Iterator end()))
  const DPosition<10>& c_i(i);
  TEST_NOT_EQUAL(c_i.end(), c_i.begin())
  std::vector<double> v;
  std::copy(c_i.begin(), c_i.end(), std::back_inserter<std::vector<double> >(v));
  TEST_EQUAL(v.size(), 10)
  ABORT_IF(v.size() != 10)
  TEST_REAL_SIMILAR(v[0], 11.0)
  TEST_REAL_SIMILAR(v[1], 2.0)
  TEST_REAL_SIMILAR(v[2], 3.0)
  TEST_REAL_SIMILAR(v[3], 4.0)
  TEST_REAL_SIMILAR(v[4], 5.0)
  TEST_REAL_SIMILAR(v[5], 6.0)
  TEST_REAL_SIMILAR(v[6], 7.0)
  TEST_REAL_SIMILAR(v[7], 8.0)
  TEST_REAL_SIMILAR(v[8], 9.0)
  TEST_REAL_SIMILAR(v[9], 10.0)
END_SECTION

START_SECTION((static Size size()))
	TEST_EQUAL(DPosition<777>::size(), 777)
	DPosition<3> p3;
	TEST_EQUAL(p3.size(), 3)
	DPosition<1> p1;
	TEST_EQUAL(p1.size(), 1)
	DPosition<123> p123;
	TEST_EQUAL(p123.size(), 123)
END_SECTION

START_SECTION((void clear()))
	DPosition<3> p;
	p[0] = 1.2;
	p[1] = 2.3;
	p[2] = 3.4;
	TEST_REAL_SIMILAR(p[0], 1.2)
	TEST_REAL_SIMILAR(p[1], 2.3)
	TEST_REAL_SIMILAR(p[2], 3.4)
	p.clear();
	TEST_REAL_SIMILAR(p[0], 0.0)
	TEST_REAL_SIMILAR(p[1], 0.0)
	TEST_REAL_SIMILAR(p[2], 0.0)
END_SECTION

START_SECTION((bool operator==(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1==p2, true)

	p1[0]=1.234;
	TEST_EQUAL(p1==p2, false)
	p2[0]=1.234;
	TEST_EQUAL(p1==p2, true)

	p1[1]=1.345;
	TEST_EQUAL(p1==p2, false)
	p2[1]=1.345;
	TEST_EQUAL(p1==p2, true)

	p1[2]=1.456;
	TEST_EQUAL(p1==p2, false)
	p2[2]=1.456;
	TEST_EQUAL(p1==p2, true)
END_SECTION

START_SECTION((bool operator!=(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1!=p2, false)

	p1[0]=1.234;
	TEST_EQUAL(p1!=p2, true)
	p2[0]=1.234;
	TEST_EQUAL(p1!=p2, false)

	p1[1]=1.345;
	TEST_EQUAL(p1!=p2, true)
	p2[1]=1.345;
	TEST_EQUAL(p1!=p2, false)

	p1[2]=1.456;
	TEST_EQUAL(p1!=p2, true)
	p2[2]=1.456;
	TEST_EQUAL(p1!=p2, false)
END_SECTION

START_SECTION((bool operator<(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1<p2, false)

	p1[0]=p2[0]-0.1;
	TEST_EQUAL(p1<p2, true)
	p2[0]=p1[0]-0.1;
	TEST_EQUAL(p1<p2, false)
	p2[0]=p1[0];

	p1[1]=p2[1]-0.1;
	TEST_EQUAL(p1<p2, true)
	p2[1]=p1[1]-0.1;
	TEST_EQUAL(p1<p2, false)
	p2[1]=p1[1];

	p1[2]=p2[2]-0.1;
	TEST_EQUAL(p1<p2, true)
	p2[2]=p1[2]-0.1;
	TEST_EQUAL(p1<p2, false)
	p2[2]=p1[2];
END_SECTION

START_SECTION((bool operator>(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1>p2, false)

	p1[0]=p2[0]-0.1;
	TEST_EQUAL(p1>p2, false)
	p2[0]=p1[0]-0.1;
	TEST_EQUAL(p1>p2, true)
	p2[0]=p1[0];
END_SECTION

START_SECTION((bool operator>=(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1>=p2, true)

	p1[0]=p2[0]-0.1;
	TEST_EQUAL(p1>=p2, false)
	p2[0]=p1[0]-0.1;
	TEST_EQUAL(p1>=p2, true)
	p2[0]=p1[0];
END_SECTION

START_SECTION((bool operator<=(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1<=p2, true)

	p1[0]=p2[0]-0.1;
	TEST_EQUAL(p1<=p2, true)
	p2[0]=p1[0]-0.1;
	TEST_EQUAL(p1<=p2, false)
END_SECTION

START_SECTION((DPosition operator-() const))
  DPosition<3> p1, p2;
  p1[0] = 5.0;
	p2 = -p1;
  TEST_EQUAL(p1!=p2, true);
	p2 = -p2;
	TEST_EQUAL(p1==p2, true);
END_SECTION

START_SECTION((DPosition operator-(const DPosition &point) const))
  DPosition<3> p1, p2, p3;
  p1[0] = 1.234;
  p1[1] = 2.234;
  p1[2] = 3.234;
  p2[0] = 0.234;
  p2[1] = 0.234;
  p2[2] = 0.234;
  p3[0] = 1.0;
  p3[1] = 2.0;
  p3[2] = 3.0;
  TEST_REAL_SIMILAR(p1[0] - p2[0], p3[0]);
  TEST_REAL_SIMILAR(p1[1] - p2[1], p3[1]);
  TEST_REAL_SIMILAR(p2[0] - p1[0], -p3[0]);
  TEST_REAL_SIMILAR(p2[1] - p1[1], -p3[1]);
END_SECTION

START_SECTION((DPosition operator+(const DPosition &point) const))
  DPosition<3> p1, p2, p3;
  p1[0] = -1.0;
  p1[1] = -2.0;
  p1[2] = -3.0;
  p2[0] = 1.0;
  p2[1] = 2.0;
  p2[2] = 3.0;
  TEST_EQUAL((p1 + p2) ==  p3, true);
END_SECTION

START_SECTION((DPosition(CoordinateType x, CoordinateType y)))
	DPosition<2> p1(11.0f,12.1f);
	TEST_REAL_SIMILAR(p1[0],11.0f);
	TEST_REAL_SIMILAR(p1[1],12.1f);
  DPosition<2> p(12.34,56.78);
  TEST_REAL_SIMILAR(p[0], 12.34)
  TEST_REAL_SIMILAR(p[1], 56.78)
END_SECTION

START_SECTION((CoordinateType getX() const))
	DPosition<2> p1(11.0f,12.1f);
	TEST_REAL_SIMILAR(p1.getX(),11.0f);
END_SECTION

START_SECTION((CoordinateType getY() const))
	DPosition<2> p1(11.0f,12.1f);
	TEST_REAL_SIMILAR(p1.getY(),12.1f);
END_SECTION

START_SECTION((void setX(CoordinateType c)))
	DPosition<2> p1(11.0f,12.1f);
	p1.setX(5.0f);
	TEST_REAL_SIMILAR(p1[0],5.0f);
	TEST_REAL_SIMILAR(p1[1],12.1f);
END_SECTION

START_SECTION((void setY(CoordinateType c)))
	DPosition<2> p1(11.0f,12.1f);
	p1.setY(5.0f);
	TEST_REAL_SIMILAR(p1[0],11.0f);
	TEST_REAL_SIMILAR(p1[1],5.0f);
END_SECTION

START_SECTION((DPosition& operator *=(CoordinateTypescalar)))
	DPosition<2> p1(3,4);
  p1 *= 5;
  DPosition<2> const p2(15,20);
  TEST_REAL_SIMILAR(p1[0],p2[0]);
  TEST_REAL_SIMILAR(p1[1],p2[1]);
END_SECTION

START_SECTION((DPosition& operator+=(const DPosition &point)))
	DPosition<2> p1(3,4);
  DPosition<2> const p2(15,20);
  p1 += p2;
  DPosition<2> const p3(18,24);
  TEST_REAL_SIMILAR(p1[0],p3[0]);
  TEST_REAL_SIMILAR(p1[1],p3[1]);
END_SECTION

START_SECTION((DPosition& operator-=(const DPosition &point)))
	DPosition<2> p1(3,4);
  DPosition<2> const p2(18,24);
  p1 -= p2;
  DPosition<2> const p3(-15,-20);
  TEST_REAL_SIMILAR(p1[0],p3[0]);
  TEST_REAL_SIMILAR(p1[1],p3[1]);
END_SECTION

START_SECTION((DPosition& operator/=(CoordinateTypescalar)))
  DPosition<2> p1(15,20);
  p1 /= 5;
	DPosition<2> const p2(3,4);
  TEST_REAL_SIMILAR(p1[0],p2[0]);
  TEST_REAL_SIMILAR(p1[1],p2[1]);
END_SECTION

START_SECTION((bool spatiallyGreaterEqual(const DPosition &point) const))
	DPosition<2> const p00(0,0), p01(0,1), p10(1,0), p11(1,1);
	TEST_EQUAL(p00.spatiallyGreaterEqual(p00), true )
	TEST_EQUAL(p00.spatiallyGreaterEqual(p01), false)
	TEST_EQUAL(p00.spatiallyGreaterEqual(p10), false)
	TEST_EQUAL(p00.spatiallyGreaterEqual(p11), false)

	TEST_EQUAL(p01.spatiallyGreaterEqual(p00), true )
	TEST_EQUAL(p01.spatiallyGreaterEqual(p01), true )
	TEST_EQUAL(p01.spatiallyGreaterEqual(p10), false)
	TEST_EQUAL(p01.spatiallyGreaterEqual(p11), false)

	TEST_EQUAL(p10.spatiallyGreaterEqual(p00), true )
	TEST_EQUAL(p10.spatiallyGreaterEqual(p01), false)
	TEST_EQUAL(p10.spatiallyGreaterEqual(p10), true )
	TEST_EQUAL(p10.spatiallyGreaterEqual(p11), false)

	TEST_EQUAL(p11.spatiallyGreaterEqual(p00), true )
	TEST_EQUAL(p11.spatiallyGreaterEqual(p01), true )
	TEST_EQUAL(p11.spatiallyGreaterEqual(p10), true )
	TEST_EQUAL(p11.spatiallyGreaterEqual(p11), true )
END_SECTION

START_SECTION((bool spatiallyLessEqual(const DPosition &point) const))
	DPosition<2> const p00(0,0), p01(0,1), p10(1,0), p11(1,1);
	TEST_EQUAL(p00.spatiallyLessEqual(p00), true )
	TEST_EQUAL(p00.spatiallyLessEqual(p01), true )
	TEST_EQUAL(p00.spatiallyLessEqual(p10), true )
	TEST_EQUAL(p00.spatiallyLessEqual(p11), true )

	TEST_EQUAL(p01.spatiallyLessEqual(p00), false)
	TEST_EQUAL(p01.spatiallyLessEqual(p01), true )
	TEST_EQUAL(p01.spatiallyLessEqual(p10), false)
	TEST_EQUAL(p01.spatiallyLessEqual(p11), true )

	TEST_EQUAL(p10.spatiallyLessEqual(p00), false)
	TEST_EQUAL(p10.spatiallyLessEqual(p01), false)
	TEST_EQUAL(p10.spatiallyLessEqual(p10), true )
	TEST_EQUAL(p10.spatiallyLessEqual(p11), true )

	TEST_EQUAL(p11.spatiallyLessEqual(p00), false)
	TEST_EQUAL(p11.spatiallyLessEqual(p01), false)
	TEST_EQUAL(p11.spatiallyLessEqual(p10), false)
	TEST_EQUAL(p11.spatiallyLessEqual(p11), true )
END_SECTION

START_SECTION((static const DPosition zero()))
  typedef DPosition<1> D1;
  TEST_EQUAL(D1::zero()[0],0);
END_SECTION

START_SECTION((static const DPosition minPositive()))
  typedef DPosition<1> D1;
  TEST_EQUAL(D1::minPositive()[0],std::numeric_limits<D1::CoordinateType>::min());
END_SECTION

START_SECTION((static const DPosition minNegative()))
  typedef DPosition<1> D1;
  TEST_EQUAL(D1::minNegative()[0],-std::numeric_limits<D1::CoordinateType>::max());
END_SECTION

START_SECTION((static const DPosition maxPositive()))
  typedef DPosition<1> D1;
  TEST_EQUAL(D1::maxPositive()[0],std::numeric_limits<D1::CoordinateType>::max());
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
