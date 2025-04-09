// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <iterator>

/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

using namespace OpenMS;

START_TEST(DPosition<D>, "$Id$")

/////////////////////////////////////////////////////////////

std::cout.precision(writtenDigits<>(double()));
std::cerr.precision(writtenDigits<>(double()));

DPosition<10>* d10_ptr = nullptr;
DPosition<10>* d10_nullPointer = nullptr;

START_SECTION((DPosition()))
	d10_ptr = new DPosition<10>;
  TEST_NOT_EQUAL(d10_ptr, d10_nullPointer)
END_SECTION

START_SECTION((~DPosition()))
	delete d10_ptr;
END_SECTION

START_SECTION(void swap(DPosition& rhs) noexcept)
{
  DPosition<3> i(1, 2, 3);
  DPosition<3> j(4, 5, 6);
  i.swap(j);
  TEST_REAL_SIMILAR(i[0], 4)
  TEST_REAL_SIMILAR(i[1], 5)
  TEST_REAL_SIMILAR(i[2], 6)
  TEST_REAL_SIMILAR(j[0], 1)
  TEST_REAL_SIMILAR(j[1], 2)
  TEST_REAL_SIMILAR(j[2], 3)
}
END_SECTION

START_SECTION(DPosition& abs() noexcept)
{
  // a bit of fuzz, just to make sure we call the correct std::abs() function for the appropriate data type
  constexpr const auto weird_negative_int = std::numeric_limits<int64_t>::lowest() + 2; // this value cannot be accurately represented by a double
  constexpr const auto weird_positive_int = -weird_negative_int;
  constexpr const double inaccutate_double(weird_negative_int);
  static_assert(int64_t(inaccutate_double) != weird_negative_int); // make sure its inaccurate
  DPosition<3, Int64> i(weird_negative_int, -5, weird_positive_int); // test if we call the correct abs() function, i.e. the one for int, not double
  i.abs();
  TEST_EQUAL(i[0], weird_positive_int)
  TEST_EQUAL(i[1], 5)
  TEST_EQUAL(i[2], weird_positive_int)
  // test we call abs() for double, not for float
  const auto small_negative_double = -std::numeric_limits<double>::epsilon();
  DPosition<3, double> j(-1.4444, -small_negative_double, small_negative_double);
  j.abs();
  TEST_EQUAL(j[0], 1.4444)
  TEST_EQUAL(j[1], -small_negative_double)  // test equal, not similar!
  TEST_EQUAL(j[2], -small_negative_double)  // test equal, not similar!
}
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

START_SECTION((DPosition(CoordinateType x, CoordinateType y)))
  DPosition<2> p(1, 2);
  TEST_REAL_SIMILAR(p[0], 1)
  TEST_REAL_SIMILAR(p[1], 2)
END_SECTION

START_SECTION(DPosition(CoordinateType x, CoordinateType y, CoordinateType z))
  DPosition<3> p(1, 2, 3);
  TEST_REAL_SIMILAR(p[0], 1)
  TEST_REAL_SIMILAR(p[1], 2)
  TEST_REAL_SIMILAR(p[2], 3)
END_SECTION

START_SECTION((CoordinateType operator*(const DPosition& point) const))
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
	TEST_TRUE(p1 == p2)

	p1[0]=1.234;
	TEST_EQUAL(p1==p2, false)
	p2[0]=1.234;
	TEST_TRUE(p1 == p2)

	p1[1]=1.345;
	TEST_EQUAL(p1==p2, false)
	p2[1]=1.345;
	TEST_TRUE(p1 == p2)

	p1[2]=1.456;
	TEST_EQUAL(p1==p2, false)
	p2[2]=1.456;
	TEST_TRUE(p1 == p2)
END_SECTION

START_SECTION((bool operator!=(const DPosition &point) const))
	DPosition<3> p1,p2;
	TEST_EQUAL(p1!=p2, false)

	p1[0]=1.234;
	TEST_FALSE(p1 == p2)
	p2[0]=1.234;
	TEST_EQUAL(p1!=p2, false)

	p1[1]=1.345;
	TEST_FALSE(p1 == p2)
	p2[1]=1.345;
	TEST_EQUAL(p1!=p2, false)

	p1[2]=1.456;
	TEST_FALSE(p1 == p2)
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
  TEST_FALSE(p1 == p2);
	p2 = -p2;
	TEST_TRUE(p1 == p2);
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

START_SECTION(([EXTRA] Test int DPosition))
{
  DPosition<2,Int> const p00(0,0), p01(0,1), p10(1,0), p11(1,1);
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
}
END_SECTION

START_SECTION(([EXTRA] Test char DPosition))
{
  DPosition<3,char> pa1;
  DPosition<3,char> pb2;
  pa1[0] = 'a';
  pb2 = -pa1;
  TEST_FALSE(pa1 == pb2)
  pb2 = -pb2;
  TEST_TRUE(pa1 == pb2)

  DPosition<1,char> pa('a');
  DPosition<1,char> pb('b');

  TEST_EQUAL(pa < pb, true)
}
END_SECTION

START_SECTION(([EXTRA] Test scalar multiplication))
{
  DPosition<2> p1(3,4);
  DPosition<2> p2 = p1 * 5;
  DPosition<2> p3 = 5 * p1;

  DPosition<2> const pResult(15,20);

  TEST_REAL_SIMILAR(p2[0],pResult[0]);
  TEST_REAL_SIMILAR(p2[1],pResult[1]);

  TEST_REAL_SIMILAR(p3[0],pResult[0]);
  TEST_REAL_SIMILAR(p3[1],pResult[1]);
}
END_SECTION

START_SECTION(([EXTRA] Test scalar division))
{
  DPosition<2> p1(15,20);
  DPosition<2> p2 = p1 / 5;
  DPosition<2> const pResult(3,4);

  TEST_REAL_SIMILAR(p2[0],pResult[0]);
  TEST_REAL_SIMILAR(p2[1],pResult[1]);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

