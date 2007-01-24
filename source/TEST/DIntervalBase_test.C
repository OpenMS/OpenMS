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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>

///////////////////////////

START_TEST(DIntervalBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;
	
//1D check
DIntervalBase<1,FloatKernelTraits>* ptr1 = 0;

CHECK((DIntervalBase()))
	ptr1 = new DIntervalBase<1,FloatKernelTraits>;
	TEST_NOT_EQUAL(ptr1, 0)
RESULT

CHECK((~DIntervalBase()))
	delete ptr1;
RESULT

//2D check
DIntervalBase<2>* ptr2 = 0;

CHECK([EXTRA] DIntervalBase())
	ptr2 = new DIntervalBase<2>;
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK([EXTRA] ~DIntervalBase())
	delete ptr2;
RESULT

//misc stuff for testing
typedef DIntervalBase<2> I2;
typedef DIntervalBase<2>::PositionType I2Pos;
const I2 empty;
I2Pos p1;
p1[0]=5.0;
p1[1]=17.5;
I2Pos p2;
p2[0]=65.0;
p2[1]=-57.5;

CHECK((PositionType const& max() const))
  TEST_EQUAL(empty.max(),I2Pos::min_negative);
RESULT

CHECK((PositionType const& min() const))
  TEST_EQUAL(empty.min(),I2Pos::max);
RESULT

CHECK([EXTRA] empty)
  TEST_EQUAL(empty.max(),I2Pos::min_negative);
  TEST_EQUAL(empty.min(),I2Pos::max);
RESULT

CHECK([EXTRA] zero)
  TEST_EQUAL(I2::zero.max(),I2Pos::zero);
  TEST_EQUAL(I2::zero.min(),I2Pos::zero);
RESULT

CHECK((void setMinMax(PositionType const & min, PositionType const & max)))
  I2 tmp(empty);
  tmp.setMinMax(p1,p2);
  TEST_REAL_EQUAL(tmp.min()[0],5.0);
  TEST_REAL_EQUAL(tmp.min()[1],-57.5);
  TEST_REAL_EQUAL(tmp.max()[0],65.0);
  TEST_REAL_EQUAL(tmp.max()[1],17.5);
RESULT

CHECK((void setMin(PositionType const & position)))
  I2 tmp(empty);
  tmp.setMin(p1);
  TEST_EQUAL(tmp.min(),p1);
  TEST_EQUAL(tmp.max(),p1);
  tmp.setMin(p2);
  TEST_REAL_EQUAL(tmp.min()[0],65.0);
  TEST_REAL_EQUAL(tmp.min()[1],-57.5);
  TEST_REAL_EQUAL(tmp.max()[0],65.0);
  TEST_REAL_EQUAL(tmp.max()[1],17.5);
RESULT

CHECK((void setMax(PositionType const & position)))
  I2 tmp(empty);
  tmp.setMax(p1);
  TEST_EQUAL(tmp.min(),p1);
  TEST_EQUAL(tmp.max(),p1);
  tmp.setMax(p2);
  TEST_REAL_EQUAL(tmp.min()[0],5.0);
  TEST_REAL_EQUAL(tmp.min()[1],-57.5);
  TEST_REAL_EQUAL(tmp.max()[0],65.0);
  TEST_REAL_EQUAL(tmp.max()[1],-57.5);
RESULT

CHECK((bool operator==(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp==tmp,true);
	TEST_EQUAL(tmp==empty,true);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp==empty,false);
RESULT

CHECK((bool operator!=(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp!=tmp,false);
	TEST_EQUAL(tmp!=empty,false);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp!=empty,true);
RESULT

CHECK((DIntervalBase(const DIntervalBase& rhs)))
	I2 tmp(p1,p2);
	I2 tmp2(tmp);
	TEST_EQUAL(tmp==tmp2,true);
RESULT

CHECK((DIntervalBase( PositionType const & minimum, PositionType const & maximum )))
	I2 tmp(p1,p2);
	I2 tmp2(tmp.min(), tmp.max());
	TEST_EQUAL(tmp==tmp2,true);
RESULT

CHECK((DIntervalBase& operator=(const DIntervalBase & rhs)))
	I2 tmp(p1,p2);
	I2 tmp2;
	TEST_EQUAL(tmp==tmp2,false);
	tmp2 = tmp;
	TEST_EQUAL(tmp==tmp2,true);
	tmp2 = tmp = empty;
	TEST_EQUAL(tmp==tmp2,true);
	TEST_EQUAL(tmp==empty,true);
RESULT

CHECK((void clear()))
	I2 tmp;
	TEST_EQUAL(tmp==empty,true);
	tmp.setMax(p1);
	TEST_EQUAL(tmp==empty,false);
	tmp.clear();
	TEST_EQUAL(tmp==empty,true);
RESULT

CHECK((PositionType center() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.center());
  TEST_REAL_EQUAL(pos[0],35.0);
  TEST_REAL_EQUAL(pos[1],-20.0);
RESULT

CHECK((PositionType diagonal() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.diagonal());
  TEST_REAL_EQUAL(pos[0],60.0);
  TEST_REAL_EQUAL(pos[1],75.0);
RESULT

CHECK((CoordinateType width() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.width(),60.0)
RESULT

CHECK((CoordinateType height() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.height(),75.0)
RESULT

CHECK((CoordinateType maxX() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.maxX(),65.0)
RESULT

CHECK((CoordinateType maxY() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.maxY(),17.5)
RESULT

CHECK((CoordinateType minX() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.minX(),5.0)
RESULT

CHECK((CoordinateType minY() const))
	I2 tmp(p1,p2);
	TEST_REAL_EQUAL(tmp.minY(),-57.5)
RESULT

CHECK((void setMinX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinX(57.67);
	TEST_REAL_EQUAL(tmp.minX(),57.67)
RESULT

CHECK((void setMaxX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxX(57.67);
	TEST_REAL_EQUAL(tmp.maxX(),57.67)
RESULT

CHECK((void setMinY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinY(57.67);
	TEST_REAL_EQUAL(tmp.minY(),57.67)
RESULT

CHECK((void setMaxY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxY(57.67);
	TEST_REAL_EQUAL(tmp.maxY(),57.67)
RESULT



CHECK((template < D2> void assign(const DIntervalBase< D2 > rhs)))
DIntervalBase<2>::PositionType p1;
p1[0]=5.0;
p1[1]=17.5;
DIntervalBase<2>::PositionType p2;
p2[0]=65.0;
p2[1]=-57.5;
DIntervalBase<2> i2(p1,p2);

DIntervalBase<3> tmp;
tmp.assign(i2);
TEST_REAL_EQUAL(tmp.min()[0],5.0);
TEST_REAL_EQUAL(tmp.min()[1],-57.5);
TEST_REAL_EQUAL(tmp.max()[0],65.0);
TEST_REAL_EQUAL(tmp.max()[1],17.5);

DIntervalBase<1> tmp2;
tmp2.assign(i2);
TEST_REAL_EQUAL(tmp2.min()[0],5.0);
TEST_REAL_EQUAL(tmp2.max()[0],65.0);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


