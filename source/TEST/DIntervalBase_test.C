// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

START_TEST(DIntervalBase, "$Id$")

/////////////////////////////////////////////////////////////
	
//1D check
DIntervalBase<1>* ptr1 = 0;
DIntervalBase<1>* nullPointer1 = 0;

START_SECTION((DIntervalBase()))
	ptr1 = new DIntervalBase<1>;
  TEST_NOT_EQUAL(ptr1, nullPointer1)
END_SECTION

START_SECTION((~DIntervalBase()))
	delete ptr1;
END_SECTION

//2D check
DIntervalBase<2>* ptr2 = 0;
DIntervalBase<2>* nullPointer2 = 0;

START_SECTION([EXTRA] DIntervalBase())
	ptr2 = new DIntervalBase<2>;
  TEST_NOT_EQUAL(ptr2, nullPointer2)
END_SECTION

START_SECTION([EXTRA] ~DIntervalBase())
	delete ptr2;
END_SECTION

//misc stuff for testing
typedef DIntervalBase<2> I2;
typedef DIntervalBase<2>::PositionType I2Pos;
I2Pos p1;
p1[0]=5.0;
p1[1]=17.5;
I2Pos p2;
p2[0]=65.0;
p2[1]=-57.5;

START_SECTION((PositionType const& maxPosition() const))
  TEST_EQUAL(I2::empty.maxPosition()==I2Pos::minNegative(), true);
  TEST_EQUAL(I2::zero.maxPosition()==I2Pos::zero(), true);
END_SECTION

START_SECTION((PositionType const& minPosition() const))
  TEST_EQUAL(I2::empty.minPosition()==I2Pos::maxPositive(), true);
  TEST_EQUAL(I2::zero.minPosition()==I2Pos::zero(), true);
END_SECTION

START_SECTION((void setMinMax(PositionType const & min, PositionType const & max)))
  I2 tmp(I2::empty);
  tmp.setMinMax(p1,p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);
END_SECTION

START_SECTION((void setMin(PositionType const & position)))
  I2 tmp(I2::empty);
  tmp.setMin(p1);
  TEST_EQUAL(tmp.minPosition(),p1);
  TEST_EQUAL(tmp.maxPosition(),p1);
  tmp.setMin(p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);
END_SECTION

START_SECTION((void setMax(PositionType const & position)))
  I2 tmp(I2::empty);
  tmp.setMax(p1);
  TEST_EQUAL(tmp.minPosition(),p1);
  TEST_EQUAL(tmp.maxPosition(),p1);
  tmp.setMax(p2);
  TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
  TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
  TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
  TEST_REAL_SIMILAR(tmp.maxPosition()[1],-57.5);
END_SECTION

START_SECTION((bool operator==(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp==tmp,true);
	TEST_EQUAL(tmp==I2::empty,true);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp==I2::empty,false);
END_SECTION

START_SECTION((bool operator!=(const DIntervalBase &rhs) const ))
	I2 tmp;
	TEST_EQUAL(tmp!=tmp,false);
	TEST_EQUAL(tmp!=I2::empty,false);
	
	tmp.setMax(p1);
	TEST_EQUAL(tmp!=I2::empty,true);
END_SECTION

START_SECTION((DIntervalBase(const DIntervalBase& rhs)))
	I2 tmp(p1,p2);
	I2 tmp2(tmp);
	TEST_EQUAL(tmp==tmp2,true);
END_SECTION

START_SECTION((DIntervalBase( PositionType const & minimum, PositionType const & maximum )))
	I2 tmp(p1,p2);
	I2 tmp2(tmp.minPosition(), tmp.maxPosition());
	TEST_EQUAL(tmp==tmp2,true);
END_SECTION

START_SECTION((DIntervalBase& operator=(const DIntervalBase & rhs)))
	I2 tmp(p1,p2);
	I2 tmp2;
	TEST_EQUAL(tmp==tmp2,false);
	tmp2 = tmp;
	TEST_EQUAL(tmp==tmp2,true);
	tmp2 = tmp = I2::empty;
	TEST_EQUAL(tmp==tmp2,true);
	TEST_EQUAL(tmp==I2::empty,true);
END_SECTION

START_SECTION((void clear()))
	I2 tmp;
	TEST_EQUAL(tmp==I2::empty,true);
	tmp.setMax(p1);
	TEST_EQUAL(tmp==I2::empty,false);
	tmp.clear();
	TEST_EQUAL(tmp==I2::empty,true);
  TEST_EQUAL(tmp.maxPosition()==I2Pos::minNegative(), true);
  TEST_EQUAL(tmp.minPosition()==I2Pos::maxPositive(), true);
END_SECTION

START_SECTION((PositionType center() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.center());
  TEST_REAL_SIMILAR(pos[0],35.0);
  TEST_REAL_SIMILAR(pos[1],-20.0);
END_SECTION

START_SECTION((PositionType diagonal() const))
  I2 tmp(p1,p2);
  I2Pos pos(tmp.diagonal());
  TEST_REAL_SIMILAR(pos[0],60.0);
  TEST_REAL_SIMILAR(pos[1],75.0);
END_SECTION

START_SECTION((CoordinateType width() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.width(),60.0)
END_SECTION

START_SECTION((CoordinateType height() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.height(),75.0)
END_SECTION

START_SECTION((CoordinateType maxX() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.maxX(),65.0)
END_SECTION

START_SECTION((CoordinateType maxY() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.maxY(),17.5)
END_SECTION

START_SECTION((CoordinateType minX() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.minX(),5.0)
END_SECTION

START_SECTION((CoordinateType minY() const))
	I2 tmp(p1,p2);
	TEST_REAL_SIMILAR(tmp.minY(),-57.5)
END_SECTION

START_SECTION((void setMinX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinX(57.67);
	TEST_REAL_SIMILAR(tmp.minX(),57.67)
END_SECTION

START_SECTION((void setMaxX(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxX(57.67);
	TEST_REAL_SIMILAR(tmp.maxX(),57.67)
END_SECTION

START_SECTION((void setMinY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMinY(57.67);
	TEST_REAL_SIMILAR(tmp.minY(),57.67)
END_SECTION

START_SECTION((void setMaxY(CoordinateType const c)))
	I2 tmp(p1,p2);
	tmp.setMaxY(57.67);
	TEST_REAL_SIMILAR(tmp.maxY(),57.67)
END_SECTION



START_SECTION((template <UInt D2> void assign(const DIntervalBase< D2 > rhs)))
DIntervalBase<2>::PositionType p1;
p1[0]=5.0;
p1[1]=17.5;
DIntervalBase<2>::PositionType p2;
p2[0]=65.0;
p2[1]=-57.5;
DIntervalBase<2> i2(p1,p2);

DIntervalBase<3> tmp;
tmp.assign(i2);
TEST_REAL_SIMILAR(tmp.minPosition()[0],5.0);
TEST_REAL_SIMILAR(tmp.minPosition()[1],-57.5);
TEST_REAL_SIMILAR(tmp.maxPosition()[0],65.0);
TEST_REAL_SIMILAR(tmp.maxPosition()[1],17.5);

DIntervalBase<1> tmp2;
tmp2.assign(i2);
TEST_REAL_SIMILAR(tmp2.minPosition()[0],5.0);
TEST_REAL_SIMILAR(tmp2.maxPosition()[0],65.0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


