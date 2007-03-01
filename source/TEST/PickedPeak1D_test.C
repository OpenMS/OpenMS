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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/PickedPeak1D.h>

///////////////////////////

START_TEST(PickedPeak1D<D>, "$Id: PickedPeak1D_test.C 1300 2007-01-18 07:27:04Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

PickedPeak1D* d10_ptr = 0;
CHECK((PickedPeak1D()))
	d10_ptr = new PickedPeak1D;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~PickedPeak1D()))
	delete d10_ptr;
RESULT

CHECK((const DoubleReal& getRValue() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
RESULT

CHECK((const DoubleReal& getArea() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getArea(), 0.0)
RESULT

CHECK((const DoubleReal& getFWHM() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
RESULT

CHECK((const DoubleReal& getLeftWidthParameter() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
RESULT

CHECK((const DoubleReal& getRightWidthParameter() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
RESULT

CHECK((const PeakShapeType::Enum& getPeakShape() const))
	const PickedPeak1D p;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::UNDEFINED)
RESULT

CHECK((const DoubleReal& getSN() const))
	const PickedPeak1D p;
	TEST_REAL_EQUAL(p.getSN(), 0.0)
RESULT

CHECK((const SignedInt& getCharge() const))
	const PickedPeak1D	p;
	TEST_EQUAL(p.getCharge(), 0)
RESULT

CHECK((const DoubleReal& getRValue() const))
	PickedPeak1D	p;
	TEST_EQUAL(p.getRValue(), 0)
RESULT

CHECK((DoubleReal& getRValue()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
	p.getRValue() = 0.456;
	TEST_REAL_EQUAL(p.getRValue(), 0.456)
	p.getRValue() = -0.12345;
	TEST_REAL_EQUAL(p.getRValue(), -0.12345)
	p.getRValue() = 0.0;
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
RESULT

CHECK((DoubleReal& getArea()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getArea(), 0.0)
	p.getArea() = 123.456;
	TEST_REAL_EQUAL(p.getArea(), 123.456)
	p.getArea() = -0.12345;
	TEST_REAL_EQUAL(p.getArea(), -0.12345)
	p.getArea() = 0.0;
	TEST_REAL_EQUAL(p.getArea(), 0.0)
RESULT

CHECK((DoubleReal& getFWHM()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
	p.getFWHM() = 123.456;
	TEST_REAL_EQUAL(p.getFWHM(), 123.456)
	p.getFWHM() = -0.12345;
	TEST_REAL_EQUAL(p.getFWHM(), -0.12345)
	p.getFWHM() = 0.0;
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
RESULT

CHECK((DoubleReal& getLeftWidthParameter()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
	p.getLeftWidthParameter() = 123.456;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 123.456)
	p.getLeftWidthParameter() = -0.12345;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), -0.12345)
	p.getLeftWidthParameter() = 0.0;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
RESULT

CHECK((DoubleReal& getRightWidthParameter()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
	p.getRightWidthParameter() = 123.456;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 123.456)
	p.getRightWidthParameter() = -0.12345;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), -0.12345)
	p.getRightWidthParameter() = 0.0;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
RESULT

CHECK((PeakShapeType::Enum& getPeakShape()))
	PickedPeak1D p;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::UNDEFINED)
	p.getPeakShape() = PeakShapeType::LORENTZ_PEAK;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::LORENTZ_PEAK)
	p.getPeakShape() = PeakShapeType::SECH_PEAK;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::SECH_PEAK)
RESULT

CHECK((DoubleReal& getSN()))
	PickedPeak1D p;
	TEST_REAL_EQUAL(p.getSN(), 0.0)
	p.getSN() = 123.456;
	TEST_REAL_EQUAL(p.getSN(), 123.456)
	p.getSN() = -0.12345;
	TEST_REAL_EQUAL(p.getSN(), -0.12345)
	p.getSN() = 0.0;
	TEST_REAL_EQUAL(p.getSN(), 0.0)
RESULT

CHECK((SignedInt& getCharge()))
	PickedPeak1D p;
	TEST_EQUAL(p.getCharge(), 0.0)
	p.getCharge() = 12;
	TEST_EQUAL(p.getCharge(), 12)
	p.getCharge() = -3;
	TEST_EQUAL(p.getCharge(), -3)
	p.getCharge() = 0;
	TEST_EQUAL(p.getCharge(), 0)
RESULT


CHECK((PickedPeak1D(PickedPeak1D const& p)))
	PickedPeak1D::PositionType pos,pos2;

	pos[0] = 21.21;
	PickedPeak1D p;
	p.getRValue() = 0.523;
	p.getArea() = 10.234;
	p.getFWHM() = 23.543;
  p.getSN() = 12.212;
  p.getLeftWidthParameter() = 2.35;
  p.getRightWidthParameter() = p.getLeftWidthParameter();
  p.getPeakShape() = PeakShapeType::LORENTZ_PEAK;
	p.getIntensity() = 123.456;
	p.getCharge() = 1234;
	p.getPos() = pos;
	p.setMetaValue("cluster_id",4711);
	
	DoubleReal rv2;
	DoubleReal a2;
	DoubleReal fwhm2;
  DoubleReal sn2;
  PeakShapeType::Enum type2;
	DoubleReal left_w2,right_w2;
	DoubleReal i2;

	PickedPeak1D copy_of_p(p);
	
	rv2 = copy_of_p.getRValue();
	a2 = copy_of_p.getArea();
	fwhm2 = copy_of_p.getFWHM();
  sn2 = copy_of_p.getSN();
  type2 = copy_of_p.getPeakShape();
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	left_w2 = copy_of_p.getLeftWidthParameter();
  right_w2 = copy_of_p.getRightWidthParameter();

	TEST_REAL_EQUAL(rv2, 0.523)
	TEST_REAL_EQUAL(a2, 10.234)
	TEST_REAL_EQUAL(fwhm2, 23.543)
	TEST_REAL_EQUAL(i2, 123.456)
	TEST_EQUAL(copy_of_p.getCharge(), 1234)
	TEST_REAL_EQUAL(sn2, 12.212)
	TEST_REAL_EQUAL(left_w2,2.35)
  TEST_REAL_EQUAL(right_w2,left_w2)
	TEST_REAL_EQUAL(pos2[0], 21.21)

  TEST_EQUAL(type2,PeakShapeType::LORENTZ_PEAK)
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4711));
RESULT

CHECK((PickedPeak1D& operator = (const PickedPeak1D& rhs)))
  DoubleReal w;
  PickedPeak1D::PositionType pos, pos2;
	w = 13.13;
	pos[0] = 21.21;
	
  PickedPeak1D p;
	p.getRValue() = 0.523;
	p.getArea() = 10.234;
	p.getFWHM() = 23.543;
  p.getSN() = 23.523;
  p.getPeakShape() = PeakShapeType::SECH_PEAK;
	p.getIntensity() = 123.456;
	p.getCharge() = 1234;
	p.getPos() = pos;
	p.getLeftWidthParameter() = w;
  p.getRightWidthParameter() = w;
	p.setMetaValue("cluster_id",4712);
	
  DoubleReal rv2;
	DoubleReal a2;
	DoubleReal fwhm2;
  DoubleReal sn2;
  PeakShapeType::Enum type2;
	DoubleReal left_w2,right_w2;
	DoubleReal i2;

	PickedPeak1D copy_of_p;
	copy_of_p = p;
		
	
	rv2 = copy_of_p.getRValue();
	a2 = copy_of_p.getArea();
	fwhm2 = copy_of_p.getFWHM();
  sn2 = copy_of_p.getSN();
  type2 = copy_of_p.getPeakShape();
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPos();
	left_w2 = copy_of_p.getLeftWidthParameter();
  right_w2 = copy_of_p.getRightWidthParameter();

	TEST_REAL_EQUAL(rv2, 0.523)
	TEST_REAL_EQUAL(a2, 10.234)
	TEST_REAL_EQUAL(fwhm2, 23.543)
	TEST_REAL_EQUAL(i2, 123.456)
	TEST_EQUAL(copy_of_p.getCharge(), 1234)
	TEST_REAL_EQUAL(sn2, 23.523)
	TEST_REAL_EQUAL(left_w2,13.13)
  TEST_REAL_EQUAL(right_w2,left_w2)
	TEST_REAL_EQUAL(pos2[0], 21.21)

  TEST_EQUAL(type2,PeakShapeType::SECH_PEAK)
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4712));
RESULT

CHECK((bool operator == (const PickedPeak1D& rhs) const))
	PickedPeak1D p1;
	PickedPeak1D p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)

	p1.getRValue()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getRValue()=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getArea()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getArea()=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getFWHM()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getFWHM()=5;
	TEST_REAL_EQUAL(p1==p2, true)	

	p1.getSN()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getSN()=5;
	TEST_REAL_EQUAL(p1==p2, true)

	p1.getPeakShape()=PeakShapeType::SECH_PEAK;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPeakShape()=PeakShapeType::SECH_PEAK;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getCharge()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getCharge()=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getLeftWidthParameter()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getLeftWidthParameter()=5;
	TEST_REAL_EQUAL(p1==p2, true)		

	p1.getRightWidthParameter()=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getRightWidthParameter()=5;
	TEST_REAL_EQUAL(p1==p2, true)	
RESULT

CHECK((bool operator != (const PickedPeak1D& rhs) const))
	PickedPeak1D p1;
	PickedPeak1D p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getRValue()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getRValue()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getArea()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getArea()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getFWHM()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getFWHM()=5;
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.getSN()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getSN()=5;
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.getPeakShape()=PeakShapeType::SECH_PEAK;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPeakShape()=PeakShapeType::SECH_PEAK;
  TEST_REAL_EQUAL(p1!=p2, false) 

	p1.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getIntensity()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getCharge()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getCharge()=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPos()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getLeftWidthParameter()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getLeftWidthParameter()=5;
	TEST_REAL_EQUAL(p1!=p2, false)	

	p1.getRightWidthParameter()=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getRightWidthParameter()=5;
	TEST_REAL_EQUAL(p1!=p2, false)	
RESULT

CHECK(([EXTRA]meta info with copy constructor))
	PickedPeak1D p;
	p.setMetaValue(2,std::string("bla"));
 	PickedPeak1D p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(([EXTRA]meta info with assignment))
	PickedPeak1D p;
	p.setMetaValue(2,std::string("bla"));
 	PickedPeak1D p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,std::string("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
