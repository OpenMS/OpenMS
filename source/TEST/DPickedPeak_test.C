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

#include <OpenMS/KERNEL/DPickedPeak.h>
#include <vector>

///////////////////////////

START_TEST(DPickedPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

DPickedPeak<10>* d10_ptr = 0;
CHECK((DPickedPeak()))
	d10_ptr = new DPickedPeak<10>;
	TEST_NOT_EQUAL(d10_ptr, 0)
RESULT

CHECK((~DPickedPeak()))
	delete d10_ptr;
RESULT

CHECK((const RValueType& getRValue() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
RESULT

CHECK((const AreaType& getArea() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getArea(), 0.0)
RESULT

CHECK((const FullWidthHalfMaxType& getFWHM() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
RESULT

CHECK((const WidthType& getLeftWidthParameter() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
RESULT

CHECK((const WidthType& getRightWidthParameter() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
RESULT

CHECK((const PeakShapeType::Enum& getPeakShape() const))
	const DPickedPeak<10> p;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::UNDEFINED)
RESULT

CHECK((const SignalToNoiseType& getSN() const))
	const DPickedPeak<10> p;
	TEST_REAL_EQUAL(p.getSN(), 0.0)
RESULT

CHECK((const ChargeType& getCharge() const))
	const DPickedPeak<10>	p;
	TEST_EQUAL(p.getCharge(), 0)
RESULT

CHECK((const RValueType& getRValue() const))
	DPickedPeak<10>	p;
	TEST_EQUAL(p.getRValue(), 0)
RESULT

CHECK((RValueType& getRValue()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
	p.setRValue(0.456);
	TEST_REAL_EQUAL(p.getRValue(), 0.456)
	p.setRValue(-0.12345);
	TEST_REAL_EQUAL(p.getRValue(), -0.12345)
	p.setRValue(0.0);
	TEST_REAL_EQUAL(p.getRValue(), 0.0)
RESULT

CHECK((AreaType& getArea()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getArea(), 0.0)
	p.setArea(123.456);
	TEST_REAL_EQUAL(p.getArea(), 123.456)
	p.setArea(-0.12345);
	TEST_REAL_EQUAL(p.getArea(), -0.12345)
	p.setArea(0.0);
	TEST_REAL_EQUAL(p.getArea(), 0.0)
RESULT

CHECK((FullWidthHalfMaxType& getFWHM()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
	p.setFWHM(123.456);
	TEST_REAL_EQUAL(p.getFWHM(), 123.456)
	p.setFWHM(-0.12345);
	TEST_REAL_EQUAL(p.getFWHM(), -0.12345)
	p.setFWHM(0.0);
	TEST_REAL_EQUAL(p.getFWHM(), 0.0)
RESULT

CHECK((WidthType& getLeftWidthParameter()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
	p.setLeftWidthParameter(123.456);
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 123.456)
	p.setLeftWidthParameter(-0.12345);
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), -0.12345)
	p.setLeftWidthParameter(0.0);
	TEST_REAL_EQUAL(p.getLeftWidthParameter(), 0.0)
RESULT

CHECK((WidthType& getRightWidthParameter()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
	p.setRightWidthParameter(123.456);
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 123.456)
	p.setRightWidthParameter(-0.12345);
	TEST_REAL_EQUAL(p.getRightWidthParameter(), -0.12345)
	p.setRightWidthParameter(0.0);
	TEST_REAL_EQUAL(p.getRightWidthParameter(), 0.0)
RESULT

CHECK((PeakShapeType::Enum& getPeakShape()))
	DPickedPeak<3> p;
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::UNDEFINED)
	p.setPeakShape(PeakShapeType::LORENTZ_PEAK);
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::LORENTZ_PEAK)
	p.setPeakShape(PeakShapeType::SECH_PEAK);
	TEST_EQUAL(p.getPeakShape(), PeakShapeType::SECH_PEAK)
RESULT

CHECK((SignalToNoiseType& getSN()))
	DPickedPeak<3> p;
	TEST_REAL_EQUAL(p.getSN(), 0.0)
	p.setSN(123.456);
	TEST_REAL_EQUAL(p.getSN(), 123.456)
	p.setSN(-0.12345);
	TEST_REAL_EQUAL(p.getSN(), -0.12345)
	p.setSN(0.0);
	TEST_REAL_EQUAL(p.getSN(), 0.0)
RESULT

CHECK((ChargeType& getCharge()))
	DPickedPeak<3> p;
	TEST_EQUAL(p.getCharge(), 0.0)
	p.setCharge(12);
	TEST_EQUAL(p.getCharge(), 12)
	p.setCharge(-3);
	TEST_EQUAL(p.getCharge(), -3)
	p.setCharge(0);
	TEST_EQUAL(p.getCharge(), 0)
RESULT


CHECK((DPickedPeak(DPickedPeak const& p)))
	DPickedPeak<3>::PositionType pos,pos2;

	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPickedPeak<3> p;
	p.setRValue(0.523);
	p.setArea(10.234);
	p.setFWHM(23.543);
  p.setSN(12.212);
  p.setLeftWidthParameter(2.35);
  p.setRightWidthParameter(p.getLeftWidthParameter());
  p.setPeakShape(PeakShapeType::LORENTZ_PEAK);
	p.setIntensity(123.456);
	p.setCharge(1234);
	p.setPosition(pos);
	p.setMetaValue("cluster_id",4711);
	
	DPickedPeak<3>::RValueType rv2;
	DPickedPeak<3>::AreaType a2;
	DPickedPeak<3>::FullWidthHalfMaxType fwhm2;
  DPickedPeak<3>::SignalToNoiseType sn2;
  PeakShapeType::Enum type2;
	DPickedPeak<3>::WidthType left_w2,right_w2;
	DPickedPeak<3>::IntensityType i2;

	DPickedPeak<3> copy_of_p(p);
	
	rv2 = copy_of_p.getRValue();
	a2 = copy_of_p.getArea();
	fwhm2 = copy_of_p.getFWHM();
  sn2 = copy_of_p.getSN();
  type2 = copy_of_p.getPeakShape();
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
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
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	

  TEST_EQUAL(type2,PeakShapeType::LORENTZ_PEAK)
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4711));
RESULT

CHECK((DPickedPeak& operator = (const DPickedPeak& rhs)))
DPickedPeak<3>::WidthType w;
DPickedPeak<3>::PositionType pos, pos2;
	w = 13.13;
	pos[0] = 21.21;
	pos[1] = 22.22;
	pos[2] = 23.23;
	DPickedPeak<3> p;
	p.setRValue(0.523);
	p.setArea(10.234);
	p.setFWHM(23.543);
  p.setSN(23.523);
  p.setPeakShape(PeakShapeType::SECH_PEAK);
	p.setIntensity(123.456);
	p.setCharge(1234);
	p.setPosition(pos);
	p.setLeftWidthParameter(w);
  p.setRightWidthParameter(w);
	p.setMetaValue("cluster_id",4712);
	
  DPickedPeak<3>::RValueType rv2;
	DPickedPeak<3>::AreaType a2;
	DPickedPeak<3>::FullWidthHalfMaxType fwhm2;
  DPickedPeak<3>::SignalToNoiseType sn2;
  PeakShapeType::Enum type2;
	DPickedPeak<3>::WidthType left_w2,right_w2;
	DPickedPeak<3>::IntensityType i2;

	DPickedPeak<3> copy_of_p;
	copy_of_p = p;
		
	
	rv2 = copy_of_p.getRValue();
	a2 = copy_of_p.getArea();
	fwhm2 = copy_of_p.getFWHM();
  sn2 = copy_of_p.getSN();
  type2 = copy_of_p.getPeakShape();
	i2 = copy_of_p.getIntensity();
	pos2 = copy_of_p.getPosition();
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
	TEST_REAL_EQUAL(pos2[1], 22.22)
	TEST_REAL_EQUAL(pos2[2], 23.23)	

  TEST_EQUAL(type2,PeakShapeType::SECH_PEAK)
	TEST_EQUAL(p.getMetaValue("cluster_id"),DataValue(4712));
RESULT

CHECK((bool operator == (const DPickedPeak& rhs) const))
	DPickedPeak<1> p1;
	DPickedPeak<1> p2(p1);
	TEST_REAL_EQUAL(p1==p2, true)

	p1.setRValue(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setRValue(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setArea(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setArea(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setFWHM(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setFWHM(5);
	TEST_REAL_EQUAL(p1==p2, true)	

	p1.setSN(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setSN(5);
	TEST_REAL_EQUAL(p1==p2, true)

	p1.setPeakShape(PeakShapeType::SECH_PEAK);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setPeakShape(PeakShapeType::SECH_PEAK);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setCharge(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setCharge(5);
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, false)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1==p2, true)
	
	p1.setLeftWidthParameter(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setLeftWidthParameter(5);
	TEST_REAL_EQUAL(p1==p2, true)		

	p1.setRightWidthParameter(5);
	TEST_REAL_EQUAL(p1==p2, false)
	p2.setRightWidthParameter(5);
	TEST_REAL_EQUAL(p1==p2, true)	
RESULT

CHECK((bool operator != (const DPickedPeak& rhs) const))
	DPickedPeak<1> p1;
	DPickedPeak<1> p2(p1);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setRValue(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setRValue(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setArea(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setArea(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setFWHM(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setFWHM(5);
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.setSN(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setSN(5);
	TEST_REAL_EQUAL(p1!=p2, false)

	p1.setPeakShape(PeakShapeType::SECH_PEAK);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setPeakShape(PeakShapeType::SECH_PEAK);
  TEST_REAL_EQUAL(p1!=p2, false) 

	p1.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setIntensity(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setCharge(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setCharge(5);
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.getPosition()[0]=5;
	TEST_REAL_EQUAL(p1!=p2, false)
	
	p1.setLeftWidthParameter(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setLeftWidthParameter(5);
	TEST_REAL_EQUAL(p1!=p2, false)	

	p1.setRightWidthParameter(5);
	TEST_REAL_EQUAL(p1!=p2, true)
	p2.setRightWidthParameter(5);
	TEST_REAL_EQUAL(p1!=p2, false)	
RESULT

CHECK(([EXTRA]meta info with copy constructor))
	DPickedPeak<1> p;
	p.setMetaValue(2,String("bla"));
 	DPickedPeak<1> p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(([EXTRA]meta info with assignment))
	DPickedPeak<1> p;
	p.setMetaValue(2,String("bla"));
 	DPickedPeak<1> p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
