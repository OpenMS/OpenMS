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
// $Maintainer: Ole Schulz-Trieglaff  $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>


///////////////////////////

START_TEST(BaseFeaFiTraits, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeaFiTraits* ptr = 0;
CHECK(FeaFiTraits())
	ptr = new FeaFiTraits();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~FeaFiTraits())
	delete ptr;
RESULT

//create dummy MSExperiment
MSExperiment<> exp;
exp.resize(2);
exp[0].setMSLevel(1);
exp[0].setRetentionTime(1.1);
exp[1].setMSLevel(1);
exp[1].setRetentionTime(2.2);

DPeak<1> p;

p.setPos(500.0);
p.setIntensity(501.0);
exp[0].push_back(p);
p.setPos(700.0);
p.setIntensity(701.0);
exp[0].push_back(p);
p.setPos(900.0);
p.setIntensity(901.0);
exp[0].push_back(p);

p.setPos(600.0);
p.setIntensity(601.0);
exp[1].push_back(p);
p.setPos(1000.0);
p.setIntensity(1001.0);
exp[1].push_back(p);


CHECK(inline const MapType& getData() const)
	FeaFiTraits t;
	TEST_EQUAL(t.getData()==FeaFiTraits::MapType(),true)
RESULT

CHECK(template <class SpectrumIteratorType> void setData(const SpectrumIteratorType& begin, const SpectrumIteratorType& end, UnsignedInt buffer_size))
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_EQUAL(t.getData().getSize(),5)
	TEST_EQUAL(t.getData().size(),2)
	TEST_EQUAL(t.getData()[0].size(),3)
	TEST_EQUAL(t.getData()[1].size(),2)
	TEST_REAL_EQUAL(t.getData()[0][0].getPos(),500.0)
	TEST_REAL_EQUAL(t.getData()[0][0].getIntensity(),501.0)
	TEST_REAL_EQUAL(t.getData()[0][1].getPos(),700.0)
	TEST_REAL_EQUAL(t.getData()[0][1].getIntensity(),701.0)
	TEST_REAL_EQUAL(t.getData()[0][2].getPos(),900.0)
	TEST_REAL_EQUAL(t.getData()[0][2].getIntensity(),901.0)
	TEST_REAL_EQUAL(t.getData()[1][0].getPos(),600.0)
	TEST_REAL_EQUAL(t.getData()[1][0].getIntensity(),601.0)
	TEST_REAL_EQUAL(t.getData()[1][1].getPos(),1000.0)
	TEST_REAL_EQUAL(t.getData()[1][1].getIntensity(),1001.0)
RESULT

CHECK(inline const Flag& getPeakFlag(const IDX& index) const)
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_EQUAL(t.getPeakFlag(make_pair(0,0)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(0,1)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(0,2)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(1,0)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(1,1)), FeaFiTraits::UNUSED)
RESULT

CHECK(inline Flag& getPeakFlag(const IDX& index))
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	t.getPeakFlag(make_pair(0,0)) = FeaFiTraits::SEED;
	TEST_EQUAL(t.getPeakFlag(make_pair(0,0)), FeaFiTraits::SEED)
	TEST_EQUAL(t.getPeakFlag(make_pair(0,1)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(0,2)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(1,0)), FeaFiTraits::UNUSED)
	TEST_EQUAL(t.getPeakFlag(make_pair(1,1)), FeaFiTraits::UNUSED)
RESULT

CHECK(inline const IntensityType& getPeakIntensity(const IDX& index) const)
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_REAL_EQUAL(t.getPeakIntensity(make_pair(0,0)), 501.0)
	TEST_REAL_EQUAL(t.getPeakIntensity(make_pair(0,1)), 701.0)
	TEST_REAL_EQUAL(t.getPeakIntensity(make_pair(0,2)), 901.0)
	TEST_REAL_EQUAL(t.getPeakIntensity(make_pair(1,0)), 601.0)
	TEST_REAL_EQUAL(t.getPeakIntensity(make_pair(1,1)), 1001.0)
RESULT

CHECK(inline const CoordinateType& getPeakMz(const IDX& index) const)
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_REAL_EQUAL(t.getPeakMz(make_pair(0,0)), 500.0)
	TEST_REAL_EQUAL(t.getPeakMz(make_pair(0,1)), 700.0)
	TEST_REAL_EQUAL(t.getPeakMz(make_pair(0,2)), 900.0)
	TEST_REAL_EQUAL(t.getPeakMz(make_pair(1,0)), 600.0)
	TEST_REAL_EQUAL(t.getPeakMz(make_pair(1,1)), 1000.0)
RESULT

CHECK(inline const CoordinateType& getPeakRt(const IDX& index) const)
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_REAL_EQUAL(t.getPeakRt(make_pair(0,0)), 1.1)
	TEST_REAL_EQUAL(t.getPeakRt(make_pair(0,1)), 1.1)
	TEST_REAL_EQUAL(t.getPeakRt(make_pair(0,2)), 1.1)
	TEST_REAL_EQUAL(t.getPeakRt(make_pair(1,0)), 2.2)
	TEST_REAL_EQUAL(t.getPeakRt(make_pair(1,1)), 2.2)
RESULT

CHECK(inline PositionType2D getPeakPos(const IDX& index) const)
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	TEST_EQUAL(t.getPeakPos(make_pair(0,0)),DPosition<2>(1.1,500.0))
	TEST_EQUAL(t.getPeakPos(make_pair(0,1)),DPosition<2>(1.1,700.0))
	TEST_EQUAL(t.getPeakPos(make_pair(0,2)),DPosition<2>(1.1,900.0))
	TEST_EQUAL(t.getPeakPos(make_pair(1,0)),DPosition<2>(2.2,600.0))
	TEST_EQUAL(t.getPeakPos(make_pair(1,1)),DPosition<2>(2.2,1000.0))
RESULT

CHECK(inline void getNextMz(IDX& index) const throw (NoSuccessor))
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	
	//scan one
	FeaFiTraits::IDX i = make_pair(0,0);
	t.getNextMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	t.getNextMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextMz(i));
	
	//scan two
	i = make_pair(1,0);
	t.getNextMz(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextMz(i));

	//test for corrupt index
#ifdef OPENMS_DEBUG
	i = make_pair(5,0);
	TEST_EXCEPTION(Exception::Precondition, t.getNextMz(i));
	i = make_pair(1,5);
	TEST_EXCEPTION(Exception::Precondition, t.getNextMz(i));
#endif
RESULT

CHECK(inline void getPrevMz(IDX& index) const throw (NoSuccessor))
	FeaFiTraits t;
	t.setData(exp.begin(),exp.end(),2);
	//scan one
	FeaFiTraits::IDX i = make_pair(0,2);
	t.getPrevMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	t.getPrevMz(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevMz(i));
	//scan two
	i = make_pair(1,1);
	t.getPrevMz(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevMz(i));

	//test for corrupt index
#ifdef OPENMS_DEBUG
	i = make_pair(5,0);
	TEST_EXCEPTION(Exception::Precondition, t.getPrevMz(i));
	i = make_pair(1,5);
	TEST_EXCEPTION(Exception::Precondition, t.getPrevMz(i));
#endif
RESULT

CHECK(void getNextRt(IDX& index) throw (NoSuccessor))
	FeaFiTraits t;
	MSExperiment<> exp2 = exp;
	exp2.resize(3);
	exp2[2].resize(1);
	exp2[2][0].setPos(800.0);
	exp2[0].resize(5);
	exp2[0][2].setPos(799.0);
	exp2[0][3].setPos(801.0);
	exp2[0][4].setPos(900.0);
	
	t.setData(exp2.begin(),exp2.end(),2);
	FeaFiTraits::IDX i;
	
	//peak one
	i = make_pair(0,0);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextRt(i));

	//peak two
	i = make_pair(0,1);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextRt(i));

	//peak three
	i = make_pair(0,2);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,0)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextRt(i));

	//peak four
	i = make_pair(0,3);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextRt(i));

	//peak five
	i = make_pair(0,4);
	t.getNextRt(i);
	TEST_EQUAL(i.first,1)
	TEST_EQUAL(i.second,1)
	t.getNextRt(i);
	TEST_EQUAL(i.first,2)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getNextRt(i));

#ifdef OPENMS_DEBUG
	//test for corrupt index
	i = make_pair(5,0);
	TEST_EXCEPTION(Exception::Precondition, t.getNextRt(i));
	i = make_pair(1,5);
	TEST_EXCEPTION(Exception::Precondition, t.getNextRt(i));
#endif
RESULT

CHECK(void getPrevRt(IDX& index) throw (NoSuccessor))
	FeaFiTraits t;
	MSExperiment<> exp2 = exp;
	exp2[1].resize(4);
	exp2[1][0].setPos(599.0);
	exp2[1][1].setPos(799.0);
	exp2[1][2].setPos(801.0);
	exp2[1][3].setPos(1000.0);
	t.setData(exp2.begin(),exp2.end(),2);
	FeaFiTraits::IDX i;
	
	//peak one
	i = make_pair(1,0);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,0)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevRt(i));

	//peak two
	i = make_pair(1,1);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,1)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevRt(i));

	//peak three
	i = make_pair(1,2);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevRt(i));

	//peak four
	i = make_pair(1,3);
	t.getPrevRt(i);
	TEST_EQUAL(i.first,0)
	TEST_EQUAL(i.second,2)
	TEST_EXCEPTION(FeaFiTraits::NoSuccessor, t.getPrevRt(i));

#ifdef OPENMS_DEBUG
	//test for corrupt index
	i = make_pair(5,0);
	TEST_EXCEPTION(Exception::Precondition, t.getPrevRt(i));
	i = make_pair(1,5);
	TEST_EXCEPTION(Exception::Precondition, t.getPrevRt(i));
#endif
RESULT

CHECK(void addConvexHull(const IndexSet& set, DFeature<2>& f) const)
	FeaFiTraits t;
	DPeak<2> p;
	DPeakArray<2> peak_array;
	p.getPosition()[0] = 1240.54;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1241.81;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 693.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 695.2;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.07;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1247.09;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1249.76;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 687.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 689.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.6;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 694.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1247.75;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1249.12;   p.getPosition()[1] = 688;       peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 689.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 691;       peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;   p.getPosition()[1] = 692.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;   p.getPosition()[1] = 693.8;     peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;   p.getPosition()[1] = 695.2;     peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1243.78;   p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1239.9;    p.getPosition()[1] = 695.4;     peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;   p.getPosition()[1] = 692;       peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;   p.getPosition()[1] = 691;       peak_array.push_back(p);
	p.getPosition()[0] = 1237.93;   p.getPosition()[1] = 688.4;     peak_array.push_back(p);
	
	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	t.setData(exp.begin(), exp.end(),100);
	
	FeaFiTraits::IndexSet set;
	for (UnsignedInt i=0; i<exp.size(); ++i) 
	{
		for (UnsignedInt j=0; j<exp[i].size(); ++j) 
		{
			set.insert(std::make_pair(i,j));
		}
	}
	
	DFeature<2> f;
	t.addConvexHull(set,f);
	DConvexHull<2>& hull = f.getConvexHulls()[0];
	TEST_REAL_EQUAL(hull.getPoints().size(), 9);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1237.27, 691)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1237.93, 688.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1240.54, 687.6)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1252.39, 687.6)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1252.39, 693.8)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1251.73, 695.2)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1251.07, 695.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1239.9, 695.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), DConvexHull<2>::PointType(1237.27, 692)) != hull.getPoints().end(), true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
