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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(BaseFeaFiTraits, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;


// default ctor
FeaFiTraits* ptr = 0;
CHECK(FeaFiTraits())
	ptr = new FeaFiTraits();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~FeaFiTraits())
	delete ptr;
RESULT

CHECK(const ConvexHullType& FeaFiTraits::calculateConvexHull(const IndexSet& set))
	FeaFiTraits t;
	DPeak<2> p;
	DPeakArray<2> peak_array;

	p.getPosition()[0] = 1240.54; p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1241.81; p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 692.8;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 693.8;	peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;	p.getPosition()[1] = 695.2;	peak_array.push_back(p);
	p.getPosition()[0] = 1251.07;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1247.09;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;	p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1249.76;	p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;	p.getPosition()[1] = 687.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 689.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 692.6;	peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;	p.getPosition()[1] = 694.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1247.75;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1249.12;	p.getPosition()[1] = 688;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 689.8;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 691;	peak_array.push_back(p);
	p.getPosition()[0] = 1252.39;	p.getPosition()[1] = 692.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1251.73;	p.getPosition()[1] = 693.8;	peak_array.push_back(p);
	p.getPosition()[0] = 1250.41;	p.getPosition()[1] = 695.2;	peak_array.push_back(p);
	p.getPosition()[0] = 1248.41;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1243.78;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1239.9;	p.getPosition()[1] = 695.4;	peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;	p.getPosition()[1] = 692;	peak_array.push_back(p);
	p.getPosition()[0] = 1237.27;	p.getPosition()[1] = 691;	peak_array.push_back(p);
	p.getPosition()[0] = 1237.93;	p.getPosition()[1] = 688.4;	peak_array.push_back(p);

	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	t.setData(exp);
	
	IndexSet set;
	set.add(0,27);
	FeaFiTraits::ConvexHullType hull = t.calculateConvexHull(set);
	TEST_REAL_EQUAL(hull.getPoints().size(), 9);
	
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1237.27, 691)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1237.93, 688.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1240.54, 687.6)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1252.39, 687.6)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1252.39, 693.8)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1251.73, 695.2)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1251.07, 695.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1239.9, 695.4)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1237.27, 692)) != hull.getPoints().end(), true);
	
RESULT

CHECK(const ConvexHullType& FeaFiTraits::calculateConvexHull(const IndexSet& set))
	FeaFiTraits t;
	DPeak<2> p;
	DPeakArray<2> peak_array;
	p.getPosition()[0] = 61.14; p.getPosition()[1] = 429.242;	peak_array.push_back(p);
	p.getPosition()[0] = 61.14; p.getPosition()[1] = 429.266;	peak_array.push_back(p);
	p.getPosition()[0] = 61.14;	p.getPosition()[1] = 429.291;	peak_array.push_back(p);
	p.getPosition()[0] = 61.14;	p.getPosition()[1] = 429.315;	peak_array.push_back(p);
	p.getPosition()[0] = 64.36;	p.getPosition()[1] = 429.242;	peak_array.push_back(p);
	p.getPosition()[0] = 64.36;	p.getPosition()[1] = 429.266;	peak_array.push_back(p);
	p.getPosition()[0] = 64.36;	p.getPosition()[1] = 429.315;	peak_array.push_back(p);
	p.getPosition()[0] = 64.36;	p.getPosition()[1] = 429.389;	peak_array.push_back(p);
	p.getPosition()[0] = 64.36;	p.getPosition()[1] = 429.437;	peak_array.push_back(p);
	
	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	t.setData(exp);

	IndexSet set;
	set.add(0,8);
	FeaFiTraits::ConvexHullType hull = t.calculateConvexHull(set);
	TEST_REAL_EQUAL(hull.getPoints().size(), 4);

	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(61.14, 429.242)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(64.36, 429.242)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(64.36, 429.437)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(61.14, 429.315)) != hull.getPoints().end(), true);

RESULT

CHECK(const ConvexHullType& FeaFiTraits::calculateConvexHull(const IndexSet& set))
	FeaFiTraits t;
	DPeak<2> p;
	DPeakArray<2> peak_array;
	p.getPosition()[0] = 51.51; p.getPosition()[1] = 428.778;	peak_array.push_back(p);
	p.getPosition()[0] = 51.51; p.getPosition()[1] = 428.802;	peak_array.push_back(p);
	p.getPosition()[0] = 51.51; p.getPosition()[1] = 428.851;	peak_array.push_back(p);
	p.getPosition()[0] = 51.51; p.getPosition()[1] = 428.876;	peak_array.push_back(p);
	p.getPosition()[0] = 51.51; p.getPosition()[1] = 428.9;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.729;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.754;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.778;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.827;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.876;	peak_array.push_back(p);
	p.getPosition()[0] = 54.72; p.getPosition()[1] = 428.924;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.754;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.778;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.802;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.827;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.851;	peak_array.push_back(p);
	p.getPosition()[0] = 57.93; p.getPosition()[1] = 428.9;	peak_array.push_back(p);

	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	t.setData(exp);
	
	IndexSet set;
	set.add(0,16);
	FeaFiTraits::ConvexHullType hull = t.calculateConvexHull(set);
	TEST_REAL_EQUAL(hull.getPoints().size(), 6);

	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(51.51, 428.778)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(54.72, 428.729)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(57.93, 428.754)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(57.93, 428.9)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(54.72, 428.924)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(51.51, 428.9)) != hull.getPoints().end(), true);	
RESULT

CHECK(const ConvexHullType& BaseFeaFiTraits::calculateConvexHull(const IndexSet& set))
	FeaFiTraits t;
	DPeak<2> p;
	DPeakArray<2> peak_array;
	p.getPosition()[0] = 1.0; p.getPosition()[1] = 3.0;	peak_array.push_back(p);
	p.getPosition()[0] = 1.0; p.getPosition()[1] = 0.0;	peak_array.push_back(p);
	p.getPosition()[0] = 0.0;	p.getPosition()[1] = 1.0;	peak_array.push_back(p);
	p.getPosition()[0] = 2.0;	p.getPosition()[1] = 0.0;	peak_array.push_back(p);
	p.getPosition()[0] = 2.0;	p.getPosition()[1] = 2.0;	peak_array.push_back(p);
	
	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	t.setData(exp);

	IndexSet set;
	set.add(0,4);
	FeaFiTraits::ConvexHullType hull = t.calculateConvexHull(set);
	TEST_REAL_EQUAL(hull.getPoints().size(), 5);

	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(0.0, 1.0)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1.0, 0.0)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(2.0, 0.0)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(2.0, 2.0)) != hull.getPoints().end(), true);
	TEST_EQUAL(find(hull.getPoints().begin(), hull.getPoints().end(), FeaFiTraits::ConvexHullType::PointType(1.0, 3.0)) != hull.getPoints().end(), true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
