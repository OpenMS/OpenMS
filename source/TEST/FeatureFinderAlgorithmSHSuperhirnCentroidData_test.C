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
// $Maintainer: Florian Zeller$
// $Authors: Florian Zeller$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>

///////////////////////////

START_TEST(CentroidData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CentroidData* ptr;

vector<double>* centroidMasses = new vector<double>(); // Centroided masses
vector<double>* centroidIntens = new vector<double>(); // Centroided intensities
RawData* raw = new RawData(*centroidMasses, *centroidIntens);

START_SECTION((CentroidData()))
	ptr = new CentroidData(1, *raw, true);
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((~CentroidData()))
	delete ptr;
END_SECTION

ptr = new CentroidData(1, *raw, true);


START_SECTION((setAndGet()))
ptr = new CentroidData(1, *raw, true);

centroidMasses->push_back(1.0);
centroidIntens->push_back(2.0);

ptr->set(*centroidMasses, *centroidIntens);

list<CentroidPeak>* centroidPeaks = new list<CentroidPeak>();
ptr->get(*centroidPeaks);
TEST_EQUAL(centroidPeaks->size(), 1);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
