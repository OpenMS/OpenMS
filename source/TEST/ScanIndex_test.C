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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/ScanIndex.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/DPeakArray.h>

///////////////////////////

START_TEST(ScanIndex, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPeakArray<2> PeakVector;

/////////////////////////////////////////////////////////////

ScanIndex<PeakVector>* index_ptr = 0;

CHECK(ScanIndex())
	index_ptr = new ScanIndex<PeakVector>;
	TEST_NOT_EQUAL(index_ptr, 0)
	TEST_EQUAL(index_ptr->size(), 0)
RESULT

CHECK(~ScanIndex())
	delete index_ptr;
RESULT

PeakVector dpa;

DPeak<2> peak1;
DPosition<2> pos1;
pos1[0] = 1.0;
pos1[1] = 1.0;
peak1.getIntensity() = 1.0;
peak1.getPosition() = pos1;
dpa.push_back(peak1);

DPeak<2> peak2;
DPosition<2> pos2;
pos2[0] = 2.0;
pos2[1] = 2.0;
peak2.getIntensity() = 1.0;
peak2.getPosition() = pos2;
dpa.push_back(peak2);

DPeak<2> peak3;
DPosition<2> pos3;
pos3[0] = 3.0;
pos3[1] = 3.0;
peak3.getIntensity() = 2.0;
peak3.getPosition() = pos3;
dpa.push_back(peak3);

DPeak<2> peak4;
DPosition<2> pos4;
pos4[0] = 4.0;
pos4[1] = 4.0;
peak4.getIntensity() = 3.0;
peak4.getPosition() = pos4;
dpa.push_back(peak4);

ScanIndex<PeakVector> scani;
	
CHECK(init())
	scani.init(dpa.begin(),dpa.end());
	TEST_EQUAL(scani.size(),5);	
RESULT

CHECK(getRank())
	TEST_EQUAL(scani.getRank(pos1[0]),0);	
	TEST_EQUAL(scani.getRank(pos2[0]),1);	
	TEST_EQUAL(scani.getRank(pos3[0]),2);
	TEST_EQUAL(scani.getRank(pos4[0]),3);
RESULT

/////////////////////////////////////////////////////////////
END_TEST


