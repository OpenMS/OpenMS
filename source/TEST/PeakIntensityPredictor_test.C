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
// $Maintainer: Alexandra Scherbart $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakIntensityPredictor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.001)

AASequence seq1("LTSEAR");
AASequence seq2("AEAQIR");
AASequence seq3("TLEDAR");

vector<AASequence> vec;
vec.push_back(seq1);
vec.push_back(seq2);
vec.push_back(seq3);

PeakIntensityPredictor* ptr;

START_SECTION(PeakIntensityPredictor())
	ptr = new PeakIntensityPredictor();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~PeakIntensityPredictor()))
	delete ptr;
END_SECTION

START_SECTION(DoubleReal predict(const AASequence& sequence))
	PeakIntensityPredictor pip;
	TEST_REAL_SIMILAR(pip.predict(seq1), -0.531675)
	TEST_REAL_SIMILAR(pip.predict(seq2), 0.0171194)
	TEST_REAL_SIMILAR(pip.predict(seq3), -0.595362)
END_SECTION


START_SECTION(DoubleReal predict(const AASequence& sequence, std::vector<DoubleReal>& add_info))
	PeakIntensityPredictor pip;
	std::vector<DoubleReal> add_info;
	pip.predict(seq1,add_info);
	TEST_EQUAL(add_info.size(),3)
	TEST_REAL_SIMILAR(add_info[0],0.0)
	TEST_REAL_SIMILAR(add_info[1],1.0)
	TEST_REAL_SIMILAR(add_info[2],2.04653)
END_SECTION

START_SECTION(std::vector<DoubleReal> predict(const std::vector<AASequence>& sequences))
	PeakIntensityPredictor pip;
	vector<DoubleReal> ref = pip.predict(vec);
	TEST_REAL_SIMILAR(ref[0], -0.531675)
	TEST_REAL_SIMILAR(ref[1], 0.0171194)
	TEST_REAL_SIMILAR(ref[2], -0.595362)
END_SECTION

START_SECTION(std::vector<DoubleReal> predict(const std::vector<AASequence>& sequences, std::vector<std::vector<DoubleReal> >& add_info))
	PeakIntensityPredictor pip;
	vector<vector<DoubleReal> > add_info;
	pip.predict(vec,add_info);
	TEST_EQUAL(add_info.size(),3)
	TEST_EQUAL(add_info[0].size(),3)
	TEST_EQUAL(add_info[1].size(),3)
	TEST_EQUAL(add_info[2].size(),3)
	TEST_REAL_SIMILAR(add_info[0][0],0.0)
	TEST_REAL_SIMILAR(add_info[0][1],1.0)
	TEST_REAL_SIMILAR(add_info[0][2],2.04653)
	TEST_REAL_SIMILAR(add_info[1][0],0.0)
	TEST_REAL_SIMILAR(add_info[1][1],1.0)
	TEST_REAL_SIMILAR(add_info[1][2],2.30648)
	TEST_REAL_SIMILAR(add_info[2][0],0.0)
	TEST_REAL_SIMILAR(add_info[2][1],1.0)
	TEST_REAL_SIMILAR(add_info[2][2],2.24984)
END_SECTION


END_TEST


