// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//
	
	
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>

using namespace OpenMS;
using namespace std;
	
///////////////////////////

AASequence seq1("ALEGDEK");
AASequence seq2("GTVVTGR");
AASequence seq3("EHVLLAR");


START_TEST(AASequenceIndeces, "$Id: AASequenceIndeces_test.C 0001 2008-08-14 11:44:24 ascherbart $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//sequence spec_id protein_id mass GB500 arginin_count KHAG800101 VASM830103 NADH010106 NADH010107 WILM950102 ROBB760107 OOBM850104 FAUJ880111 FINA770101 ARGP820102 M F H Q Y target_log
//ALEGDEK 15 0587  761.368 1337.53 0  129.3 1.145   31  565  1.5200000 -6.60000e+00  -3.240000 1  7.18  5.23 0 0 0 0 0 2.08623342
//GTVVTGR 15 0587  689.394 1442.70 1  383.2 1.042  241  403  7.1800000 -3.00000e-01 -16.010000 1  5.55  5.02 0 0 0 0 0 1.35346120
//EHVLLAR 15 0587  837.494 1442.70 1  318.5 1.259  171  190 18.1300000  3.00000e-01  -9.970000 2  7.73  9.34 0 0 1 0 0 5.22075034

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(static vector<DoubleReal> getPropertyVector(AASequence& sequence))
	vector<DoubleReal> calculated = AAIndex::getPropertyVector(seq3);
	TEST_REAL_SIMILAR(calculated[0], 1.0);
	TEST_REAL_SIMILAR(calculated[1], 9.34);
	TEST_REAL_SIMILAR(calculated[2], 0.0);
	TEST_REAL_SIMILAR(calculated[3], 2.0);
	TEST_REAL_SIMILAR(calculated[4], 7.73);
	TEST_REAL_SIMILAR(calculated[5], 1442.7);
	TEST_REAL_SIMILAR(calculated[6], 1.0);
	TEST_REAL_SIMILAR(calculated[7], 318.5);
	TEST_REAL_SIMILAR(calculated[8], 0.0);
	TEST_REAL_SIMILAR(calculated[9], 836.978);		
	TEST_REAL_SIMILAR(calculated[10], 171.0);
	TEST_REAL_SIMILAR(calculated[11], 190.0);
	TEST_REAL_SIMILAR(calculated[12], -9.97);
	TEST_REAL_SIMILAR(calculated[13], 0.0);
	TEST_REAL_SIMILAR(calculated[14], 0.3);
	TEST_REAL_SIMILAR(calculated[15], 1.259);
	TEST_REAL_SIMILAR(calculated[16], 18.13);
	TEST_REAL_SIMILAR(calculated[17], 0.0);
		
END_SECTION

START_SECTION(static DoubleReal calculateGB(AASequence& seq, DoubleReal T=500.0) )
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq1), 1337.53)
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq2), 1442.70)
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq3), 1442.70)

	TEST_NOT_EQUAL(AAIndex::calculateGB(seq1,100.0), 1337.53)
	TEST_NOT_EQUAL(AAIndex::calculateGB(seq2,100.0), 1442.70)
	TEST_NOT_EQUAL(AAIndex::calculateGB(seq3,100.0), 1442.70)
END_SECTION


END_TEST
	
