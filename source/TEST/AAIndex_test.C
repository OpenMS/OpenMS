// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
	
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>

using namespace OpenMS;
using namespace std;
	
///////////////////////////

AASequence seq1("ALEGDEK");
AASequence seq2("GTVVTGR");
AASequence seq3("EHVLLAR");


START_TEST(AASequenceIndeces, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//sequence spec_id protein_id mass GB500 arginin_count KHAG800101 VASM830103 NADH010106 NADH010107 WILM950102 ROBB760107 OOBM850104 FAUJ880111 FINA770101 ARGP820102 M F H Q Y target_log
//ALEGDEK 15 0587  761.368 1337.53 0  129.3 1.145   31  565  1.5200000 -6.60000e+00  -3.240000 1  7.18  5.23 0 0 0 0 0 2.08623342
//GTVVTGR 15 0587  689.394 1442.70 1  383.2 1.042  241  403  7.1800000 -3.00000e-01 -16.010000 1  5.55  5.02 0 0 0 0 0 1.35346120
//EHVLLAR 15 0587  837.494 1442.70 1  318.5 1.259  171  190 18.1300000  3.00000e-01  -9.970000 2  7.73  9.34 0 0 1 0 0 5.22075034

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(static DoubleReal calculateGB(const AASequence& seq, DoubleReal T=500.0) )
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq1), 1337.53)
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq2), 1442.70)
	TEST_REAL_SIMILAR(AAIndex::calculateGB(seq3), 1442.70)

	TEST_NOT_EQUAL(AAIndex::calculateGB(seq1,100.0), 1337.53)
	TEST_NOT_EQUAL(AAIndex::calculateGB(seq2,100.0), 1442.70)
	TEST_NOT_EQUAL(AAIndex::calculateGB(seq3,100.0), 1442.70)
END_SECTION

START_SECTION(static DoubleReal aliphatic(char aa))
	TEST_REAL_SIMILAR(AAIndex::aliphatic('A'),1.0)
	TEST_REAL_SIMILAR(AAIndex::aliphatic('B'),0.0)
END_SECTION

START_SECTION(static DoubleReal acidic(char aa))
	TEST_REAL_SIMILAR(AAIndex::acidic('D'),1.0)
	TEST_REAL_SIMILAR(AAIndex::acidic('A'),0.0)
END_SECTION

START_SECTION(static DoubleReal basic(char aa))
	TEST_REAL_SIMILAR(AAIndex::basic('K'),1.0)
	TEST_REAL_SIMILAR(AAIndex::basic('A'),0.0)
END_SECTION

START_SECTION(static DoubleReal polar(char aa))
	TEST_REAL_SIMILAR(AAIndex::polar('S'),1.0)
	TEST_REAL_SIMILAR(AAIndex::polar('A'),0.0)
END_SECTION

START_SECTION(static DoubleReal getKHAG800101(char aa))
 TEST_REAL_SIMILAR(AAIndex::getKHAG800101('A'),49.1)
END_SECTION

START_SECTION(static DoubleReal getVASM830103(char aa))
 TEST_REAL_SIMILAR(AAIndex::getVASM830103('A'),0.159)
END_SECTION

START_SECTION(static DoubleReal getNADH010106(char aa))
 TEST_REAL_SIMILAR(AAIndex::getNADH010106('A'),5.0)
END_SECTION

START_SECTION(static DoubleReal getNADH010107(char aa))
 TEST_REAL_SIMILAR(AAIndex::getNADH010107('A'),-2.0)
END_SECTION

START_SECTION(static DoubleReal getWILM950102(char aa))
 TEST_REAL_SIMILAR(AAIndex::getWILM950102('A'),2.62)
END_SECTION

START_SECTION(static DoubleReal getROBB760107(char aa))
 TEST_REAL_SIMILAR(AAIndex::getROBB760107('A'),0.0)
END_SECTION

START_SECTION(static DoubleReal getOOBM850104(char aa))
 TEST_REAL_SIMILAR(AAIndex::getOOBM850104('A'),-2.49)
END_SECTION

START_SECTION(static DoubleReal getFAUJ880111(char aa))
 TEST_REAL_SIMILAR(AAIndex::getFAUJ880111('A'),0.0)
END_SECTION

START_SECTION(static DoubleReal getFINA770101(char aa))
 TEST_REAL_SIMILAR(AAIndex::getFINA770101('A'),1.08)
END_SECTION

START_SECTION(static DoubleReal getARGP820102(char aa))
 TEST_REAL_SIMILAR(AAIndex::getARGP820102('A'),1.18)
END_SECTION


END_TEST
