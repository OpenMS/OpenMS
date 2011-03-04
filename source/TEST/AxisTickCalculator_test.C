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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;

START_TEST(AxisTickCalculator, "$Id$")

/////////////////////////////////////////////////////////////


START_SECTION((static void calcGridLines(DoubleReal x1, DoubleReal x2, GridVector& grid)))
	std::vector<std::vector<DoubleReal> > vector1;
  AxisTickCalculator::calcGridLines(1.0,4.0,vector1);
	
  TEST_EQUAL(2,vector1.size());
	TEST_EQUAL(4,vector1[0].size());
	
	TEST_REAL_SIMILAR(1.0,vector1[0][0]);
	TEST_REAL_SIMILAR(2.0,vector1[0][1]);
	TEST_REAL_SIMILAR(3.0,vector1[0][2]);
	TEST_REAL_SIMILAR(4.0,vector1[0][3]);
		
END_SECTION        

START_SECTION((static void calcLogGridLines(DoubleReal x1, DoubleReal x2, GridVector& grid)))
	std::vector<std::vector<DoubleReal> > vector1;
	AxisTickCalculator::calcLogGridLines(log10(10.0),log10(100.0),vector1);

	TEST_EQUAL(1,vector1[0].size());
	TEST_EQUAL(8,vector1[1].size());
	
	TEST_REAL_SIMILAR(1.0,vector1[0][0]);
	
	TEST_REAL_SIMILAR( 1.30103 ,vector1[1][0]);
	TEST_REAL_SIMILAR( 1.47712 ,vector1[1][1]);
	TEST_REAL_SIMILAR( 1.60206 ,vector1[1][2]);
	TEST_REAL_SIMILAR( 1.69897 ,vector1[1][3]);
	TEST_REAL_SIMILAR( 1.77815 ,vector1[1][4]);
	TEST_REAL_SIMILAR( 1.8451 ,vector1[1][5]);
	TEST_REAL_SIMILAR( 1.90309 ,vector1[1][6]);
	TEST_REAL_SIMILAR( 1.95425 ,vector1[1][7]);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



