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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ConsensusIsotopePattern.h>

///////////////////////////

START_TEST(ConsensusIsotopePattern, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ConsensusIsotopePattern* ptr;
START_SECTION((ConsensusIsotopePattern()))
	ptr = new ConsensusIsotopePattern();
	TEST_NOT_EQUAL(ptr,0)
END_SECTION

START_SECTION((~ConsensusIsotopePattern()))
	delete ptr;
END_SECTION

ptr = new ConsensusIsotopePattern();

START_SECTION(TODO)
// 			// constructs the consensus pattern:
// 			void constructConsusPattern();
// 			// order an isotope trace in the correct cluster:
// 			void addIsotopeTrace(double, double);
// 			// condenses the pattern, make average peaks from the traces:
// 			void condensIsotopePattern(std::pair<std::vector<double>, std::vector<double> >*);
// 
// 			///////////////////////////////
// 			// start here all the get / set
// 			// function to access the
// 			// variables of the class
// 
// 			std::map<double, double>::iterator getConsensIsotopeIteratorStart();
// 			std::map<double, double>::iterator getConsensIsotopeIteratorEnd();
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
