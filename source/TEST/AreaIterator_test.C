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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/AreaIterator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AreaIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef MSExperiment<> Map;
typedef Internal::AreaIterator<Map::PeakType, Map::PeakType&, Map::PeakType*, Map::Iterator, Map::SpectrumType::Iterator> AI;
typedef Internal::AreaIterator<Map::PeakType, const Map::PeakType&, const Map::PeakType*, Map::ConstIterator, Map::SpectrumType::ConstIterator> CAI;

AI* ptr1 = 0, *ptr2 = 0;

Map exp;
exp.resize(5);
exp[0].resize(2);
exp[0].setRetentionTime(2.0);
exp[0].setMSLevel(1);
exp[0][0].setPos(502);
exp[0][1].setPos(510);

exp[1].resize(2);
exp[1].setRetentionTime(4.0);
exp[1].setMSLevel(1);
exp[1][0].setPos(504);
exp[1][1].setPos(506);

exp[2].setRetentionTime(6.0);
exp[2].setMSLevel(1);

exp[3].resize(2);
exp[3].setRetentionTime(8.0);
exp[3].setMSLevel(1);
exp[3][0].setPos(504.1);
exp[3][1].setPos(506.1);

exp[4].resize(2);
exp[4].setRetentionTime(10.0);
exp[4].setMSLevel(1);
exp[4][0].setPos(502.1);
exp[4][1].setPos(510.1);

CHECK(AreaIterator(const SpectrumIteratorType& spectrum_end, const PeakIteratorType& peak_end ))
	ptr1 = new AI(exp.end(),exp.back().end());
RESULT

CHECK(AreaIterator(const SpectrumIteratorType& begin, const SpectrumIteratorType& end, CoordinateType low_mz, CoordinateType high_mz))
	ptr2 = new AI(exp.RTBegin(0), exp.RTEnd(0), 0, 0);
RESULT

CHECK(~AreaIterator())
	delete ptr1;
	delete ptr2;
RESULT

CHECK(Tut ConstAreaIterator was er soll )
	//whole area
	AI it = AI(exp.RTBegin(0), exp.RTEnd(15), 500, 520);
	TEST_REAL_EQUAL(it->getPos(),502.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),510.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),504.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),506.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),504.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),506.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),502.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),510.1);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);
	
	//center peaks
	it = AI(exp.RTBegin(3), exp.RTEnd(9), 503, 509);
	TEST_REAL_EQUAL(it->getPos(),504.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),506.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),504.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),506.1);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);
	
	//upper left area
	it = AI(exp.RTBegin(0), exp.RTEnd(7), 505, 520);
	TEST_REAL_EQUAL(it->getPos(),510.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),506.0);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);
	
	//upper right area
	it = AI(exp.RTBegin(5), exp.RTEnd(11), 505, 520);
	TEST_REAL_EQUAL(it->getPos(),506.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),510.1);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);
	
	//lower right
	it = AI(exp.RTBegin(5), exp.RTEnd(11), 500, 505);
	TEST_REAL_EQUAL(it->getPos(),504.1);
	++it;
	TEST_REAL_EQUAL(it->getPos(),502.1);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);

	//lower left
	it = AI(exp.RTBegin(0), exp.RTEnd(7), 500, 505);
	TEST_REAL_EQUAL(it->getPos(),502.0);
	++it;
	TEST_REAL_EQUAL(it->getPos(),504.0);
	++it;
	TEST_EQUAL(it==exp.areaEnd(),true);

	//Test with empty RT range
	it = AI(exp.RTBegin(5), exp.RTEnd(5.5), 500, 520);
	TEST_EQUAL(it==exp.areaEnd(),true);

	//Test with empty MZ range
	it = AI(exp.RTBegin(0), exp.RTEnd(15), 505, 505.5);
	TEST_EQUAL(it==exp.areaEnd(),true);

	//Test with empty RT + MZ range
	it = AI(exp.RTBegin(5), exp.RTEnd(5.5), 505, 505.5);
	TEST_EQUAL(it==exp.areaEnd(),true);

	//Test with empty (no MS level 1) experiment
	exp[0].setMSLevel(2);
	exp[1].setMSLevel(2);
	exp[2].setMSLevel(2);
	exp[3].setMSLevel(2);
	exp[4].setMSLevel(2);
	it = AI(exp.RTBegin(0), exp.RTEnd(15), 500, 520);
	TEST_EQUAL(it==exp.areaEnd(),true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
