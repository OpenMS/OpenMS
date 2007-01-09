// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(AreaIterator, "$Id: PeakFileOptions_test.C 1077 2006-12-15 13:39:54Z kadauke $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef DPeak<1> Peak;
typedef MSExperiment<Peak> Experiment;
typedef Experiment::AIterator AIterator;

AIterator* ptr1 = 0, *ptr2 = 0;
Experiment exp;

CHECK(AreaIterator())
	// initialize experiment
	Experiment::SpectrumType spec;
	spec.push_back(Peak());
	spec.setRetentionTime(0);
	exp.push_back(spec);

	ptr1 = new AIterator();
RESULT

CHECK(AreaIterator(SpectrumIteratorType begin, SpectrumIteratorType end, CoordinateType low_mz, CoordinateType high_mz))
	ptr2 = new AIterator(exp.RTBegin(0), exp.RTEnd(0), 0, 0);
RESULT

CHECK(~AreaIterator())
	delete ptr1;
	delete ptr2;
RESULT

CHECK(AreaIterator(const AreaIterator& rhs))
	AIterator iter = exp.areaEnd();
	AIterator iter2(iter);
	TEST_EQUAL(iter == iter2, true);
RESULT

CHECK(AreaIterator& operator=(const AreaIterator& rhs))
	AIterator iter = exp.areaEnd();
	AIterator iter2;
	iter2 = iter;
	TEST_EQUAL(iter == iter2, true);
RESULT

CHECK(bool operator==(const AreaIterator& rhs))
	AIterator iter, iter2;
	TEST_EQUAL(iter == iter2, true);
RESULT

CHECK(bool operator!=(const AreaIterator& rhs))
	AIterator iter, iter2;
	TEST_EQUAL(iter != iter2, false);
RESULT

CHECK(AreaIterator& operator++())
	DPosition<2> lower(0), upper(10);
	DRange<2> range(lower, upper);
	AIterator iter = exp.areaBegin(range);
	AIterator iter2(iter);
	++iter;
	TEST_EQUAL(iter == iter2, false);
	TEST_EQUAL(iter == exp.areaEnd(), true);
	TEST_EQUAL(iter2 == exp.areaEnd(), false);
RESULT

CHECK(AreaIterator operator++(int))
	DPosition<2> lower(0), upper(10);
	DRange<2> range(lower, upper);
	AIterator iter = exp.areaBegin(range);
	AIterator iter2(iter);
	iter++;
	TEST_EQUAL(iter == iter2, false);
	TEST_EQUAL(iter == exp.areaEnd(), true);
	TEST_EQUAL(iter2 == exp.areaEnd(), false);
RESULT

CHECK(reference operator*())
	DPosition<2> lower(0), upper(10);
	DRange<2> range(lower, upper);
	AIterator iter = exp.areaBegin(range);
	TEST_EQUAL(*iter == Peak(), true);
RESULT

CHECK(pointer operator->())
	DPosition<2> lower(0), upper(10);
	DRange<2> range(lower, upper);
	AIterator iter = exp.areaBegin(range);
	TEST_EQUAL(iter->getPos() == 0, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
