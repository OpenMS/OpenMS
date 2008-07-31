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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/TimeStamp.h>
#include <fstream>

///////////////////////////

START_TEST(PreciseTime, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// tests for class PreciseTime::

#define BUSY_WAIT \
	double x = 0.0;  for (int i = 0; i < 20000000; i++, x += rand());

PreciseTime* t_ptr;
CHECK(PreciseTime::PreciseTime())
	t_ptr = new PreciseTime;
	TEST_NOT_EQUAL(t_ptr, 0)
RESULT

CHECK(PreciseTime::~PreciseTime())
	delete t_ptr;
RESULT

CHECK(PreciseTime::getSeconds() const )
	PreciseTime t;
	TEST_EQUAL(t.getSeconds(), 0)
RESULT


CHECK(PreciseTime::getMicroSeconds() const )
	PreciseTime t;
	TEST_EQUAL(t.getMicroSeconds(), 0)
RESULT


CHECK(PreciseTime::set(long secs, long usecs) )
	PreciseTime t;
	TEST_EQUAL(t.getSeconds(), 0)
	TEST_EQUAL(t.getMicroSeconds(), 0)
	t.set(1,1);
	TEST_EQUAL(t.getSeconds(), 1)
	TEST_EQUAL(t.getMicroSeconds(), 1)
	t.set(9999,12345);
	TEST_EQUAL(t.getSeconds(), 9999)
	TEST_EQUAL(t.getMicroSeconds(), 12345)
RESULT


CHECK(PreciseTime::PreciseTime(const PreciseTime& time))
	PreciseTime t1;
	t1.set(12345678, 456789);
	PreciseTime t2(t1);
	TEST_EQUAL(t2, t1)
	TEST_EQUAL(t2.getSeconds(), 12345678)
	TEST_EQUAL(t2.getMicroSeconds(), 456789)
RESULT

CHECK(PreciseTime::set(const PreciseTime& time) )	
	PreciseTime t1, t2;
	t1.set(12345678, 456789);
	t2.set(t1);
	TEST_EQUAL(t2, t1)
	TEST_EQUAL(t2.getSeconds(), 12345678)
	TEST_EQUAL(t2.getMicroSeconds(), 456789)
RESULT


CHECK(PreciseTime::PreciseTime& operator = (const PreciseTime& time) )
	PreciseTime t1, t2;
	t1.set(12345678, 456789);
	t2 = t1;
	TEST_EQUAL(t2, t1)
	TEST_EQUAL(t2.getSeconds(), 12345678)
	TEST_EQUAL(t2.getMicroSeconds(), 456789)
RESULT

CHECK(void PreciseTime::clear() )
	PreciseTime t1;
	PreciseTime t2;
	TEST_EQUAL(t1, t2)
	TEST_EQUAL(t1.getSeconds(), 0)
	TEST_EQUAL(t1.getMicroSeconds(), 0)
	t1.set(12345, 23456);
	TEST_EQUAL(t1.getSeconds(), 12345)
	TEST_EQUAL(t1.getMicroSeconds(), 23456)
	t1.clear();
	TEST_EQUAL(t1, t2)
RESULT


CHECK(PreciseTime::bool operator < (const PreciseTime& time) const  )
	PreciseTime t1, t2;
	t1.set(12345678, 456789);
	t2.set(12345679, 456789);
	TEST_EQUAL((t2 < t1), false)
	TEST_EQUAL((t1 < t2), true)
	t2.set(12345678, 456789);
	TEST_EQUAL((t2 < t1), false)
	TEST_EQUAL((t1 < t2), false)
	t2.set(12345678, 2345);
	TEST_EQUAL((t2 < t1), true)
	TEST_EQUAL((t1 < t2), false)
RESULT


CHECK(PreciseTime::bool operator > (const PreciseTime& time) const  )
	PreciseTime t1, t2;
	t1.set(12345678, 456789);
	t2.set(12345679, 456789);
	TEST_EQUAL((t2 > t1), true)
	TEST_EQUAL((t1 > t2), false)
	t2.set(12345678, 456789);
	TEST_EQUAL((t2 > t1), false)
	TEST_EQUAL((t1 > t2), false)
	t2.set(12345678, 2345);
	TEST_EQUAL((t2 > t1), false)
	TEST_EQUAL((t1 > t2), true)
RESULT


CHECK(PreciseTime::bool operator == (const PreciseTime& time) const  )
	PreciseTime t1, t2;
	t1.set(12345678, 456789);
	t2.set(12345679, 456789);
	TEST_EQUAL((t2 == t1), false)
	TEST_EQUAL((t1 == t2), false)
	t2.set(12345678, 456789);
	TEST_EQUAL((t2 == t1), true)
	TEST_EQUAL((t1 == t2), true)
	t2.set(12345678, 2345);
	TEST_EQUAL((t2 == t1), false)
	TEST_EQUAL((t1 == t2), false)
RESULT


CHECK(PreciseTime::now())
	PreciseTime t1(PreciseTime::now());
	TEST_NOT_EQUAL(t1.getSeconds(), 0)
	TEST_NOT_EQUAL(t1.getMicroSeconds(), 0)
	BUSY_WAIT
	PreciseTime t2(PreciseTime::now());
	TEST_NOT_EQUAL(t2.getSeconds(), 0)
	TEST_NOT_EQUAL(t2.getMicroSeconds(), 0)
	STATUS(t1.getSeconds() << "/" << t1.getMicroSeconds())
	STATUS(t2.getSeconds() << "/" << t2.getMicroSeconds())
	TEST_EQUAL((t1 < t2), true)
	TEST_EQUAL((t1 == t2), false)
RESULT

CHECK(ostream& operator << (ostream& os, const PreciseTime& time))
	PreciseTime t(12345678, 456789);
	string filename;
	NEW_TMP_FILE(filename);
	ofstream of(filename.c_str(), std::ios::out);
	of << t << std::endl;
	of.close();
	// ???? This still has to be ported from BALL
	// TEST_FILE_REGEXP(filename.c_str(), "data/PreciseTime_test.txt")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
