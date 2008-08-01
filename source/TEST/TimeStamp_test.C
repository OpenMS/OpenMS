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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CONCEPT/TimeStamp.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <fstream>

///////////////////////////

START_TEST(TimeStamp, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#define BUSY_WAIT { double x = 0.0; for (int i = 0; i < 2000000; i++, x += rand()); } 

using namespace OpenMS;
using namespace std;

TimeStamp* ts = 0;
CHECK(TimeStamp())
	ts = new TimeStamp;
	TEST_NOT_EQUAL(ts, 0)
RESULT


CHECK(~TimeStamp())
	delete ts;
RESULT


CHECK(getTime() const  )
  TimeStamp* t1 = new TimeStamp;	
	t1->stamp();
	STATUS(*t1)
	BUSY_WAIT
  TimeStamp* t2 = new TimeStamp;
	t2->stamp();
	STATUS(*t2)
	TEST_NOT_EQUAL(t1->getTime(), t2->getTime())
	TEST_EQUAL((t1->getTime() < t2->getTime()), true)
	TEST_EQUAL((t1->getTime() > t2->getTime()), false)
	delete t1;
	delete t2;
RESULT

CHECK(isNewerThan(const Time& time) const  )
	TimeStamp* ts1 = new TimeStamp;
	ts1->stamp();
	STATUS(*ts1)
	BUSY_WAIT
	TimeStamp* ts2 = new TimeStamp;
	ts2->stamp();
	STATUS(*ts2)
	TEST_EQUAL(ts1->isNewerThan(ts1->getTime()), false)
	TEST_EQUAL(ts1->isNewerThan(ts2->getTime()), false)
	TEST_EQUAL(ts2->isNewerThan(ts1->getTime()), true)
	TEST_EQUAL(ts2->isNewerThan(ts2->getTime()), false)
	delete ts1;
	delete ts2;
RESULT


CHECK(isOlderThan(const Time& time) const  )
	TimeStamp* ts1 = new TimeStamp;
	ts1->stamp();
	STATUS(*ts1)
	BUSY_WAIT
	TimeStamp* ts2 = new TimeStamp;
	ts2->stamp();
	STATUS(*ts2)
	TEST_EQUAL(ts1->isOlderThan(ts1->getTime()), false)
	TEST_EQUAL(ts1->isOlderThan(ts2->getTime()), true)
	TEST_EQUAL(ts2->isOlderThan(ts1->getTime()), false)
	TEST_EQUAL(ts2->isOlderThan(ts2->getTime()), false)
	delete ts1;
	delete ts2;
RESULT


CHECK(isNewerThan(const TimeStamp& stamp) const  )
	TimeStamp* ts1 = new TimeStamp;
	ts1->stamp();
	STATUS(*ts1)
	BUSY_WAIT
	TimeStamp* ts2 = new TimeStamp;
	ts2->stamp();
	STATUS(*ts2)
	TEST_EQUAL(ts1->isNewerThan(*ts1), false)
	TEST_EQUAL(ts1->isNewerThan(*ts2), false)
	TEST_EQUAL(ts2->isNewerThan(*ts1), true)
	TEST_EQUAL(ts2->isNewerThan(*ts2), false)
	delete ts1;
	delete ts2;
RESULT


CHECK(isOlderThan(const TimeStamp& stamp) const  )
	TimeStamp* ts1 = new TimeStamp;
	ts1->stamp();
	STATUS(*ts1)
	BUSY_WAIT
	TimeStamp* ts2 = new TimeStamp;
	ts2->stamp();
	STATUS(*ts2)
	TEST_EQUAL(ts1->isOlderThan(*ts1), false)
	TEST_EQUAL(ts1->isOlderThan(*ts2), true)
	TEST_EQUAL(ts2->isOlderThan(*ts1), false)
	TEST_EQUAL(ts2->isOlderThan(*ts2), false)
	delete ts1;
	delete ts2;
RESULT


CHECK(stamp(const Time& time = ZERO) )
  TimeStamp* ts1 = new TimeStamp;
	ts1->stamp();
	STATUS(*ts1)
	BUSY_WAIT
  TimeStamp* ts2 = new TimeStamp;
	ts2->stamp();
	STATUS(*ts2)
	TEST_EQUAL(ts1->isNewerThan(*ts1), false)
	TEST_EQUAL(ts1->isNewerThan(*ts2), false)
	TEST_EQUAL(ts2->isNewerThan(*ts1), true)
	TEST_EQUAL(ts2->isNewerThan(*ts2), false)
	BUSY_WAIT
	ts1->stamp();
	TEST_EQUAL(ts1->isNewerThan(*ts1), false)
	TEST_EQUAL(ts1->isNewerThan(*ts2), true)
	TEST_EQUAL(ts2->isNewerThan(*ts1), false)
	TEST_EQUAL(ts2->isNewerThan(*ts2), false)
	delete ts1;
	delete ts2;
RESULT

CHECK(operator << (std::ostream& os, const TimeStamp& ts))
	TimeStamp t;
	// a very nasty way to break the encapsulation, but simplifies
	// things a great deal....!
	PreciseTime& t_ref = const_cast<PreciseTime&>(t.getTime());
	t_ref.set(12345678, 456789);
	stringstream str;
	str << t;
	String time_str = str.str();
	TEST_EQUAL(time_str.hasPrefix("1970052"), true)
	TEST_EQUAL(time_str.hasSuffix("456789"), true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
