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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/SYSTEM/StopWatch.h>

///////////////////////////

START_TEST(StopWatch, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

CHECK(StopWatch& operator = (const StopWatch& stop_watch))
  // ???
RESULT

CHECK(StopWatch())
  // ???
RESULT

CHECK(StopWatch(const StopWatch& stop_watch))
  // ???
RESULT

CHECK(bool isRunning() const throw())
  // ???
RESULT

CHECK(bool operator != (const StopWatch& stop_watch) const throw())
  // ???
RESULT

CHECK(bool operator < (const StopWatch& stop_watch) const throw())
  // ???
RESULT

CHECK(bool operator <= (const StopWatch& stop_watch) const throw())
  // ???
RESULT

CHECK(bool operator == (const StopWatch& stop_watch) const)
  // ???
RESULT

CHECK(bool operator > (const StopWatch& stop_watch) const throw())
  // ???
RESULT

CHECK(bool operator >= (const StopWatch& stop_watch) const throw())
  // ???
RESULT

CHECK(bool start())
  // ???
RESULT

CHECK(bool stop())
  // ???
RESULT

CHECK(double getCPUTime() const throw())
  // ???
RESULT

CHECK(double getClockTime() const)
  // ???
RESULT

CHECK(double getSystemTime() const)
  // ???
RESULT

CHECK(double getUserTime() const)
  // ???
RESULT

CHECK(void clear())
  // ???
RESULT

CHECK(void reset())
  // ???
RESULT

CHECK(~StopWatch())
  // ???
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
