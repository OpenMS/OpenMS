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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------



#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/RunningAveragePosition.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(RunningAverage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


RunningAveragePosition< DPosition<2> >* ptr = 0;
CHECK(RunningAverage())
	ptr = new RunningAveragePosition< DPosition<2> >();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Class())
	delete ptr;
RESULT

RunningAveragePosition< DPosition<4> > run;
	
DPosition<4> pos1;
pos1[0] = 1.0;
pos1[1] = 2.0;
pos1[2] = 3.0;
pos1[3] = 4.0;
	
DPosition<4> pos2;
pos2[0] = 1.0;
pos2[1] = 2.0;
pos2[2] = 3.0;
pos2[3] = 4.0;
	
CHECK(add())
		
	run.add(pos1,2);
	run.add(pos2,2);
	
	TEST_EQUAL(run.getPosition()[0],1)
	TEST_EQUAL(run.getPosition()[1],2)	
	TEST_EQUAL(run.getPosition()[2],3)
	TEST_EQUAL(run.getPosition()[3],4)
	
RESULT

CHECK(substract())
		
	run.subtract(pos2,2);
		
	TEST_EQUAL(run.getPosition()[0],1)
	TEST_EQUAL(run.getPosition()[1],2)	
	TEST_EQUAL(run.getPosition()[2],3)
	TEST_EQUAL(run.getPosition()[3],4)
	
	run.subtract(pos1,2);
		
	TEST_EQUAL(run.getPosition()[0],0)
	TEST_EQUAL(run.getPosition()[1],0)	
	TEST_EQUAL(run.getPosition()[2],0)
	TEST_EQUAL(run.getPosition()[3],0)
	
	TEST_EQUAL(run.getWeight(),0)
	
RESULT

CHECK(clear())
		
	run.add(pos2,2);
	run.add(pos1,2);
	run.clear();
	
	TEST_EQUAL(run.getWeight(),0)
	
	TEST_EQUAL(run.getPosition()[0],0)
	TEST_EQUAL(run.getPosition()[1],0)	
	TEST_EQUAL(run.getPosition()[2],0)
	TEST_EQUAL(run.getPosition()[3],0)
	
RESULT

END_TEST


