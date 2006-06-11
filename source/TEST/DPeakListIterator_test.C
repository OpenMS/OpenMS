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
// $Id: DPeakListIterator_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeakList.h>

///////////////////////////

START_TEST(DPeakList<D>, "$Id: DPeakListIterator_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

//construct a peak list to test on
DPeakList<1> dpl;
for(unsigned int i=1;i<11;i++)
{
	DPeak<1> peak;
	peak.getPosition()[0]=i;
	peak.getIntensity()=pow((float)-2,(int)i);
	dpl.push_back(peak);
}

CHECK( operator * () / operator ++ ())
DPeakList<1>::Iterator it;
it= dpl.begin();
TEST_REAL_EQUAL((*it).getPosition()[0], 1)
TEST_REAL_EQUAL((*it).getIntensity(), -2.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 2)
TEST_REAL_EQUAL((*it).getIntensity(), 4.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 3)
TEST_REAL_EQUAL((*it).getIntensity(), -8.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 4)
TEST_REAL_EQUAL((*it).getIntensity(), 16.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 5)
TEST_REAL_EQUAL((*it).getIntensity(), -32.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 6)
TEST_REAL_EQUAL((*it).getIntensity(), 64.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 7)
TEST_REAL_EQUAL((*it).getIntensity(), -128.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 8)
TEST_REAL_EQUAL((*it).getIntensity(), 256.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 9)
TEST_REAL_EQUAL((*it).getIntensity(), -512.0)
++it;
TEST_REAL_EQUAL((*it).getPosition()[0], 10)
TEST_REAL_EQUAL((*it).getIntensity(), 1024.0)
RESULT

CHECK( operator -> () / operator -- ())
DPeakList<1>::Iterator it;
it= dpl.end();
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 10)
TEST_REAL_EQUAL(it->getIntensity(), 1024.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 9)
TEST_REAL_EQUAL(it->getIntensity(), -512.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 8)
TEST_REAL_EQUAL(it->getIntensity(), 256.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 7)
TEST_REAL_EQUAL(it->getIntensity(), -128.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 6)
TEST_REAL_EQUAL(it->getIntensity(), 64.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 5)
TEST_REAL_EQUAL(it->getIntensity(), -32.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 4)
TEST_REAL_EQUAL(it->getIntensity(), 16.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 3)
TEST_REAL_EQUAL(it->getIntensity(), -8.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 2)
TEST_REAL_EQUAL(it->getIntensity(), 4.0)
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 1)
TEST_REAL_EQUAL(it->getIntensity(), -2.0)
RESULT

CHECK( default constructor / operator = )
DPeakList<1>::Iterator it;
it = dpl.begin();
TEST_REAL_EQUAL(it->getPosition()[0], 1)
TEST_REAL_EQUAL(it->getIntensity(), -2.0)
RESULT

CHECK( copy constructor )
DPeakList<1>::Iterator it = dpl.begin();
DPeakList<1>::Iterator it2(it);
TEST_REAL_EQUAL(it2->getPosition()[0], 1)
TEST_REAL_EQUAL(it2->getIntensity(), -2.0)
RESULT

CHECK( operator ++ (int) )
DPeakList<1>::Iterator it = dpl.begin();
it++;
TEST_REAL_EQUAL(it->getPosition()[0], 2)
DPeakList<1>::Iterator it2(it++);
TEST_REAL_EQUAL(it2->getPosition()[0], 2)
TEST_REAL_EQUAL(it->getPosition()[0], 3)
RESULT

CHECK( operator -- (int) )
DPeakList<1>::Iterator it = dpl.end();
it--;
TEST_REAL_EQUAL(it->getPosition()[0], 10)
DPeakList<1>::Iterator it2(it--);
TEST_REAL_EQUAL(it2->getPosition()[0], 10)
TEST_REAL_EQUAL(it->getPosition()[0], 9)
RESULT

CHECK( operator * () assignment )
DPeakList<1>::Iterator it = dpl.begin();
++it;
++it;
(*(it++)).getPosition()[0]=37;
TEST_REAL_EQUAL(it->getPosition()[0],4)
it = dpl.begin();
++it;
TEST_REAL_EQUAL(it->getPosition()[0],2)
++it;
TEST_REAL_EQUAL(it->getPosition()[0],37)
++it;
TEST_REAL_EQUAL(it->getPosition()[0],4)
RESULT

CHECK( operator -> () assignment )
DPeakList<1>::Iterator it = dpl.begin();
++it;
++it;
(it++)->getPosition()[0]=37;
TEST_REAL_EQUAL(it->getPosition()[0],4)
it = dpl.begin();
++it;
TEST_REAL_EQUAL(it->getPosition()[0],2)
++it;
TEST_REAL_EQUAL(it->getPosition()[0],37)
++it;
TEST_REAL_EQUAL(it->getPosition()[0],4)
RESULT

CHECK( operator -> () const )
DPeakList<1>::Iterator it = dpl.begin();
++it;
TEST_REAL_EQUAL(it->getPosition()[0],2)
RESULT

CHECK( operator == () )
DPeakList<1>::Iterator it = dpl.begin();
DPeakList<1>::Iterator it2 = dpl.begin();
TEST_REAL_EQUAL(it==it2,true)
++it2;
TEST_REAL_EQUAL(it==it2,false)
++it;
TEST_REAL_EQUAL(it==it2,true)
RESULT

CHECK( operator != () )
DPeakList<1>::Iterator it = dpl.begin();
DPeakList<1>::Iterator it2 = dpl.begin();
TEST_REAL_EQUAL(it!=it2,false)
++it2;
TEST_REAL_EQUAL(it!=it2,true)
++it;
TEST_REAL_EQUAL(it!=it2,false)
RESULT

CHECK( swap(i1,i2) )
DPeakList<1>::Iterator it = dpl.begin();
DPeakList<1>::Iterator it2 = dpl.end();
--it2;
TEST_REAL_EQUAL(it->getPosition()[0],1)
TEST_REAL_EQUAL(it2->getPosition()[0],10)
swap(it,it2);
TEST_REAL_EQUAL(it->getPosition()[0],10)
TEST_REAL_EQUAL(it2->getPosition()[0],1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

