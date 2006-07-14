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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeakArray.h>

///////////////////////////

START_TEST(DPeakArray<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

//construct a peak array to test on
DPeakArray<1> dpa;
for(unsigned int i=1;i<11;i++)
{
	DPeak<1> peak;
	peak.getPosition()[0]=i;
	peak.getIntensity()=pow((float)-2,(int)i);
	dpa.push_back(peak);
}

CHECK( operator * () / operator + (size_type) )
TEST_REAL_EQUAL((*dpa.begin()).getPosition()[0], 1)
TEST_REAL_EQUAL((*dpa.begin()).getIntensity(), -2.0)
TEST_REAL_EQUAL((*(dpa.begin()+1)).getPosition()[0], 2)
TEST_REAL_EQUAL((*(dpa.begin()+1)).getIntensity(), 4.0)
TEST_REAL_EQUAL((*(dpa.begin()+2)).getPosition()[0], 3)
TEST_REAL_EQUAL((*(dpa.begin()+2)).getIntensity(), -8.0)
TEST_REAL_EQUAL((*(dpa.begin()+3)).getPosition()[0], 4)
TEST_REAL_EQUAL((*(dpa.begin()+3)).getIntensity(), 16.0)
RESULT

CHECK( operator -> () / operator - (size_type) )
TEST_REAL_EQUAL((dpa.end()-1)->getPosition()[0], 10)
TEST_REAL_EQUAL((dpa.end()-1)->getIntensity(), 1024.0)
TEST_REAL_EQUAL((dpa.end()-2)->getPosition()[0], 9)
TEST_REAL_EQUAL((dpa.end()-2)->getIntensity(), -512.0)
TEST_REAL_EQUAL((dpa.end()-3)->getPosition()[0], 8)
TEST_REAL_EQUAL((dpa.end()-3)->getIntensity(), 256.0)
TEST_REAL_EQUAL((dpa.end()-4)->getPosition()[0], 7)
TEST_REAL_EQUAL((dpa.end()-4)->getIntensity(), -128.0)
RESULT

CHECK( default constructor / operator = )
DPeakArray<1>::Iterator it;
it = dpa.begin();
TEST_REAL_EQUAL(it->getPosition()[0], 1)
TEST_REAL_EQUAL(it->getIntensity(), -2.0)
RESULT

CHECK( copy constructor )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2(it);
TEST_REAL_EQUAL(it2->getPosition()[0], 1)
TEST_REAL_EQUAL(it2->getIntensity(), -2.0)
RESULT

CHECK( operator ++ )
DPeakArray<1>::Iterator it = dpa.begin();
++it;
TEST_REAL_EQUAL(it->getPosition()[0], 2)
DPeakArray<1>::Iterator it2(++it);
TEST_REAL_EQUAL(it2->getPosition()[0], 3)
TEST_REAL_EQUAL(it->getPosition()[0], 3)
RESULT

CHECK( operator ++ (int) )
DPeakArray<1>::Iterator it = dpa.begin();
it++;
TEST_REAL_EQUAL(it->getPosition()[0], 2)
DPeakArray<1>::Iterator it2(it++);
TEST_REAL_EQUAL(it2->getPosition()[0], 2)
TEST_REAL_EQUAL(it->getPosition()[0], 3)
RESULT

CHECK( operator -- )
DPeakArray<1>::Iterator it = dpa.end();
--it;
TEST_REAL_EQUAL(it->getPosition()[0], 10)
DPeakArray<1>::Iterator it2(--it);
TEST_REAL_EQUAL(it2->getPosition()[0], 9)
TEST_REAL_EQUAL(it->getPosition()[0], 9)
RESULT

CHECK( operator -- (int) )
DPeakArray<1>::Iterator it = dpa.end();
it--;
TEST_REAL_EQUAL(it->getPosition()[0], 10)
DPeakArray<1>::Iterator it2(it--);
TEST_REAL_EQUAL(it2->getPosition()[0], 10)
TEST_REAL_EQUAL(it->getPosition()[0], 9)
RESULT

CHECK( friend operator + (size_type , Iterator) )
DPeakArray<1>::Iterator it = 1+dpa.begin();
TEST_REAL_EQUAL(it->getPosition()[0],2)
TEST_REAL_EQUAL(it->getIntensity(), 4.0)
RESULT

CHECK( operator += (size_type) )
DPeakArray<1>::Iterator it = dpa.begin();
it +=4;
TEST_REAL_EQUAL(it->getPosition()[0],5)
TEST_REAL_EQUAL(it->getIntensity(), -32.0)
RESULT

CHECK( operator -= (size_type) )
DPeakArray<1>::Iterator it = dpa.end();
it -=6;
TEST_REAL_EQUAL(it->getPosition()[0],5)
TEST_REAL_EQUAL(it->getIntensity(), -32.0)
RESULT

CHECK( friend operator - (Iterator,Iterator) )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
TEST_REAL_EQUAL(it2-it,10)
TEST_REAL_EQUAL((it2-2)-(it+2),6)
TEST_REAL_EQUAL((it+6)-(it2-6),2)
RESULT

CHECK( operator [] (size_type) )
DPeakArray<1>::Iterator it = dpa.begin();
TEST_REAL_EQUAL(it[1].getPosition()[0],2)
TEST_REAL_EQUAL(it[5].getPosition()[0],6)
DPeakArray<1>::Iterator it2 = dpa.end();
TEST_REAL_EQUAL(it2[-1].getPosition()[0],10)
TEST_REAL_EQUAL(it2[-5].getPosition()[0],6)
RESULT

CHECK( operator [] (size_type) assignment)
DPeakArray<1>::Iterator it = dpa.begin();
it[3].getPosition()[0]=4711;
TEST_REAL_EQUAL(dpa[2].getPosition()[0],3)
TEST_REAL_EQUAL(dpa[3].getPosition()[0],4711)
TEST_REAL_EQUAL(dpa[4].getPosition()[0],5)
RESULT

CHECK( operator * () assignment )
DPeakArray<1>::Iterator it = dpa.begin();
it+=3;
(*(it++)).getPosition()[0]=45;
TEST_REAL_EQUAL(dpa[2].getPosition()[0],3)
TEST_REAL_EQUAL(dpa[3].getPosition()[0],45)
TEST_REAL_EQUAL(dpa[4].getPosition()[0],5)
RESULT

CHECK( operator -> () assignment )
DPeakArray<1>::Iterator it = dpa.begin();
it+=3;
(it++)->getPosition()[0]=47;
TEST_REAL_EQUAL(dpa[2].getPosition()[0],3)
TEST_REAL_EQUAL(dpa[3].getPosition()[0],47)
TEST_REAL_EQUAL(dpa[4].getPosition()[0],5)
RESULT

CHECK( operator -> () const )
DPeakArray<1>::Iterator it = dpa.begin();
it+=3;
TEST_REAL_EQUAL(it->getPosition()[0],47)
RESULT

CHECK( operator < () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it<it2,true)
--it2;
TEST_REAL_EQUAL(it<it2,false)
--it2;
TEST_REAL_EQUAL(it<it2,false)
RESULT

CHECK( operator > () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it>it2,false)
--it2;
TEST_REAL_EQUAL(it>it2,false)
--it2;
TEST_REAL_EQUAL(it>it2,true)
RESULT

CHECK( operator <= () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it<=it2,true)
--it2;
TEST_REAL_EQUAL(it<=it2,true)
--it2;
TEST_REAL_EQUAL(it<=it2,false)
RESULT

CHECK( operator >= () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it>=it2,false)
--it2;
TEST_REAL_EQUAL(it>=it2,true)
--it2;
TEST_REAL_EQUAL(it>=it2,true)
RESULT

CHECK( operator == () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it==it2,false)
--it2;
TEST_REAL_EQUAL(it==it2,true)
--it2;
TEST_REAL_EQUAL(it==it2,false)
RESULT

CHECK( operator != () )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end();
it+=5;
it2-=4;
TEST_REAL_EQUAL(it!=it2,true)
--it2;
TEST_REAL_EQUAL(it!=it2,false)
--it2;
TEST_REAL_EQUAL(it!=it2,true)
RESULT

CHECK( swap(i1,i2) )
DPeakArray<1>::Iterator it = dpa.begin();
DPeakArray<1>::Iterator it2 = dpa.end()-1;
TEST_REAL_EQUAL(it->getPosition()[0],1)
TEST_REAL_EQUAL(it2->getPosition()[0],10)
swap(it,it2);
TEST_REAL_EQUAL(it->getPosition()[0],10)
TEST_REAL_EQUAL(it2->getPosition()[0],1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

