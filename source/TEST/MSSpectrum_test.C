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
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum<>* ptr = 0;
CHECK((MSSpectrum()))
	ptr = new MSSpectrum<>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MSSpectrum()))
	delete ptr;
RESULT

CHECK((MSSpectrum(const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().setMzRangeStart(5.1);
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	tmp.getContainer().push_back(peak);
	
	MSSpectrum<> tmp2(tmp);
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings().getMzRangeStart(),5.1);
	TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_EQUAL(tmp2.getContainer()[0].getPosition()[0],47.11);
RESULT

CHECK((MSSpectrum& operator= (const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().setMzRangeStart(5.1);
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	tmp.getContainer().push_back(peak);
	
	//normal assignment
	MSSpectrum<> tmp2;
	tmp2 = tmp;
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings().getMzRangeStart(),5.1);
	TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_EQUAL(tmp2.getContainer()[0].getPosition()[0],47.11);
	
	//Assignment of empty object
	//normal assignment
	tmp2 = MSSpectrum<>();
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings().getMzRangeStart(),0.0);
	TEST_EQUAL(tmp2.size(),0);
RESULT

CHECK((bool operator== (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.getInstrumentSettings().setMzRangeStart(5.1);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	edit.getContainer().push_back(peak);
	TEST_EQUAL(edit==empty,false);
RESULT

CHECK((bool operator!= (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.getInstrumentSettings().setMzRangeStart(5.1);
	TEST_EQUAL(edit!=empty,true);
	
	edit = empty;
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	edit.getContainer().push_back(peak);
	TEST_EQUAL(edit!=empty,true);
RESULT

CHECK(([EXTRA] MSSpectrum<DRawDataPoint<1> >))
	MSSpectrum<DRawDataPoint<1> > tmp;
	MSSpectrum<DRawDataPoint<1> >::PeakType rdp;
	rdp.getPosition()[0] = 47.11;
	tmp.getContainer().push_back(rdp);
	TEST_EQUAL(tmp.size(),1);
	TEST_REAL_EQUAL(tmp.getContainer()[0].getPosition()[0],47.11);	
RESULT

CHECK((typename BaseSpectrum::Iterator MZEnd(double mz)))
	MSSpectrum<DRawDataPoint<1> > tmp;
	MSSpectrum<DRawDataPoint<1> >::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	MSSpectrum<DRawDataPoint<1> >::Iterator it;
	
	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((typename BaseSpectrum::Iterator MZBegin(double mz)))
	MSSpectrum<DRawDataPoint<1> > tmp;
	MSSpectrum<DRawDataPoint<1> >::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	MSSpectrum<DRawDataPoint<1> >::Iterator it;
	
	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((const typename BaseSpectrum::ConstIterator MZEnd(double mz) const))
	MSSpectrum<DRawDataPoint<1> > tmp;
	MSSpectrum<DRawDataPoint<1> >::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	MSSpectrum<DRawDataPoint<1> >::ConstIterator it;
	
	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((const typename BaseSpectrum::ConstIterator MZBegin(double mz) const))
	MSSpectrum<DRawDataPoint<1> > tmp;
	MSSpectrum<DRawDataPoint<1> >::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	MSSpectrum<DRawDataPoint<1> >::ConstIterator it;
	
	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



