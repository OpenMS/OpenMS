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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

///////////////////////////

START_TEST(RangeUtils<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

//InRTRange

InRTRange<DSpectrum<> >* ptr = 0;
START_SECTION((InRTRange(double min, double max, bool reverse = false)))
	ptr = new InRTRange<DSpectrum<> >(5,10,false);
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(([EXTRA]~InRTRange()))
	delete ptr;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	InRTRange<DSpectrum<> > r(5,10,false);
	InRTRange<DSpectrum<> > r2(5,10,true);
	DSpectrum<> s;
	s.setRT(4.9);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);
	s.setRT(5.0);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRT(7.5);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRT(10.0);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRT(10.1);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);	
END_SECTION


//MSLevelRange

InMSLevelRange<DSpectrum<> >* ptr2 = 0;
START_SECTION((MSLevelRange(const IntList& levels, bool reverse = false)))
	IntList tmp;
	ptr2 = new InMSLevelRange<DSpectrum<> >(tmp,false);
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION(([EXTRA]~InMSLevelRange()))
	delete ptr2;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	vector<UInt> tmp;
	tmp.push_back(2);
	tmp.push_back(3);
	tmp.push_back(4);
	InMSLevelRange<DSpectrum<> > r(tmp,false);
	InMSLevelRange<DSpectrum<> > r2(tmp,true);
	DSpectrum<> s;
	s.setMSLevel(1);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);
	s.setMSLevel(2);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setMSLevel(3);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setMSLevel(4);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setMSLevel(5);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);	
END_SECTION

//HasScanMode

HasScanMode<MSSpectrum<> >* ptr2_1 = 0;
START_SECTION((HasScanMode(Int mode, bool reverse = false)))
	ptr2_1 = new HasScanMode<MSSpectrum<> >(1,false);
	TEST_NOT_EQUAL(ptr2, 0)
END_SECTION

START_SECTION(([EXTRA]~HasScanMode()))
	delete ptr2_1;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	HasScanMode<MSSpectrum<> > r(InstrumentSettings::ZOOM,false);
	HasScanMode<MSSpectrum<> > r2(InstrumentSettings::FULL,true);
	MSSpectrum<> s;
	s.getInstrumentSettings().setScanMode(InstrumentSettings::ZOOM);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), true);
	s.getInstrumentSettings().setScanMode(InstrumentSettings::FULL);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), false);
END_SECTION

//InMzRange

InMzRange<Peak1D >* ptr3 = 0;
START_SECTION((InMzRange(double min, double max, bool reverse = false)))
	ptr3 = new InMzRange<Peak1D >(5.0,10.0,false);
	TEST_NOT_EQUAL(ptr3, 0)
END_SECTION

START_SECTION(([EXTRA]~InMzRange()))
	delete ptr3;
END_SECTION

START_SECTION((bool operator()(const PeakType& p) const))
	InMzRange<Peak1D > r(5.0,10.0,false);
	InMzRange<Peak1D > r2(5.0,10.0,true);
	Peak1D p;
	p.getPosition()[0] = 4.9;
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);
	p.getPosition()[0] = 5.0;
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.getPosition()[0] = 7.5;
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.getPosition()[0] = 10.0;
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.getPosition()[0] = 10.1;
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);	
END_SECTION

//IntensityRange

InIntensityRange<Peak1D >* ptr4 = 0;
START_SECTION((IntensityRange(double min, double max, bool reverse = false)))
	ptr4 = new InIntensityRange<Peak1D >(5.0,10.0,false);
	TEST_NOT_EQUAL(ptr4, 0)
END_SECTION

START_SECTION(([EXTRA]~InIntensityRange()))
	delete ptr4;
END_SECTION

START_SECTION((bool operator()(const PeakType& p) const))
	InIntensityRange<Peak1D > r(5.0,10.0,false);
	InIntensityRange<Peak1D > r2(5.0,10.0,true);
	Peak1D p;
	p.setIntensity(4.9);
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);
	p.setIntensity(5.0);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(7.5);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(10.0);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(10.1);
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);	
END_SECTION


//IsEmptySpectrum

IsEmptySpectrum<DSpectrum<> >* ptr47 = 0;
START_SECTION((IsEmptySpectrum(bool reverse = false)))
	ptr47 = new IsEmptySpectrum<DSpectrum<> >();
	TEST_NOT_EQUAL(ptr47, 0)
END_SECTION

START_SECTION(([EXTRA]~IsEmptySpectrum()))
	delete ptr47;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IsEmptySpectrum<DSpectrum<> > s;
	IsEmptySpectrum<DSpectrum<> > s2(true);
	DSpectrum<> spec;
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
	spec.resize(5);
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
