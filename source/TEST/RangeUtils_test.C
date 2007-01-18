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

//RTRange

RTRange<DSpectrum<1> >* ptr = 0;
CHECK((RTRange(double min, double max, bool reverse = false)))
	ptr = new RTRange<DSpectrum<1> >(5,10,false);
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(([EXTRA]~RTRange()))
	delete ptr;
RESULT

CHECK((bool operator()(const SpectrumType& s) const))
	RTRange<DSpectrum<1> > r(5,10,false);
	RTRange<DSpectrum<1> > r2(5,10,true);
	DSpectrum<1> s;
	s.setRetentionTime(4.9);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);
	s.setRetentionTime(5.0);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRetentionTime(7.5);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRetentionTime(10.0);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), false);
	s.setRetentionTime(10.1);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), true);	
RESULT


//MSLevelRange

MSLevelRange<DSpectrum<1> >* ptr2 = 0;
CHECK((MSLevelRange(const std::vector<UnsignedInt>& levels, bool reverse = false)))
	vector<UnsignedInt> tmp;
	ptr2 = new MSLevelRange<DSpectrum<1> >(tmp,false);
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK(([EXTRA]~MSLevelRange()))
	delete ptr2;
RESULT

CHECK((bool operator()(const SpectrumType& s) const))
	vector<UnsignedInt> tmp;
	tmp.push_back(2);
	tmp.push_back(3);
	tmp.push_back(4);
	MSLevelRange<DSpectrum<1> > r(tmp,false);
	MSLevelRange<DSpectrum<1> > r2(tmp,true);
	DSpectrum<1> s;
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
RESULT

//ScanModePredicate

ScanModePredicate<MSSpectrum<> >* ptr2_1 = 0;
CHECK((ScanModePredicate(SignedInt mode, bool reverse = false)))
	ptr2_1 = new ScanModePredicate<MSSpectrum<> >(1,false);
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK(([EXTRA]~ScanModePredicate()))
	delete ptr2_1;
RESULT

CHECK((bool operator()(const SpectrumType& s) const))
	ScanModePredicate<MSSpectrum<> > r(1,false);
	ScanModePredicate<MSSpectrum<> > r2(2,true);
	MSSpectrum<> s;
	s.getInstrumentSettings().setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), true);
	s.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSCAN);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), false);
RESULT

//MZRange

MzRange<DPeak<1> >* ptr3 = 0;
CHECK((MzRange(double min, double max, bool reverse = false)))
	ptr3 = new MzRange<DPeak<1> >(5.0,10.0,false);
	TEST_NOT_EQUAL(ptr3, 0)
RESULT

CHECK(([EXTRA]~MzRange()))
	delete ptr3;
RESULT

CHECK((bool operator()(const PeakType& p) const))
	MzRange<DPeak<1> > r(5.0,10.0,false);
	MzRange<DPeak<1> > r2(5.0,10.0,true);
	DPeak<1> p;
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
RESULT

//IntensityRange

IntensityRange<DPeak<1> >* ptr4 = 0;
CHECK((IntensityRange(double min, double max, bool reverse = false)))
	ptr4 = new IntensityRange<DPeak<1> >(5.0,10.0,false);
	TEST_NOT_EQUAL(ptr4, 0)
RESULT

CHECK(([EXTRA]~IntensityRange()))
	delete ptr4;
RESULT

CHECK((bool operator()(const PeakType& p) const))
	IntensityRange<DPeak<1> > r(5.0,10.0,false);
	IntensityRange<DPeak<1> > r2(5.0,10.0,true);
	DPeak<1> p;
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
RESULT


//SpectrumEmptyPredicate

SpectrumEmptyPredicate<DSpectrum<1> >* ptr47 = 0;
CHECK((SpectrumEmptyPredicate(bool reverse = false)))
	ptr47 = new SpectrumEmptyPredicate<DSpectrum<1> >();
	TEST_NOT_EQUAL(ptr47, 0)
RESULT

CHECK(([EXTRA]~SpectrumEmptyPredicate()))
	delete ptr47;
RESULT

CHECK((bool operator()(const SpectrumType& s) const))
	SpectrumEmptyPredicate<DSpectrum<1> > s;
	SpectrumEmptyPredicate<DSpectrum<1> > s2(true);
	DSpectrum<1> spec;
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
	spec.resize(5);
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
