// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
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

InRTRange<MSSpectrum<> >* ptr = 0;
START_SECTION((InRTRange(double min, double max, bool reverse = false)))
	ptr = new InRTRange<MSSpectrum<> >(5,10,false);
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(([EXTRA]~InRTRange()))
	delete ptr;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	InRTRange<MSSpectrum<> > r(5,10,false);
	InRTRange<MSSpectrum<> > r2(5,10,true);
	MSSpectrum<> s;
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

InMSLevelRange<MSSpectrum<> >* ptr2 = 0;
START_SECTION((MSLevelRange(const IntList& levels, bool reverse = false)))
	IntList tmp;
	ptr2 = new InMSLevelRange<MSSpectrum<> >(tmp,false);
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
	InMSLevelRange<MSSpectrum<> > r(tmp,false);
	InMSLevelRange<MSSpectrum<> > r2(tmp,true);
	MSSpectrum<> s;
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
	HasScanMode<MSSpectrum<> > r(InstrumentSettings::SIM,false);
	HasScanMode<MSSpectrum<> > r2(InstrumentSettings::MASSSPECTRUM,true);
	MSSpectrum<> s;
	s.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), true);
	s.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
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
	p.setIntensity(4.9f);
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);
	p.setIntensity(5.0f);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(7.5f);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(10.0f);
	TEST_EQUAL(r(p), true);
	TEST_EQUAL(r2(p), false);
	p.setIntensity(10.1f);
	TEST_EQUAL(r(p), false);
	TEST_EQUAL(r2(p), true);	
END_SECTION


//IsEmptySpectrum

IsEmptySpectrum<MSSpectrum<> >* ptr47 = 0;
START_SECTION((IsEmptySpectrum(bool reverse = false)))
	ptr47 = new IsEmptySpectrum<MSSpectrum<> >();
	TEST_NOT_EQUAL(ptr47, 0)
END_SECTION

START_SECTION(([EXTRA]~IsEmptySpectrum()))
	delete ptr47;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IsEmptySpectrum<MSSpectrum<> > s;
	IsEmptySpectrum<MSSpectrum<> > s2(true);
	MSSpectrum<> spec;
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
	spec.resize(5);
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);	
END_SECTION

//IsZoomSpectrum

IsZoomSpectrum<MSSpectrum<> >* ptr48 = 0;
START_SECTION((IsZoomSpectrum(bool reverse = false)))
	ptr48 = new IsZoomSpectrum<MSSpectrum<> >();
	TEST_NOT_EQUAL(ptr48, 0)
END_SECTION

START_SECTION(([EXTRA]~IsZoomSpectrum()))
	delete ptr48;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IsZoomSpectrum<MSSpectrum<> > s;
	IsZoomSpectrum<MSSpectrum<> > s2(true);
	MSSpectrum<> spec;
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);
	spec.getInstrumentSettings().setZoomScan(true);
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);	
END_SECTION

//HasActivationMethod

HasActivationMethod<MSSpectrum<> >* ptr49 = 0;
START_SECTION((HasActivationMethod(const StringList& methods, bool reverse = false)))
	ptr49 = new HasActivationMethod<MSSpectrum<> >(StringList::create(""));
	TEST_NOT_EQUAL(ptr49, 0)
END_SECTION

START_SECTION(([EXTRA]~HasActivationMethod()))
	delete ptr49;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	HasActivationMethod<MSSpectrum<> > s(StringList::create(Precursor::NamesOfActivationMethod[1]+","+Precursor::NamesOfActivationMethod[2]));
	HasActivationMethod<MSSpectrum<> > s2(StringList::create(Precursor::NamesOfActivationMethod[1]+","+Precursor::NamesOfActivationMethod[2]),true);
	MSSpectrum<> spec;
	std::vector<Precursor> pc;
	Precursor p;
	set <Precursor::ActivationMethod> sa1;
	sa1.insert( Precursor::PSD ); //occurs
	sa1.insert( Precursor::BIRD );//just a dummy

	p.setActivationMethods(sa1);
	pc.push_back(p);
	spec.setPrecursors(pc);

	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
	
	// does not occur as activation method
	set <Precursor::ActivationMethod> sa2;
	sa2.insert( Precursor::BIRD );
	p.setActivationMethods(sa2);
	pc[0] = p;
	spec.setPrecursors(pc);

	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);	

	// multiple precursors:
	// adding another dummy
	set <Precursor::ActivationMethod> sa3;
	sa3.insert( Precursor::LCID );
	p.setActivationMethods(sa3);
	pc.push_back(p);
	spec.setPrecursors(pc);
	
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);	
	
	// adding a matching precursor
	set <Precursor::ActivationMethod> sa4;
	sa4.insert( Precursor::PD );
	p.setActivationMethods(sa4);
	pc.push_back(p);
	spec.setPrecursors(pc);
	
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);	
	

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
