// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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

InRTRange<MSSpectrum>* ptr = nullptr;
InRTRange<MSSpectrum>* nullPointer = nullptr;
START_SECTION((InRTRange(double min, double max, bool reverse = false)))
	ptr = new InRTRange<MSSpectrum>(5,10,false);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(([EXTRA]~InRTRange()))
	delete ptr;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	InRTRange<MSSpectrum> r(5,10,false);
	InRTRange<MSSpectrum> r2(5,10,true);
	MSSpectrum s;
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

InMSLevelRange<MSSpectrum>* ptr2 = nullptr;
InMSLevelRange<MSSpectrum>* nullPointer2 = nullptr;
START_SECTION((MSLevelRange(const IntList& levels, bool reverse = false)))
	IntList tmp;
	ptr2 = new InMSLevelRange<MSSpectrum>(tmp,false);
  TEST_NOT_EQUAL(ptr2, nullPointer2)
END_SECTION

START_SECTION(([EXTRA]~InMSLevelRange()))
	delete ptr2;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IntList tmp;
	tmp.push_back(2);
	tmp.push_back(3);
	tmp.push_back(4);
	InMSLevelRange<MSSpectrum> r(tmp,false);
	InMSLevelRange<MSSpectrum> r2(tmp,true);
	MSSpectrum s;
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

HasScanMode<MSSpectrum>* ptr2_1 = nullptr;
HasScanMode<MSSpectrum>* nullPointer2_1 = nullptr;
START_SECTION((HasScanMode(Int mode, bool reverse = false)))
	ptr2_1 = new HasScanMode<MSSpectrum>(1,false);
  TEST_NOT_EQUAL(ptr2_1, nullPointer2_1)
END_SECTION

START_SECTION(([EXTRA]~HasScanMode()))
	delete ptr2_1;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	HasScanMode<MSSpectrum> r(InstrumentSettings::SIM,false);
	HasScanMode<MSSpectrum> r2(InstrumentSettings::MASSSPECTRUM,true);
	MSSpectrum s;
	s.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
	TEST_EQUAL(r(s), true);
	TEST_EQUAL(r2(s), true);
	s.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
	TEST_EQUAL(r(s), false);
	TEST_EQUAL(r2(s), false);
END_SECTION

//InMzRange

InMzRange<Peak1D >* ptr3 = nullptr;
InMzRange<Peak1D >* nullPointer3 = nullptr;
START_SECTION((InMzRange(double min, double max, bool reverse = false)))
	ptr3 = new InMzRange<Peak1D >(5.0,10.0,false);
  TEST_NOT_EQUAL(ptr3, nullPointer3)
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

InIntensityRange<Peak1D >* ptr4 = nullptr;
InIntensityRange<Peak1D >* nullPointer4 = nullptr;
START_SECTION((IntensityRange(double min, double max, bool reverse = false)))
	ptr4 = new InIntensityRange<Peak1D >(5.0,10.0,false);
  TEST_NOT_EQUAL(ptr4, nullPointer4)
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

IsEmptySpectrum<MSSpectrum>* ptr47 = nullptr;
IsEmptySpectrum<MSSpectrum>* nullPointer47 = nullptr;
START_SECTION((IsEmptySpectrum(bool reverse = false)))
	ptr47 = new IsEmptySpectrum<MSSpectrum>();
  TEST_NOT_EQUAL(ptr47, nullPointer47)
END_SECTION

START_SECTION(([EXTRA]~IsEmptySpectrum()))
	delete ptr47;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IsEmptySpectrum<MSSpectrum> s;
	IsEmptySpectrum<MSSpectrum> s2(true);
	MSSpectrum spec;
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
	spec.resize(5);
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);
END_SECTION

//IsZoomSpectrum

IsZoomSpectrum<MSSpectrum>* ptr48 = nullptr;
IsZoomSpectrum<MSSpectrum>* nullPointer48 = nullptr;
START_SECTION((IsZoomSpectrum(bool reverse = false)))
	ptr48 = new IsZoomSpectrum<MSSpectrum>();
  TEST_NOT_EQUAL(ptr48, nullPointer48)
END_SECTION

START_SECTION(([EXTRA]~IsZoomSpectrum()))
	delete ptr48;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	IsZoomSpectrum<MSSpectrum> s;
	IsZoomSpectrum<MSSpectrum> s2(true);
	MSSpectrum spec;
	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);
	spec.getInstrumentSettings().setZoomScan(true);
	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);
END_SECTION

//HasActivationMethod

HasActivationMethod<MSSpectrum>* ptr49 = nullptr;
HasActivationMethod<MSSpectrum>* nullPointer49 = nullptr;
START_SECTION((HasActivationMethod(const StringList& methods, bool reverse = false)))
	ptr49 = new HasActivationMethod<MSSpectrum>(ListUtils::create<String>(""));
  TEST_NOT_EQUAL(ptr49, nullPointer49)
END_SECTION

START_SECTION(([EXTRA]~HasActivationMethod()))
	delete ptr49;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	HasActivationMethod<MSSpectrum> s(ListUtils::create<String>(Precursor::NamesOfActivationMethod[1]+","+Precursor::NamesOfActivationMethod[2]));
	HasActivationMethod<MSSpectrum> s2(ListUtils::create<String>(Precursor::NamesOfActivationMethod[1]+","+Precursor::NamesOfActivationMethod[2]),true);
	MSSpectrum spec;
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


//InPrecursorMZRange

InPrecursorMZRange<MSSpectrum>* ptr50 = nullptr;
InPrecursorMZRange<MSSpectrum>* nullPointer50 = nullptr;
START_SECTION((InPrecursorMZRange(const double& mz_left, const double& mz_right, bool reverse = false)))
	ptr50 = new InPrecursorMZRange<MSSpectrum>(100.0, 200.0);
  TEST_NOT_EQUAL(ptr50, nullPointer50)
END_SECTION

START_SECTION(([EXTRA]~InPrecursorMZRange()))
	delete ptr50;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
	InPrecursorMZRange<MSSpectrum> s(100.0, 200.0);
	InPrecursorMZRange<MSSpectrum> s2(100.0, 200.0,true);
	MSSpectrum spec;
	std::vector<Precursor> pc;
	Precursor p;
  p.setMZ(150.0);
  pc.push_back(p);
	spec.setPrecursors(pc);

	TEST_EQUAL(s(spec), true);
	TEST_EQUAL(s2(spec), false);

	// outside of allowed window
	p.setMZ(444.0);
	pc[0] = p;
	spec.setPrecursors(pc);

	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);

	// multiple precursors:
	// adding second which is within limits... but we require all of them to be...
	p.setMZ(150.0);
	pc.push_back(p);
	spec.setPrecursors(pc);

	TEST_EQUAL(s(spec), false);
	TEST_EQUAL(s2(spec), true);



END_SECTION



//InPrecursorMZRange

IsInIsolationWindow<MSSpectrum>* ptr500 = nullptr;
IsInIsolationWindow<MSSpectrum>* nullPointer500 = nullptr;
START_SECTION((IsInIsolationWindow(const double& mz_left, const double& mz_right, bool reverse = false)))
  ptr500 = new IsInIsolationWindow<MSSpectrum>(ListUtils::create<double>("100.0, 200.0"));
  TEST_NOT_EQUAL(ptr500, nullPointer500)
END_SECTION

START_SECTION(([EXTRA]~IsInIsolationWindow()))
  delete ptr500;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
  IsInIsolationWindow<MSSpectrum> s(ListUtils::create<double>("300.0, 100.0, 200.0, 400.0"));    // unsorted on purpose
  IsInIsolationWindow<MSSpectrum> s2(ListUtils::create<double>("300.0, 100.0, 200.0, 400.0"), true);
  
  MSSpectrum spec;
  spec.setMSLevel(2);
  std::vector<Precursor> pc;
  Precursor p;
  p.setMZ(200.3);
  p.setIsolationWindowLowerOffset(0.5);
  p.setIsolationWindowUpperOffset(0.5);
  pc.push_back(p);
  spec.setPrecursors(pc);

  TEST_EQUAL(s(spec), true);
  TEST_EQUAL(s2(spec), false);

  // outside of allowed window
  p.setMZ(201.1);
  pc[0] = p;
  spec.setPrecursors(pc);

  TEST_EQUAL(s(spec), false);
  TEST_EQUAL(s2(spec), true);

  // multiple precursors:
  // adding second which is within limits... so it's a hit (any PC must match)
  p.setMZ(299.9);
  pc.push_back(p);
  spec.setPrecursors(pc);

  TEST_EQUAL(s(spec), true);
  TEST_EQUAL(s2(spec), false);

END_SECTION


//HasScanPolarity

HasScanPolarity<MSSpectrum>* ptr51 = nullptr;
HasScanPolarity<MSSpectrum>* nullPointer51 = nullptr;
START_SECTION((HasScanPolarity(Int polarity,bool reverse = false)))
  ptr51 = new HasScanPolarity<MSSpectrum>(0);
  TEST_NOT_EQUAL(ptr48, nullPointer51)
END_SECTION

START_SECTION(([EXTRA]~HasScanPolarity()))
  delete ptr51;
END_SECTION

START_SECTION((bool operator()(const SpectrumType& s) const))
  HasScanPolarity<MSSpectrum> s(IonSource::POSITIVE);
  HasScanPolarity<MSSpectrum> s2(IonSource::POSITIVE, true);
  MSSpectrum spec;
  TEST_EQUAL(s(spec), false);
  TEST_EQUAL(s2(spec), true);
  spec.getInstrumentSettings().setPolarity(IonSource::POSITIVE);
  TEST_EQUAL(s(spec), true);
  TEST_EQUAL(s2(spec), false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
