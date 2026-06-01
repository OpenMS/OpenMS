// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/SpectrumLookup.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumLookup, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumLookup* ptr = nullptr;
SpectrumLookup* null_ptr = nullptr;

START_SECTION((SpectrumLookup()))
{
	ptr = new SpectrumLookup();
	TEST_NOT_EQUAL(ptr, null_ptr);
  TEST_REAL_SIMILAR(ptr->rt_tolerance, 0.01);
}
END_SECTION

START_SECTION((~SpectrumLookup()))
{
	delete ptr;
}
END_SECTION

vector<MSSpectrum> spectra;
MSSpectrum spectrum;
spectrum.setNativeID("spectrum=0");
spectrum.setRT(1.0);
spectra.push_back(spectrum);
spectrum.setNativeID("spectrum=1");
spectrum.setRT(2.0);
spectra.push_back(spectrum);
spectrum.setNativeID("spectrum=2");
spectrum.setRT(3.0);
spectra.push_back(spectrum);

SpectrumLookup lookup;

START_SECTION((bool empty() const))
{
  TEST_EQUAL(lookup.empty(), true);
}
END_SECTION


START_SECTION((template <typename SpectrumContainer> void readSpectra(const SpectrumContainer&, const String&)))
{
  lookup.readSpectra(spectra);
  TEST_EQUAL(lookup.empty(), false);
}
END_SECTION


START_SECTION((Size findByRT(double) const))
{
  TEST_EQUAL(lookup.findByRT(2.0), 1);

  TEST_EXCEPTION(Exception::ElementNotFound, lookup.findByRT(5.0));
}
END_SECTION


START_SECTION((Size findByNativeID(const String&) const))
{
  TEST_EQUAL(lookup.findByNativeID("spectrum=1"), 1);

  TEST_EXCEPTION(Exception::ElementNotFound,
                 lookup.findByNativeID("spectrum=3"));
}
END_SECTION


START_SECTION((Size findByIndex(Size, bool) const))
{
  TEST_EQUAL(lookup.findByIndex(1), 1);
  TEST_EQUAL(lookup.findByIndex(1, true), 0);

  TEST_EXCEPTION(Exception::ElementNotFound, lookup.findByIndex(0, true));
}
END_SECTION


START_SECTION((Size findByScanNumber(Size) const))
{
  TEST_EQUAL(lookup.findByScanNumber(1), 1);

  TEST_EXCEPTION(Exception::ElementNotFound, lookup.findByScanNumber(5));
}
END_SECTION


START_SECTION((void addReferenceFormat(const String&)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, lookup.addReferenceFormat("XXX"));

  // tested with other methods below:
  lookup.addReferenceFormat("scan_number=(?<SCAN>\\d+)");
  lookup.addReferenceFormat("(?<ID>spectrum=\\d+)");
}
END_SECTION


START_SECTION((Size findByReference(const String&) const))
{
  TEST_EQUAL(lookup.findByReference("scan_number=1"), 1);
  TEST_EQUAL(lookup.findByReference("name=bla,spectrum=0"), 0);

  TEST_EXCEPTION(Exception::ParseError, lookup.findByReference("test123"));
}
END_SECTION


START_SECTION((static Int extractScanNumber(const String&,
                                            const boost::regex&)))
{
  boost::regex re("spectrum=(?<SCAN>\\d+)");
  TEST_EQUAL(SpectrumLookup::extractScanNumber("spectrum=42", re), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", re, true), -1);

  TEST_EXCEPTION(Exception::ParseError, SpectrumLookup::extractScanNumber("scan=42", re));
}
END_SECTION

START_SECTION((static Int extractScanNumber(const String&,
                                            const String&)))
{
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", "MS:1000768"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", "MS:1000769"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", "MS:1000771"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", "MS:1000772"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scan=42", "MS:1000776"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("sample=1 period=1 cycle=42 experiment=1", "MS:1000770"), 42001);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("file=42", "MS:1000773"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("file=42", "MS:1000775"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("index=42", "MS:1000774"), 43);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("scanId=42", "MS:1001508"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("spectrum=42", "MS:1000777"), 42);
  TEST_EQUAL(SpectrumLookup::extractScanNumber("42", "MS:1001530"), 42);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
