// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
