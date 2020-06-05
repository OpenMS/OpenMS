// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSPGenericFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPGenericFile* ptr = nullptr;
MSPGenericFile* null_ptr = nullptr;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("MSPGenericFile_input.msp");

START_SECTION(MSPGenericFile())
{
  ptr = new MSPGenericFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSPGenericFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  MSPGenericFile msp;
  MSExperiment experiment;
  msp.load(input_filepath, experiment);
  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  TEST_EQUAL(spectra.size(), 3)

  const MSSpectrum& s1 = spectra[0];
  TEST_EQUAL(s1.size(), 14)
  TEST_EQUAL(s1.getName(), "name1 of first")

  TEST_EQUAL(s1.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Synon"), "name2 of 1st|name3 of firsttt")

  TEST_EQUAL(s1.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Formula"), "A11B22C333")

  TEST_EQUAL(s1.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("MW"), "156")

  TEST_EQUAL(s1.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("CAS#"), "0123-45-6")

  TEST_EQUAL(s1.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("NIST#"), "654321")

  TEST_EQUAL(s1.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("DB#"), "1")

  TEST_EQUAL(s1.metaValueExists("Comments"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Comments"), "Some comment")

  TEST_EQUAL(s1.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s1.getMetaValue("Num Peaks"), "14")

  TEST_EQUAL(s1[0].getPos(), 27)
  TEST_EQUAL(s1[0].getIntensity(), 29)
  TEST_EQUAL(s1[5].getPos(), 60)
  TEST_EQUAL(s1[5].getIntensity(), 41)
  TEST_EQUAL(s1[10].getPos(), 90)
  TEST_EQUAL(s1[10].getIntensity(), 168)
  TEST_EQUAL(s1[13].getPos(), 105)
  TEST_EQUAL(s1[13].getIntensity(), 36)

  const MSSpectrum& s2 = spectra[1];
  TEST_EQUAL(s2.size(), 15)
  TEST_EQUAL(s2.getName(), "name1 of second")

  TEST_EQUAL(s2.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Synon"), "name2 of 2nd|name3 of seconddd")

  TEST_EQUAL(s2.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Formula"), "A44B55C666")

  TEST_EQUAL(s2.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("MW"), "589")

  TEST_EQUAL(s2.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("CAS#"), "3210-45-6")

  TEST_EQUAL(s2.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("NIST#"), "789564")

  TEST_EQUAL(s2.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("DB#"), "2")

  TEST_EQUAL(s2.metaValueExists("Comments"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Comments"), "Some other comment")

  TEST_EQUAL(s2.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s2.getMetaValue("Num Peaks"), "15")

  TEST_EQUAL(s2[0].getPos(), 27)
  TEST_EQUAL(s2[0].getIntensity(), 29)
  TEST_EQUAL(s2[5].getPos(), 260)
  TEST_EQUAL(s2[5].getIntensity(), 41)
  TEST_EQUAL(s2[10].getPos(), 290)
  TEST_EQUAL(s2[10].getIntensity(), 168)
  TEST_EQUAL(s2[14].getPos(), 310)
  TEST_EQUAL(s2[14].getIntensity(), 20)

  const MSSpectrum& s3 = spectra[2];
  TEST_EQUAL(s3.size(), 16)
  TEST_EQUAL(s3.getName(), "name1 of third")

  TEST_EQUAL(s3.metaValueExists("Synon"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Synon"), "name2 of 3rd|name3 of thirddd")

  TEST_EQUAL(s3.metaValueExists("Formula"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Formula"), "A12B12C123")

  TEST_EQUAL(s3.metaValueExists("MW"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("MW"), "562")

  TEST_EQUAL(s3.metaValueExists("CAS#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("CAS#"), "4210-47-4")

  TEST_EQUAL(s3.metaValueExists("NIST#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("NIST#"), "749514")

  TEST_EQUAL(s3.metaValueExists("DB#"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("DB#"), "3")

  TEST_EQUAL(s3.metaValueExists("Comments"), false) // this spectrum doesn't have a comment

  TEST_EQUAL(s3.metaValueExists("Num Peaks"), true)
  TEST_STRING_EQUAL(s3.getMetaValue("Num Peaks"), "16")

  TEST_EQUAL(s3[0].getPos(), 27)
  TEST_EQUAL(s3[0].getIntensity(), 29)
  TEST_EQUAL(s3[5].getPos(), 260)
  TEST_EQUAL(s3[5].getIntensity(), 41)
  TEST_EQUAL(s3[10].getPos(), 290)
  TEST_EQUAL(s3[10].getIntensity(), 168)
  TEST_EQUAL(s3[14].getPos(), 310)
  TEST_EQUAL(s3[14].getIntensity(), 20)
  TEST_EQUAL(s3[15].getPos(), 111)
  TEST_EQUAL(s3[15].getIntensity(), 44)
}
END_SECTION

START_SECTION(void addSpectrumToLibrary(
  MSSpectrum& spectrum,
  MSExperiment& library
))
{
  MSPGenericFile_friend msp_f;
  MSExperiment lib;

  MSSpectrum spec;
  spec.setName(""); // empty name
  spec.setMetaValue("is_valid", 1);

  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setName("foo"); // Num Peaks still absent!
  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setMetaValue("Num Peaks", "2");
  // Num Peaks is set but raw data poins have not been added
  TEST_EXCEPTION(Exception::ParseError, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.push_back(Peak1D(1.0, 2.0));
  spec.push_back(Peak1D(3.0, 4.0)); // now the spectrum is valid
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 1)

  spec.setName("bar");
  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 1);
  msp_f.addSpectrumToLibrary(spec, lib); // duplicate, won't be added
  TEST_EQUAL(lib.size(), 2)

  spec.setMetaValue("is_valid", 0);
  spec.setName("not a duplicate");
  msp_f.addSpectrumToLibrary(spec, lib);
  TEST_EQUAL(lib.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
