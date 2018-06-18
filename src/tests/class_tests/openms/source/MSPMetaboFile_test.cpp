// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/FORMAT/MSPMetaboFile.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSPMetaboFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPMetaboFile* ptr = nullptr;
MSPMetaboFile* null_ptr = nullptr;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("MSPMetaboFile_input.msp");

START_SECTION(MSPMetaboFile())
{
  ptr = new MSPMetaboFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSPMetaboFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  MSPMetaboFile msp;
  MSExperiment experiment;
  msp.load(input_filepath, experiment);
  const vector<MSSpectrum>& spectra = experiment.getSpectra();
  TEST_EQUAL(spectra.size(), 3)

  const MSSpectrum& s1 = spectra[0];
  TEST_EQUAL(s1.size(), 14)
  TEST_EQUAL(s1.getName(), "name1 of first")

  const MSSpectrum::StringDataArrays& SDAs1 = s1.getStringDataArrays();
  MSSpectrum::StringDataArrays::const_iterator it;

  it = getDataArrayByName(SDAs1, "Synon");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 2)
  TEST_STRING_EQUAL((*it)[0], "name2 of 1st")
  TEST_STRING_EQUAL((*it)[1], "name3 of firsttt")

  it = getDataArrayByName(SDAs1, "Formula");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "A11B22C333")

  it = getDataArrayByName(SDAs1, "MW");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "156")

  it = getDataArrayByName(SDAs1, "CAS#");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "0123-45-6")

  it = getDataArrayByName(SDAs1, "NIST#");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "654321")

  it = getDataArrayByName(SDAs1, "DB#");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "1")

  it = getDataArrayByName(SDAs1, "Comments");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "Some comment")

  it = getDataArrayByName(SDAs1, "Num Peaks");
  TEST_EQUAL(it == SDAs1.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "14")

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

  const MSSpectrum::StringDataArrays& SDAs2 = s2.getStringDataArrays();

  it = getDataArrayByName(SDAs2, "Synon");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 2)
  TEST_STRING_EQUAL((*it)[0], "name2 of 2nd")
  TEST_STRING_EQUAL((*it)[1], "name3 of seconddd")

  it = getDataArrayByName(SDAs2, "Formula");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "A44B55C666")

  it = getDataArrayByName(SDAs2, "MW");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "589")

  it = getDataArrayByName(SDAs2, "CAS#");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "3210-45-6")

  it = getDataArrayByName(SDAs2, "NIST#");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "789564")

  it = getDataArrayByName(SDAs2, "DB#");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "2")

  it = getDataArrayByName(SDAs2, "Comments");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "Some other comment")

  it = getDataArrayByName(SDAs2, "Num Peaks");
  TEST_EQUAL(it == SDAs2.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "15")

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

  const MSSpectrum::StringDataArrays& SDAs3 = s3.getStringDataArrays();

  it = getDataArrayByName(SDAs3, "Synon");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 2)
  TEST_STRING_EQUAL((*it)[0], "name2 of 3rd")
  TEST_STRING_EQUAL((*it)[1], "name3 of thirddd")

  it = getDataArrayByName(SDAs3, "Formula");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "A12B12C123")

  it = getDataArrayByName(SDAs3, "MW");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "562")

  it = getDataArrayByName(SDAs3, "CAS#");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "4210-47-4")

  it = getDataArrayByName(SDAs3, "NIST#");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "749514")

  it = getDataArrayByName(SDAs3, "DB#");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "3")

  it = getDataArrayByName(SDAs3, "Comments");
  TEST_EQUAL(it == SDAs3.cend(), true) // this spectrum doesn't have a comment

  it = getDataArrayByName(SDAs3, "Num Peaks");
  TEST_EQUAL(it == SDAs3.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL(it->front(), "16")

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

START_SECTION(void pushParsedInfoToNamedDataArray(
  MSSpectrum& spectrum,
  const String& name,
  const String& info
) const)
{
  MSPMetaboFile_friend msp_f;
  MSSpectrum spectrum;
  MSSpectrum::StringDataArrays::const_iterator it;

  const String field_synon { "Synon" };
  const String synon1 { "foo" };
  const String synon2 { "bar" };

  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_synon, synon1);

  const MSSpectrum::StringDataArrays& SDAs = spectrum.getStringDataArrays();

  TEST_EQUAL(SDAs.size(), 1)
  it = getDataArrayByName(SDAs, field_synon);
  TEST_EQUAL(it == SDAs.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL((*it)[0], synon1)

  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_synon, synon2);

  TEST_EQUAL(SDAs.size(), 1)
  it = getDataArrayByName(SDAs, field_synon);
  TEST_EQUAL(it == SDAs.cend(), false)
  TEST_EQUAL(it->size(), 2)
  TEST_STRING_EQUAL((*it)[0], synon1)
  TEST_STRING_EQUAL((*it)[1], synon2)

  const String field_comments { "Comments" };
  const String comment { "seems to work fine" };

  msp_f.pushParsedInfoToNamedDataArray(spectrum, field_comments, comment);

  TEST_EQUAL(SDAs.size(), 2)
  it = getDataArrayByName(SDAs, field_comments);
  TEST_EQUAL(it == SDAs.cend(), false)
  TEST_EQUAL(it->size(), 1)
  TEST_STRING_EQUAL((*it)[0], comment)
}
END_SECTION

START_SECTION(void addSpectrumToLibrary(
  MSSpectrum& spectrum,
  MSExperiment& library
))
{
  MSPMetaboFile_friend msp_f;
  MSExperiment lib;

  MSSpectrum spec;
  spec.setName(""); // empty name
  spec.setMetaValue("is_valid", 1);

  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  spec.setName("foo"); // Num Peaks still absent!
  TEST_EXCEPTION(Exception::MissingInformation, msp_f.addSpectrumToLibrary(spec, lib))
  TEST_EQUAL(lib.size(), 0)

  msp_f.pushParsedInfoToNamedDataArray(spec, "Num Peaks", "2");
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
