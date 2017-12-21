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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/TXTToMzMLConverter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TXTToMzMLConverter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TXTToMzMLConverter* ptr = 0;
TXTToMzMLConverter* null_ptr = 0;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.txt");
const String output_filepath = OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.mzML");

START_SECTION(TXTToMzMLConverter())
{
  ptr = new TXTToMzMLConverter();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TXTToMzMLConverter())
{
  delete ptr;
}
END_SECTION

ptr = new TXTToMzMLConverter();

START_SECTION(MSExperiment loadInputFile(const String& filename) const)
{
  const MSExperiment experiment = ptr->loadInputFile(input_filepath);
  const vector<MSChromatogram> chromatograms = experiment.getChromatograms();
  TEST_EQUAL(chromatograms.size(), 1);
  TEST_EQUAL(chromatograms[0].size(), 3301);
  const MSChromatogram& c = chromatograms[0];
  TEST_REAL_SIMILAR(c[0].getRT(), 0.0)
  TEST_REAL_SIMILAR(c[0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(c[660].getRT(), 2.2)
  TEST_REAL_SIMILAR(c[660].getIntensity(), -0.812998)
  TEST_REAL_SIMILAR(c[1320].getRT(), 4.4)
  TEST_REAL_SIMILAR(c[1320].getIntensity(), -0.791189)
  TEST_REAL_SIMILAR(c[1980].getRT(), 6.6)
  TEST_REAL_SIMILAR(c[1980].getIntensity(), -0.285533)
  TEST_REAL_SIMILAR(c[2640].getRT(), 8.8)
  TEST_REAL_SIMILAR(c[2640].getIntensity(), -0.485941)
  TEST_REAL_SIMILAR(c[3300].getRT(), 11.0)
  TEST_REAL_SIMILAR(c[3300].getIntensity(), -0.130904)
}
END_SECTION

START_SECTION(void storeMzMLFile(const String& filename, const MSExperiment& experiment) const)
{
  const MSExperiment experiment = ptr->loadInputFile(input_filepath);
  ptr->storeMzMLFile(output_filepath, experiment);
  MzMLFile mzml;
  MSExperiment read_exp;
  mzml.load(output_filepath, read_exp);
  const MSChromatogram c1 = experiment.getChromatograms()[0];
  const MSChromatogram c2 = read_exp.getChromatograms()[0];
  TEST_EQUAL(c1.size(), c2.size())
  for (Size i = 0; i < c1.size(); ++i)
  {
    TEST_REAL_SIMILAR(c1[i].getRT(), c2[i].getRT())
    TEST_REAL_SIMILAR(c1[i].getIntensity(), c2[i].getIntensity())
  }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
