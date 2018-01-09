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
#include <OpenMS/FORMAT/ChromeleonFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChromeleonFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromeleonFile* ptr = 0;
ChromeleonFile* null_ptr = 0;
const String input_filepath = OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.txt");
const String output_filepath = OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.mzML");

START_SECTION(ChromeleonFile())
{
  ptr = new ChromeleonFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ChromeleonFile())
{
  delete ptr;
}
END_SECTION

ptr = new ChromeleonFile();

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  MSExperiment experiment;
  ptr->load(input_filepath, experiment);
  TEST_EQUAL(experiment.getMetaValue("mzml_id"), "20171013_C61_ISO_P1_GA1")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getName(), "HM_metode_ZorBax_0,02%_Acetic_acid_ver6")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getSoftware().getName(), "New ProcMethod")
  TEST_EQUAL(experiment.getMetaValue("injection_date"), "10/13/2017")
  TEST_EQUAL(experiment.getMetaValue("injection_time"), "6:28:26 PM")
  TEST_EQUAL(experiment.getMetaValue("detector"), "UV")
  TEST_EQUAL(experiment.getMetaValue("signal_quantity"), "Absorbance")
  TEST_EQUAL(experiment.getMetaValue("signal_unit"), "mAU")
  TEST_EQUAL(experiment.getMetaValue("signal_info"), "WVL:280 nm")
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

  MzMLFile mzml;
  mzml.store(output_filepath, experiment);
  MSExperiment read_exp;
  mzml.load(output_filepath, read_exp);
  TEST_EQUAL(read_exp.getChromatograms().size(), 1);
  const MSChromatogram& c1 = experiment.getChromatograms()[0];
  const MSChromatogram& c2 = read_exp.getChromatograms()[0];
  TEST_EQUAL(c1.size(), c2.size())
  TEST_REAL_SIMILAR(c1[0].getRT(), c2[0].getRT())
  TEST_REAL_SIMILAR(c1[0].getIntensity(), c2[0].getIntensity())
  TEST_REAL_SIMILAR(c1[660].getRT(), c2[660].getRT())
  TEST_REAL_SIMILAR(c1[660].getIntensity(), c2[660].getIntensity())
  TEST_REAL_SIMILAR(c1[1320].getRT(), c2[1320].getRT())
  TEST_REAL_SIMILAR(c1[1320].getIntensity(), c2[1320].getIntensity())
  TEST_REAL_SIMILAR(c1[1980].getRT(), c2[1980].getRT())
  TEST_REAL_SIMILAR(c1[1980].getIntensity(), c2[1980].getIntensity())
  TEST_REAL_SIMILAR(c1[2640].getRT(), c2[2640].getRT())
  TEST_REAL_SIMILAR(c1[2640].getIntensity(), c2[2640].getIntensity())
  TEST_REAL_SIMILAR(c1[3300].getRT(), c2[3300].getRT())
  TEST_REAL_SIMILAR(c1[3300].getIntensity(), c2[3300].getIntensity())
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
