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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerCWT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerCWT* ptr = nullptr;
PeakPickerCWT* nullPointer = nullptr;
START_SECTION((PeakPickerCWT()))
  ptr = new PeakPickerCWT();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerCWT()))
  delete ptr;
END_SECTION

// load input and output data
MzMLFile mz_ml_file;
PeakMap input, output;
mz_ml_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT_test.mzML"), input);
mz_ml_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT_test_output.mzML"), output);

// set data type (this is not stored correctly in mzData)
for (Size s=0; s<output.size(); ++s)
{
  output[s].setType(SpectrumSettings::CENTROID);
}

// set up PeakPicker
PeakPickerCWT pp;
Param param;
param.setValue("peak_width", 0.15);
param.setValue("signal_to_noise", 3.);
pp.setParameters(param);   

START_SECTION((void pick(const MSSpectrum &input, MSSpectrum &output) const))
  MSSpectrum spec;
  pp.pick(input[0],spec);
  //updating results: PeakMap sp; sp.addSpectrum(spec); MzMLFile().store(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT_test_output.mzML"), sp);

// TEST_EQUAL(spec.SpectrumSettings::operator==(output[0]), true)-> are not equal as peak picking step is written to the spectrum settings
  for (Size p=0; p<spec.size(); ++p)
  {
    TEST_REAL_SIMILAR(spec[p].getMZ(), output[0][p].getMZ())
    TEST_REAL_SIMILAR(spec[p].getIntensity(), output[0][p].getIntensity())
  }
END_SECTION


START_SECTION((void pickExperiment(const PeakMap &input, PeakMap &output)))
  PeakMap exp;
  pp.pickExperiment(input, exp);
  TEST_EQUAL(exp.ExperimentalSettings::operator==(input), true)

  ABORT_IF(exp.size() != output.size())
  for (Size s = 0; s < exp.size(); ++s)
  {
    for (Size p = 0; p<exp[s].size(); ++p)
    {
      TEST_REAL_SIMILAR(exp[s][p].getMZ(), output[s][p].getMZ())
      TEST_REAL_SIMILAR(exp[s][p].getIntensity(), output[s][p].getIntensity())
    }
  }
END_SECTION

START_SECTION(double estimatePeakWidth(const PeakMap& input))
  PeakPickerCWT pp;
  // add empty spectra.. make sure that the algorithm does not stumble
  double peak_width = pp.estimatePeakWidth(input);
  TEST_REAL_SIMILAR(peak_width, 0.155283369123936)
  // add empty spectra.. make sure that the algorithm does not stumble
  input.addSpectrum(MSSpectrum());
  input.addSpectrum(MSSpectrum());
  peak_width = pp.estimatePeakWidth(input);
  TEST_REAL_SIMILAR(peak_width, 0.155283369123936)
  // test on empty container
  input.clear(true);
  peak_width = pp.estimatePeakWidth(input);
  TEST_REAL_SIMILAR(peak_width, 0.0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



