// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerRapid.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerRapid, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerRapid* ptr = 0;
PeakPickerRapid* null_ptr = 0;
START_SECTION(PeakPickerRapid())
{
    ptr = new PeakPickerRapid();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~PeakPickerRapid()))
{
  delete ptr;
}
END_SECTION


PeakPickerRapid ppr;
MSExperiment<Peak1D> input, output;
// load Orbitrap input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"),input);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn0_out.mzML"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
	output[scan_idx].setType(SpectrumSettings::PEAKS);
}

START_SECTION((template < typename PeakType > bool computeTPG(const PeakType &p1, const PeakType &p2, const PeakType &p3, DoubleReal &mu, DoubleReal &sigma, DoubleReal &area, DoubleReal &height) const ))
{
    // test mean of gaussian if 3 non-symmetric points are given
    DoubleReal mean, s, area, height;
    Peak1D p1, p2, p3;
    p1.setMZ(100.5);
    p1.setIntensity(0.3520653);

    p2.setMZ(101.0);
    p2.setIntensity(0.3989423);

    p3.setMZ(101.6);
    p3.setIntensity(0.3332246);

    ppr.computeTPG(p1, p2, p3, mean, s, area, height);
    TEST_REAL_SIMILAR(mean, 101.0)
    TEST_REAL_SIMILAR(s, 1.0)
    TEST_REAL_SIMILAR(area, 1.0)
    TEST_REAL_SIMILAR(height, 1.0/std::sqrt(2.0*OpenMS::Constants::PI))

    // test height and area of scaled gaussian (factor = 100)
    p1.setMZ(-0.5);
    p1.setIntensity(100.0*0.3520653);

    p2.setMZ(0.0);
    p2.setIntensity(100.0*0.3989423);

    p3.setMZ(0.6);
    p3.setIntensity(100.0*0.3332246);

    ppr.computeTPG(p1, p2, p3, mean, s, area, height);
    TEST_REAL_SIMILAR(mean, 0.0)
    TEST_REAL_SIMILAR(s, 1.0)
    TEST_REAL_SIMILAR(area, 100.0)
    TEST_REAL_SIMILAR(height, 100.0/std::sqrt(2.0*OpenMS::Constants::PI))
}
END_SECTION

START_SECTION((template < typename PeakType > void pick(const MSSpectrum< PeakType > &cinput, MSSpectrum< PeakType > &output)))
{
  // should find the same peaks as spline based peak picker (PeakPickerHiRes)
  MSSpectrum<Peak1D> tmp_spec;
  ppr.pick(input[0],tmp_spec);

  TEST_EQUAL(tmp_spec.size(),output[0].size());
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx+1].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx+1].getIntensity())
  }
}
END_SECTION

START_SECTION((template < typename PeakType > void pickExperiment(MSExperiment< PeakType > &input, MSExperiment< PeakType > &output)))
{
	/*
  MSExperiment<Peak1D> tmp_exp;
  ppr.pickExperiment(input,tmp_exp);

  TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
    {
     TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
     TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
    }
  }
*/
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



