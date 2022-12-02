// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DeconvolvedSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DeconvolvedSpectrum* ptr = 0;
DeconvolvedSpectrum* null_ptr = 0;
START_SECTION(DeconvolvedSpectrum())
{
  ptr = new DeconvolvedSpectrum();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~DeconvolvedSpectrum())
{
  delete ptr;
}
END_SECTION


// load test data
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FLASHDeconv_sample_input1.mzML"), input);

/// detailed constructor
MSSpectrum test_spec = input[0];
START_SECTION((DeconvolvedSpectrum(const MSSpectrum &spectrum, const int scan_number)))
{
  DeconvolvedSpectrum tmp_spec = DeconvolvedSpectrum(1);
  tmp_spec.setOriginalSpectrum(test_spec);
  TEST_EQUAL(tmp_spec.getScanNumber(), 1);
  TEST_EQUAL(tmp_spec.getOriginalSpectrum().size(), test_spec.size());
}
END_SECTION


////////
DeconvolvedSpectrum test_deconv_spec = DeconvolvedSpectrum(1);
START_SECTION((int getScanNumber() const))
{
  test_deconv_spec.setOriginalSpectrum(test_spec);
  int tmp_num = test_deconv_spec.getScanNumber();
  TEST_EQUAL(tmp_num, 1);
}
END_SECTION

START_SECTION((const MSSpectrum& getOriginalSpectrum() const))
{
  MSSpectrum tmp_s = test_deconv_spec.getOriginalSpectrum();
  TEST_EQUAL(tmp_s.size(), test_spec.size());
}
END_SECTION


FLASHDeconvAlgorithm fd_algo = FLASHDeconvAlgorithm();
Param fd_param;
fd_param.setValue("min_charge", 5);
fd_param.setValue("max_charge", 20);
fd_algo.setParameters(fd_param);
fd_algo.calculateAveragine(false);
std::vector<DeconvolvedSpectrum> survey_specs;
const std::map<int, std::vector<std::vector<float>>> null_map;

    fd_algo.performSpectrumDeconvolution(input[1], survey_specs, 2, null_map);
DeconvolvedSpectrum prec_deconv_spec_1 = fd_algo.getDeconvolvedSpectrum();

    fd_algo.performSpectrumDeconvolution(input[3], survey_specs, 4, null_map);
DeconvolvedSpectrum prec_deconv_spec_2 = fd_algo.getDeconvolvedSpectrum();

survey_specs.push_back(prec_deconv_spec_2);
    fd_algo.performSpectrumDeconvolution(input[5], survey_specs, 6, null_map);
DeconvolvedSpectrum ms2_deconv_spec = fd_algo.getDeconvolvedSpectrum();

START_SECTION((double getCurrentMaxMass(const double max_mass) const))
{
  double ms1_max_mass = test_deconv_spec.getCurrentMaxMass(1000.);
  double ms2_max_mass = ms2_deconv_spec.getCurrentMaxMass(13673.239337872);
  TOLERANCE_ABSOLUTE(1);
  TEST_REAL_SIMILAR(ms1_max_mass, 1000.);
  TEST_REAL_SIMILAR(ms2_max_mass, 13674.);
}
END_SECTION

START_SECTION((double getCurrentMinMass(const double min_mass) const))
{
  double ms1_min_mass = test_deconv_spec.getCurrentMinMass(1000.);
  double ms2_min_mass = ms2_deconv_spec.getCurrentMinMass(1000.);
  TEST_REAL_SIMILAR(ms1_min_mass, 1000.);
  TEST_REAL_SIMILAR(ms2_min_mass, 50.);
}
END_SECTION

START_SECTION((MSSpectrum toSpectrum(const int mass_charge)))
{
  MSSpectrum peakgroup_spec = prec_deconv_spec_1.toSpectrum(9);
  TEST_EQUAL(peakgroup_spec.size(), 4);
  TEST_REAL_SIMILAR(peakgroup_spec.getRT(), 251.72280736002);
}
END_SECTION

START_SECTION((PeakGroup getPrecursorPeakGroup() const))
{
  PeakGroup tmp_precursor_pgs = ms2_deconv_spec.getPrecursorPeakGroup();

  TEST_EQUAL(tmp_precursor_pgs.size(), 65);
  TOLERANCE_ABSOLUTE(5);
  TEST_REAL_SIMILAR(tmp_precursor_pgs.getMonoMass(), 13674.2798657377);
  TEST_REAL_SIMILAR(tmp_precursor_pgs.getIntensity(), 232085.069824219);
  TEST_EQUAL(tmp_precursor_pgs.getScanNumber(), 4);
}
END_SECTION

START_SECTION((const Precursor getPrecursor() const))
{
  Precursor tmp_precursor = ms2_deconv_spec.getPrecursor();
  TEST_EQUAL(tmp_precursor.getCharge(), 9);
  TOLERANCE_ABSOLUTE(10);
  TEST_REAL_SIMILAR(tmp_precursor.getUnchargedMass(), 13682.3053614085);
  TEST_REAL_SIMILAR(tmp_precursor.getIntensity(), 12293.4);
}
END_SECTION

START_SECTION((int getPrecursorCharge() const))
{
  int prec_cs = ms2_deconv_spec.getPrecursorCharge();
  TEST_EQUAL(prec_cs, 9);
}
END_SECTION

START_SECTION((int getPrecursorScanNumber() const))
{
  int p_scan_num = ms2_deconv_spec.getPrecursorScanNumber();
  TEST_EQUAL(p_scan_num, 4);
}
END_SECTION

START_SECTION((int getCurrentMaxAbsCharge(const int max_abs_charge) const))
{
  int tmp_cs_ms1 = test_deconv_spec.getCurrentMaxAbsCharge(5);
  int tmp_cs_ms2 = ms2_deconv_spec.getCurrentMaxAbsCharge(5);

  TEST_EQUAL(tmp_cs_ms1, 5);
  TEST_EQUAL(tmp_cs_ms2, 9);
}
END_SECTION

START_SECTION(String& getActivationMethod() const)
{
  Precursor::ActivationMethod act_method = ms2_deconv_spec.getActivationMethod();
  TEST_EQUAL(Precursor::ActivationMethod::ETD, act_method); // TODO: why ETD?
}
END_SECTION
////////

/// < public methods without tests > : TODOs
/// - default constructors and operators are not used (copy, move, assignment)
/// - setters (setPrecursor, etc.)
/// - updatePeakGroupQvalues
/// - nested stuff

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST