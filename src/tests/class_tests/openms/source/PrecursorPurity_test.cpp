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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>


START_TEST(PrecursorPurity, "$Id$")

using namespace OpenMS;
using namespace std;

PeakMap spectra;
MzMLFile f;
PeakFileOptions options;
options.clearMSLevels();
options.addMSLevel(1);
options.addMSLevel(2);
f.getOptions() = options;

// the file is a copy of OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML")
// which contains two MS1 spectra and 5 MS2 spectra between them
f.load(OPENMS_GET_TEST_DATA_PATH("PrecursorPurity_input.mzML"), spectra);

START_SECTION(static PurityScores computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm))

  TEST_EQUAL(spectra.size(), 7)

  Precursor pre = spectra[2].getPrecursors()[0];
  PrecursorPurity::PurityScores score = PrecursorPurity::computePrecursorPurity(spectra[6], pre, 0.2, false);
  TEST_REAL_SIMILAR(score.total_intensity, 11557777.1875)
  TEST_REAL_SIMILAR(score.target_intensity, 10320935.25)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.89298)
  TEST_EQUAL(score.target_peak_count, 2)
  TEST_EQUAL(score.residual_peak_count, 2)

  // testing a narrower tolerance for deisotoping
  score = PrecursorPurity::computePrecursorPurity(spectra[6], pre, 10, true);
  TEST_REAL_SIMILAR(score.total_intensity, 11557777.1875)
  TEST_REAL_SIMILAR(score.target_intensity, 8923915)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.77211)
  TEST_EQUAL(score.target_peak_count, 1)
  TEST_EQUAL(score.residual_peak_count, 3)

  pre = spectra[3].getPrecursors()[0];
  score = PrecursorPurity::computePrecursorPurity(spectra[0], pre, 0.2, false);

  TEST_REAL_SIMILAR(score.total_intensity, 9098343.89062)
  TEST_REAL_SIMILAR(score.target_intensity, 8266450.1875)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.90856)
  TEST_EQUAL(score.target_peak_count, 3)
  TEST_EQUAL(score.residual_peak_count, 2)

END_SECTION

START_SECTION(static std::vector<PurityScores> computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm))

  vector<PrecursorPurity::PurityScores> purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.2, false);

  TEST_REAL_SIMILAR(purityscores[0].total_intensity, 9849578.5)
  TEST_REAL_SIMILAR(purityscores[0].target_intensity, 9849578.5)
  TEST_REAL_SIMILAR(purityscores[0].signal_proportion, 1)
  TEST_EQUAL(purityscores[0].target_peak_count, 2)
  TEST_EQUAL(purityscores[0].residual_peak_count, 0)

  TEST_REAL_SIMILAR(purityscores[1].total_intensity, 22845744.8125)
  TEST_REAL_SIMILAR(purityscores[1].target_intensity, 19139751)
  TEST_REAL_SIMILAR(purityscores[1].signal_proportion, 0.83778)
  TEST_EQUAL(purityscores[1].target_peak_count, 4)
  TEST_EQUAL(purityscores[1].residual_peak_count, 4)

  TEST_REAL_SIMILAR(purityscores[2].total_intensity, 19751783.375)
  TEST_REAL_SIMILAR(purityscores[2].target_intensity, 18752920.5)
  TEST_REAL_SIMILAR(purityscores[2].signal_proportion, 0.94942)
  TEST_EQUAL(purityscores[2].target_peak_count, 6)
  TEST_EQUAL(purityscores[2].residual_peak_count, 3)

  TEST_REAL_SIMILAR(purityscores[3].total_intensity, 23979143.35156)
  TEST_REAL_SIMILAR(purityscores[3].target_intensity, 18107037.9375)
  TEST_REAL_SIMILAR(purityscores[3].signal_proportion, 0.75511)
  TEST_EQUAL(purityscores[3].target_peak_count, 4)
  TEST_EQUAL(purityscores[3].residual_peak_count, 7)

  TEST_REAL_SIMILAR(purityscores[4].total_intensity, 11964238)
  TEST_REAL_SIMILAR(purityscores[4].target_intensity, 11964238)
  TEST_REAL_SIMILAR(purityscores[4].signal_proportion, 1)
  TEST_EQUAL(purityscores[4].target_peak_count, 2)
  TEST_EQUAL(purityscores[4].residual_peak_count, 0)


END_SECTION




END_TEST
