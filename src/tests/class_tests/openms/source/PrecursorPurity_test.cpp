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
  // the MS1 spectra are soectra[0] and spectra[6]

  Precursor pre = spectra[2].getPrecursors()[0];
  PrecursorPurity::PurityScores score = PrecursorPurity::computePrecursorPurity(spectra[6], pre, 10, true);
  TEST_REAL_SIMILAR(score.total_intensity, 11557777.1875)
  TEST_REAL_SIMILAR(score.target_intensity, 8923915)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.77211)
  TEST_EQUAL(score.target_peak_count, 1)
  TEST_EQUAL(score.residual_peak_count, 3)

  pre = spectra[3].getPrecursors()[0];
  score = PrecursorPurity::computePrecursorPurity(spectra[0], pre, 0.2, false);
  TEST_REAL_SIMILAR(score.total_intensity, 9098343.89062)
  TEST_REAL_SIMILAR(score.target_intensity, 7057944)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.77573)
  TEST_EQUAL(score.target_peak_count, 1)
  TEST_EQUAL(score.residual_peak_count, 4)

END_SECTION

START_SECTION(static computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm))

  map<String, PrecursorPurity::PurityScores> purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);

  TEST_EQUAL(purityscores.size(), 5)

  // using the ID of an MS1 spectrum, a new ID for the map, adds a new Score to the map, initialized to 0
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].signal_proportion, 0)
  TEST_EQUAL(purityscores[spectra[0].getNativeID()].target_peak_count, 0)
  TEST_EQUAL(purityscores[spectra[0].getNativeID()].residual_peak_count, 0)

  // 5 MS2 spectra between two MS1 spectra
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].total_intensity, 5517171)
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].target_intensity, 5517171)
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].signal_proportion, 1)
  TEST_EQUAL(purityscores[spectra[1].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[1].getNativeID()].residual_peak_count, 0)

  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].total_intensity, 11287967.625)
  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].target_intensity, 7390478.5)
  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].signal_proportion, 0.65472)
  TEST_EQUAL(purityscores[spectra[2].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[2].getNativeID()].residual_peak_count, 3)

  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].total_intensity, 9098343.89062)
  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].target_intensity, 7057944)
  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].signal_proportion, 0.77573)
  TEST_EQUAL(purityscores[spectra[3].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[3].getNativeID()].residual_peak_count, 4)

  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].total_intensity, 9762418.03906)
  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].target_intensity, 7029896.5)
  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].signal_proportion, 0.72009)
  TEST_EQUAL(purityscores[spectra[4].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[4].getNativeID()].residual_peak_count, 5)

  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].total_intensity, 5465177)
  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].target_intensity, 5465177)
  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].signal_proportion, 1)
  TEST_EQUAL(purityscores[spectra[5].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[5].getNativeID()].residual_peak_count, 0)

  // using the ID of an MS1 spectrum, a new ID for the map, adds a new Score to the map, initialized to 0
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].signal_proportion, 0)
  TEST_EQUAL(purityscores[spectra[6].getNativeID()].target_peak_count, 0)
  TEST_EQUAL(purityscores[spectra[6].getNativeID()].residual_peak_count, 0)

  TEST_REAL_SIMILAR(purityscores["randomString"].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores["randomString"].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores["randomString"].signal_proportion, 0)
  TEST_EQUAL(purityscores["randomString"].target_peak_count, 0)
  TEST_EQUAL(purityscores["randomString"].residual_peak_count, 0)

  PeakMap spectra_copy = spectra;

  // If the PeakMap does not start with an MS1 spectrum, PrecursorPurity is not applicable and returns an empty vector
  spectra[0].setMSLevel(2);
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  spectra[5].setMSLevel(1);
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  spectra = spectra_copy;

  // duplicate IDs, skip the computation
  spectra = spectra_copy;
  TEST_EQUAL(spectra[2].getNativeID(), "controllerType=0 controllerNumber=1 scan=4")
  spectra[2].setNativeID(spectra[4].getNativeID());
  TEST_EQUAL(spectra[2].getNativeID(), "controllerType=0 controllerNumber=1 scan=8")
  TEST_EQUAL(spectra[4].getNativeID(), "controllerType=0 controllerNumber=1 scan=8")
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  // empty ID, skip the computation
  spectra = spectra_copy;
  spectra[3].setNativeID("");
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

END_SECTION




END_TEST
