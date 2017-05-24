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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

using namespace OpenMS;

START_TEST(OPXLSpectrumProcessingAlgorithms, "$Id$")
TheoreticalSpectrumGeneratorXLMS specGen;
Param param = specGen.getParameters();
param.setValue("add_isotopes", "false");
param.setValue("add_metainfo", "false");
param.setValue("add_first_prefix_ion", "false");
param.setValue("y_intensity", 10.0, "Intensity of the y-ions");
specGen.setParameters(param);

PeakSpectrum theo_spec_1, theo_spec_2, exp_spec_1, exp_spec_2;
specGen.getCommonIonSpectrum(exp_spec_1, AASequence::fromString("PEPTIDE"), 2, true, 3);
specGen.getCommonIonSpectrum(exp_spec_2, AASequence::fromString("PEPTEDI"), 3, true, 3);

param.setValue("add_metainfo", "true");
specGen.setParameters(param);

specGen.getCommonIonSpectrum(theo_spec_1, AASequence::fromString("PEPTIDE"), 3, true, 3);
specGen.getCommonIonSpectrum(theo_spec_2, AASequence::fromString("PEPTEDI"), 4, true, 3);

START_SECTION(static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum))

  PeakSpectrum merged_spec = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theo_spec_1, theo_spec_2);

  TEST_EQUAL(merged_spec.size(), 36)
  TEST_EQUAL(merged_spec.getIntegerDataArrays().size(), 1)
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0].size(), 36)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0].size(), 36)
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0][10], 3)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0][10], "[alpha|ci$y2]")
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0][20], 2)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0][20], "[alpha|ci$y2]")
  TEST_REAL_SIMILAR(merged_spec[10].getMZ(), 83.04780)
  TEST_REAL_SIMILAR(merged_spec[20].getMZ(), 132.04732)

  for (Size i = 0; i < merged_spec.size()-1; ++i)
  {
    TEST_EQUAL(merged_spec[i].getMZ() <= merged_spec[i+1].getMZ(), true)
  }

END_SECTION

START_SECTION(static void getSpectrumAlignment(std::vector <std::pair <Size, Size> >& alignment, const PeakSpectrum & s1, const PeakSpectrum & s2, double tolerance, bool relative_tolerance, double intensity_cutoff = 0.0))
  std::vector <std::pair <Size, Size> > alignment1;
  std::vector <std::pair <Size, Size> > alignment2;

  theo_spec_1.sortByPosition();

  OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment1, theo_spec_1, exp_spec_1, 50, true);
  OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(alignment2, theo_spec_2, exp_spec_2, 50, true);

  TEST_EQUAL(alignment1.size(), 15)
  TEST_EQUAL(alignment2.size(), 15)
  for (Size i = 0; i < alignment2.size(); i++)
  {
    TEST_REAL_SIMILAR(theo_spec_2[alignment2[i].first].getMZ(), exp_spec_2[alignment2[i].second].getMZ())
  }

END_SECTION

START_SECTION(static PeakSpectrum deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = false))

  // this is basically testing, that nothing weird happens to an already deisotoped spectrum
  // TODO expand tests, use spectra with isotopic patterns
  // maybe read in one of the OPXL test files?
  PeakSpectrum deisotoped_spec = OPXLSpectrumProcessingAlgorithms::deisotopeAndSingleChargeMSSpectrum(theo_spec_2, 1, 3, 50, true);

  TEST_EQUAL(deisotoped_spec.size(), theo_spec_2.size())

END_SECTION

END_TEST
