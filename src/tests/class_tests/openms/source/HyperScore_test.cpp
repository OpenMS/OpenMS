// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

using namespace OpenMS;
using namespace std;

START_TEST(HyperScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

HyperScore* ptr = nullptr;
HyperScore* null_ptr = nullptr;

TheoreticalSpectrumGenerator tsg;
Param param = tsg.getParameters();
param.setValue("add_metainfo", "true");
tsg.setParameters(param);

START_SECTION(HyperScore())
{
  ptr = new HyperScore();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~HyperScore())
{
  delete ptr;
}
END_SECTION

START_SECTION(([EXTRA] mass - accuracy score preserves exact HyperScore and discounts dispersed matches))
{
  MSSpectrum theory;
  theory.getStringDataArrays().emplace_back();
  for (Size i = 0; i < 4; ++i)
  {
    theory.emplace_back(500.0 + i * 100.0, 1.0);
    theory.getStringDataArrays()[0].push_back(i < 2 ? "b2+" : "y2+");
  }
  HyperScore::PSMDetail detail, original;
  const double baseline = HyperScore::computeWithDetail(20.0, true, theory, theory, original);
  TEST_REAL_SIMILAR(HyperScore::computeMassAccuracy(20.0, true, theory, theory, 7.0, detail), baseline)
  TEST_EQUAL(detail.matched_prefix_ions, 2)
  TEST_EQUAL(detail.matched_suffix_ions, 2)
  TEST_REAL_SIMILAR(detail.mean_error, 0.0)
  for (double sign : {-1.0, 1.0})
  {
    MSSpectrum shifted = theory;
    for (auto& peak : shifted)
      peak.setMZ(peak.getMZ() * (1.0 + sign * 7e-6));
    const double weight = std::exp(-0.5);
    const double expected = std::log1p(4.0 * weight) + 2.0 * std::lgamma(2.0 * weight + 1.0);
    const double score = HyperScore::computeMassAccuracy(20.0, true, shifted, theory, 7.0, detail);
    TEST_REAL_SIMILAR(score, expected)
    TEST_TRUE(score < baseline)
    TEST_EQUAL(detail.matched_prefix_ions, 2)
    TEST_EQUAL(detail.matched_suffix_ions, 2)
    TEST_REAL_SIMILAR(detail.mean_error, 7.0)
    TEST_REAL_SIMILAR(HyperScore::computeMassAccuracy(0.02, false, shifted, theory, 7.0, detail), expected)
    TEST_REAL_SIMILAR(detail.mean_error, 0.00455)
  }
  MSSpectrum singleton;
  singleton.emplace_back(500.0035, 1.0);
  const double single = HyperScore::computeMassAccuracy(20.0, true, singleton, theory, 7.0, detail);
  TEST_REAL_SIMILAR(single, std::log1p(std::exp(-0.5)))
  TEST_TRUE(single > 0.0)
  TEST_REAL_SIMILAR(HyperScore::computeMassAccuracy(20.0, true, MSSpectrum {}, theory, 7.0, detail), 0.0)
  TEST_EQUAL(detail.matched_prefix_ions, 0)
  TEST_REAL_SIMILAR(detail.mean_error, 0.0)
  MSSpectrum outside;
  outside.emplace_back(100.0, 1.0);
  TEST_REAL_SIMILAR(HyperScore::computeMassAccuracy(20.0, true, outside, theory, 7.0, detail), 0.0)
  for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
  {
    TEST_EXCEPTION(Exception::InvalidParameter, HyperScore::computeMassAccuracy(20.0, true, theory, theory, invalid, detail))
    TEST_EXCEPTION(Exception::InvalidParameter, HyperScore::computeMassAccuracy(invalid, true, theory, theory, 7.0, detail))
  }
  theory.getStringDataArrays().clear();
  TEST_EXCEPTION(Exception::InvalidValue, HyperScore::computeMassAccuracy(20.0, true, theory, theory, 7.0, detail))
}
END_SECTION

START_SECTION((static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum &exp_spectrum, const RichPeakSpectrum &theo_spectrum)))
{
  PeakSpectrum exp_spectrum;
  PeakSpectrum theo_spectrum;

  AASequence peptide = AASequence::fromString("PEPTIDE");
  
  // empty spectrum
  tsg.getSpectrum(theo_spectrum, peptide, 1, 1);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 0.0);

  // full match, 11 identical masses, identical intensities (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 1);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 13.8516496);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 13.8516496);

  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // no match
  tsg.getSpectrum(exp_spectrum, peptide, 1, 3);
  tsg.getSpectrum(theo_spectrum, AASequence::fromString("YYYYYY"), 1, 3);
  TEST_REAL_SIMILAR(HyperScore::compute(1e-5, false, exp_spectrum, theo_spectrum), 0.0);
  
  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // full match, 33 identical masses, identical intensities (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 3);
  tsg.getSpectrum(theo_spectrum, peptide, 1, 3);
  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 67.8210771);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 67.8210771);

  // full match if ppm tolerance and partial match for Da tolerance
  for (Size i = 0; i < theo_spectrum.size(); ++i)
  {
    double mz = pow( theo_spectrum[i].getMZ(), 2);
    exp_spectrum[i].setMZ(mz);
    theo_spectrum[i].setMZ(mz + 9 * 1e-6 * mz); // +9 ppm error
  }

  TEST_REAL_SIMILAR(HyperScore::compute(0.1, false, exp_spectrum, theo_spectrum), 3.401197);
  TEST_REAL_SIMILAR(HyperScore::compute(10, true, exp_spectrum, theo_spectrum), 67.8210771);
}
END_SECTION

START_SECTION([EXTRA] experimental peak listed twice)
{
  // A second peak with the same m/z (as in merged or pseudo-MS/MS spectra) must not hide the peaks
  // above it: the score and the ion counts equal those of the spectrum without the copy.
  PeakSpectrum theo_spectrum, exp_spectrum;
  tsg.getSpectrum(theo_spectrum, AASequence::fromString("LISWYDNEFGYSNR"), 1, 1); // b2-b13, y1-y13
  for (const Peak1D& p : theo_spectrum)
  {
    exp_spectrum.push_back(p);
  }
  PeakSpectrum exp_dup = exp_spectrum;
  exp_dup.push_back(exp_spectrum[0]); // y1 twice, below all other ions
  exp_dup.sortByPosition();

  for (bool ppm : {true, false})
  {
    const double tol = ppm ? 10.0 : 0.01;
    HyperScore::PSMDetail d, d_dup;
    TEST_REAL_SIMILAR(HyperScore::computeWithDetail(tol, ppm, exp_spectrum, theo_spectrum, d), 45.7974749)
    TEST_REAL_SIMILAR(HyperScore::computeWithDetail(tol, ppm, exp_dup, theo_spectrum, d_dup), 45.7974749)
    TEST_EQUAL(d_dup.matched_prefix_ions, 12)
    TEST_EQUAL(d_dup.matched_suffix_ions, 13)
    TEST_REAL_SIMILAR(HyperScore::compute(tol, ppm, exp_dup, theo_spectrum), 45.7974749)
  }
}
END_SECTION

START_SECTION([EXTRA] computeWithDetail: mean fragment error in Da)
{
  // Every matched ion is 0.004 Da off. The mean error must be 0.004, not 0: an unqualified abs()
  // resolved to the integer overload where no header exported the double one (aarch64), which
  // truncated each error and made ProSE's fragment calibration keep the configured tolerance.
  PeakSpectrum theo_spectrum, exp_spectrum;
  tsg.getSpectrum(theo_spectrum, AASequence::fromString("PEPTIDER"), 1, 1);
  for (const Peak1D& p : theo_spectrum)
  {
    exp_spectrum.emplace_back(p.getMZ() + 0.004, p.getIntensity());
  }
  HyperScore::PSMDetail d;
  HyperScore::computeWithDetail(0.01, false, exp_spectrum, theo_spectrum, d);
  TEST_EQUAL(d.matched_prefix_ions + d.matched_suffix_ions, theo_spectrum.size())
  TEST_REAL_SIMILAR(d.mean_error, 0.004)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

