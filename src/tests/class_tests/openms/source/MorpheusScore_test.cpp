// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/MorpheusScore.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

using namespace OpenMS;
using namespace std;

START_TEST(MorpheusScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MorpheusScore* ptr = 0;
MorpheusScore* null_ptr = 0;

TheoreticalSpectrumGenerator tsg;
Param param = tsg.getParameters();
param.setValue("add_metainfo", "true");
tsg.setParameters(param);

START_SECTION(MorpheusScore())
{
  ptr = new MorpheusScore();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MorpheusScore())
{
  delete ptr;
}
END_SECTION


START_SECTION((static MorpheusScore::Result compute(
  double fragment_mass_tolerance, 
  bool fragment_mass_tolerance_unit_ppm, 
  const PeakSpectrum &exp_spectrum, 
  const PeakSpectrum &theo_spectrum)))
{
  PeakSpectrum exp_spectrum;
  PeakSpectrum theo_spectrum;

  AASequence peptide = AASequence::fromString("PEPTIDE");
  
  // empty spectrum
  tsg.getSpectrum(theo_spectrum, peptide, 1, 1);
  TEST_REAL_SIMILAR(MorpheusScore::compute(0.1, false, exp_spectrum, theo_spectrum).score, 0.0);

  // full match, 11 identical masses, identical intensities (=1). Total score is number of matches (11) + fraction of TIC (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 1);
  TEST_EQUAL(exp_spectrum.size(), 11);
  TEST_EQUAL(theo_spectrum.size(), 11);
  TEST_REAL_SIMILAR(MorpheusScore::compute(0.1, false, exp_spectrum, theo_spectrum).score, 11.0 + 1.0);
  TEST_REAL_SIMILAR(MorpheusScore::compute(10, true, exp_spectrum, theo_spectrum).score, 11.0 + 1.0);

  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // no match
  tsg.getSpectrum(exp_spectrum, peptide, 1, 1);
  tsg.getSpectrum(theo_spectrum, AASequence::fromString("EDITPEP"), 1, 1);
  TEST_REAL_SIMILAR(MorpheusScore::compute(1e-5, false, exp_spectrum, theo_spectrum).score, 0.0);
  
  exp_spectrum.clear(true);
  theo_spectrum.clear(true);

  // full match, 33 identical masses, identical intensities (=1)
  tsg.getSpectrum(exp_spectrum, peptide, 1, 3);
  tsg.getSpectrum(theo_spectrum, peptide, 1, 3);
  TEST_REAL_SIMILAR(MorpheusScore::compute(0.1, false, exp_spectrum, theo_spectrum).score, 33.0 + 1.0);
  TEST_REAL_SIMILAR(MorpheusScore::compute(10, true, exp_spectrum, theo_spectrum).score, 33.0 + 1.0);

  // full match if ppm tolerance and partial match for Da tolerance
  for (Size i = 0; i < theo_spectrum.size(); ++i)
  {
    double mz = pow( theo_spectrum[i].getMZ(), 2);
    exp_spectrum[i].setMZ(mz);
    theo_spectrum[i].setMZ(mz + 9 * 1e-6 * mz); // +9 ppm error
  }

  TEST_EQUAL(MorpheusScore::compute(0.1, false, exp_spectrum, theo_spectrum).matches, 4);
  TEST_REAL_SIMILAR(MorpheusScore::compute(0.1, false, exp_spectrum, theo_spectrum).score, 4.1212);
  TEST_EQUAL(MorpheusScore::compute(10, true, exp_spectrum, theo_spectrum).matches, 33);
  TEST_REAL_SIMILAR(MorpheusScore::compute(10, true, exp_spectrum, theo_spectrum).score, 33.0 + 1.0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

