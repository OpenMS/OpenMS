// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLIonLikelihood.h"
///////////////////////////

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

namespace
{
  using Model = NuXLIonLikelihood::Model;

  // b and y ions (charge 1, b1 included) as NuXLIonLikelihood generates them
  PeakSpectrum theoreticalIons(const AASequence& peptide)
  {
    TheoreticalSpectrumGenerator tsg;
    Param p = tsg.getParameters();
    p.setValue("add_first_prefix_ion", "true");
    p.setValue("add_metainfo", "true");
    tsg.setParameters(p);
    PeakSpectrum theo;
    tsg.getSpectrum(theo, peptide, 1, 1);
    return theo;
  }

  // a spectrum with one peak per ion of @p peptide; ions that contain residue @p site (if >= 0) are shifted by @p shift
  PeakMap spectrumOf(const AASequence& peptide, int site, double shift)
  {
    const PeakSpectrum theo = theoreticalIons(peptide);
    MSSpectrum spec;
    const auto& names = theo.getStringDataArrays()[0];
    for (Size i = 0; i < theo.size(); ++i)
    {
      bool prefix = false;
      Size ordinal = 0;
      Model::parseIonName(names[i], prefix, ordinal);
      const bool carries = site >= 0 && (prefix ? site < static_cast<int>(ordinal) : site >= static_cast<int>(peptide.size() - ordinal));
      spec.push_back(Peak1D(theo[i].getMZ() + (carries ? shift : 0.0), static_cast<float>(100.0 + i)));
    }
    spec.sortByPosition();
    PeakMap exp;
    exp.addSpectrum(spec);
    return exp;
  }

  Size count(const NuXLIonLikelihood::Outcomes& o, NuXLIonLikelihood::Group g, bool matched_only, Size absent)
  {
    Size n = 0;
    for (Size k = 0; k < o.ions.size(); ++k)
    {
      if (o.groups[k] == g && (!matched_only || o.ions[k].outcome != absent)) ++n;
    }
    return n;
  }
}

START_TEST(NuXLIonLikelihood, "$Id$")

const AASequence peptide = AASequence::fromString("PEPTIDEK");
const Size absent = Model().absentOutcome();

START_SECTION(static AASequence reversedSequence(const AASequence& peptide))
{
  TEST_EQUAL(NuXLIonLikelihood::reversedSequence(peptide).toString(), "EDITPEPK")
  TEST_EQUAL(NuXLIonLikelihood::reversedSequence(AASequence::fromString("PEPM(Oxidation)TIDEK")).toString(), "EDITM(Oxidation)PEPK")
  TEST_EQUAL(NuXLIonLikelihood::reversedSequence(AASequence::fromString("PK")).toString(), "PK")
}
END_SECTION

START_SECTION(Evidence match(Size scan, const AASequence& peptide, int precursor_charge, const std::vector<double>& shifts, int site) const)
{
  // linear peptide: every b and y ion is matched, the reversed sequence matches only its unchanged ions
  {
    NuXLIonLikelihood il(10.0, true);
    il.setSpectra(spectrumOf(peptide, -1, 0.0));
    const NuXLIonLikelihood::Evidence e = il.match(0, peptide, 2, {}, -1);
    TEST_EQUAL(e.signal.ions.size(), 14)
    TEST_EQUAL(count(e.signal, NuXLIonLikelihood::LINEAR, true, absent), 14)
    TEST_EQUAL(e.noise.ions.size(), 14)
    TEST_EQUAL(count(e.noise, NuXLIonLikelihood::LINEAR, true, absent) < 14, true)
  }

  // cross-link at residue 2 (P) with a fragment adduct of 100 Da: b3-b7 and y6-y7 carry it
  {
    NuXLIonLikelihood il(10.0, true);
    il.setSpectra(spectrumOf(peptide, 2, 100.0));
    const NuXLIonLikelihood::Evidence e = il.match(0, peptide, 2, {50.0, 100.0}, 2);
    TEST_EQUAL(e.signal.ions.size(), 14)
    TEST_EQUAL(count(e.signal, NuXLIonLikelihood::XL_SHIFTED, false, absent), 7)
    TEST_EQUAL(count(e.signal, NuXLIonLikelihood::XL_SHIFTED, true, absent), 7) // matched with the 100 Da adduct
    TEST_EQUAL(count(e.signal, NuXLIonLikelihood::XL_UNSHIFTED, true, absent), 7)

    // without the 100 Da adduct, the shifted ions are absent
    const NuXLIonLikelihood::Evidence other = il.match(0, peptide, 2, {50.0}, 2);
    TEST_EQUAL(count(other.signal, NuXLIonLikelihood::XL_SHIFTED, true, absent), 0)
    TEST_EQUAL(count(other.signal, NuXLIonLikelihood::XL_UNSHIFTED, true, absent), 7)

    // not localized: every ion may carry any adduct or none
    const NuXLIonLikelihood::Evidence unlocalized = il.match(0, peptide, 2, {50.0, 100.0}, -1);
    TEST_EQUAL(count(unlocalized.signal, NuXLIonLikelihood::XL_UNLOCALIZED, true, absent), 14)
  }

  // no spectrum: no ions
  {
    NuXLIonLikelihood il(10.0, true);
    il.setSpectra(PeakMap());
    TEST_EQUAL(il.match(0, peptide, 2, {}, -1).signal.ions.size(), 0)
  }
}
END_SECTION

START_SECTION(Features score(Size fold, const Outcomes& outcomes) const)
{
  NuXLIonLikelihood il(10.0, true);
  il.setSpectra(spectrumOf(peptide, -1, 0.0));
  const NuXLIonLikelihood::Evidence e = il.match(0, peptide, 2, {}, -1);

  // untrained (no observations): flat models, the ions carry no evidence
  il.finalize();
  TEST_REAL_SIMILAR(il.score(0, e.signal).llr, 0.0)

  // trained on fold 0 only: fold 0 favours the fully matched peptide over its reversed sequence, fold 1 is flat
  for (Size i = 0; i < 100; ++i) { il.addTraining(0, e); }
  il.finalize();
  const NuXLIonLikelihood::Features signal = il.score(0, e.signal);
  const NuXLIonLikelihood::Features noise = il.score(0, e.noise);
  TEST_EQUAL(signal.llr > 0.0, true)
  TEST_EQUAL(signal.llr > noise.llr, true)
  TEST_REAL_SIMILAR(signal.explained_presence, 1.0)
  TEST_REAL_SIMILAR(signal.top_predicted_observed, 1.0)
  TEST_EQUAL(noise.explained_presence < 1.0, true)
  TEST_REAL_SIMILAR(signal.llr_shifted, 0.0)
  TEST_REAL_SIMILAR(il.score(1, e.signal).llr, 0.0)

  // resetFold clears the training
  il.resetFold(0);
  il.finalize();
  TEST_REAL_SIMILAR(il.score(0, e.signal).llr, 0.0)
}
END_SECTION

END_TEST
