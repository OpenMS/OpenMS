// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <algorithm>
#include <cmath>

using namespace OpenMS;
using namespace std;

namespace
{
  // add the averagine envelope of a precursor (mono m/z, charge) with the given total intensity
  void addEnvelope(MSSpectrum& s, double mono_mz, int charge, double total)
  {
    const double mass = (mono_mz - Constants::PROTON_MASS_U) * charge;
    const IsotopeDistribution dist = CoarseIsotopePatternGenerator(5).estimateFromPeptideWeight(mass);
    for (Size k = 0; k < dist.size(); ++k)
    {
      s.emplace_back(mono_mz + k * Constants::C13C12_MASSDIFF_U / charge, total * dist[k].getIntensity());
    }
  }
}

START_TEST(CoIsolationDetector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static double isotopeKL(const MSSpectrum& ms1, double mono_mz, int charge, double tol_da, Size n_isotopes = 4)))
{
  MSSpectrum s;
  addEnvelope(s, 600.3, 2, 1e6);
  s.sortByPosition();
  // the right charge fits averagine, the wrong one does not, and no peaks gives 10
  TEST_EQUAL(CoIsolationDetector::isotopeKL(s, 600.3, 2, 0.01) < 0.01, true)
  TEST_EQUAL(CoIsolationDetector::isotopeKL(s, 600.3, 3, 0.01) > 0.3, true)
  TEST_REAL_SIMILAR(CoIsolationDetector::isotopeKL(s, 700.0, 2, 0.01), 10.0)
  TEST_REAL_SIMILAR(CoIsolationDetector::isotopeKL(MSSpectrum(), 600.3, 2, 0.01), 10.0)
}
END_SECTION

START_SECTION((static std::vector<CoIsolated> detect(const MSSpectrum& ms1, double window_lower, double window_upper, double selected_mz, int min_charge, int max_charge, double tol_da, double max_kl, Size max_n)))
{
  MSSpectrum s;
  addEnvelope(s, 600.30, 2, 1e6);  // selected precursor
  addEnvelope(s, 600.95, 3, 4e5);  // co-isolated
  addEnvelope(s, 605.00, 2, 9e6);  // outside the window
  s.sortByPosition();

  auto found = CoIsolationDetector::detect(s, 599.8, 601.8, 600.30, 1, 4, 0.01, 0.3, 4);
  TEST_EQUAL(found.size(), 1)
  ABORT_IF(found.size() != 1)
  TEST_REAL_SIMILAR(found[0].mono_mz, 600.95)
  TEST_EQUAL(found[0].charge, 3)

  // max_n 0 finds nothing; an empty window finds nothing
  TEST_EQUAL(CoIsolationDetector::detect(s, 599.8, 601.8, 600.30, 1, 4, 0.01, 0.3, 0).size(), 0)
  TEST_EQUAL(CoIsolationDetector::detect(s, 620.0, 621.0, 600.30, 1, 4, 0.01, 0.3, 4).size(), 0)

  // a heavy precursor whose M+1 is its most intense peak: reported at its monoisotopic peak
  MSSpectrum heavy;
  addEnvelope(heavy, 500.00, 2, 1e6);   // selected precursor
  addEnvelope(heavy, 1250.60, 2, 5e6);  // 2.5 kDa, M+1 > M
  heavy.sortByPosition();
  found = CoIsolationDetector::detect(heavy, 1249.0, 1252.5, 500.00, 1, 4, 0.01, 0.3, 4);
  TEST_EQUAL(found.empty(), false)
  ABORT_IF(found.empty())
  TEST_REAL_SIMILAR(found[0].mono_mz, 1250.60)
  TEST_EQUAL(found[0].charge, 2)
}
END_SECTION

START_SECTION((static Size annotate(PeakMap& exp, int min_charge, int max_charge, double tolerance, bool tolerance_ppm, double max_kl, Size max_n)))
{
  PeakMap exp;
  MSSpectrum ms1;
  ms1.setMSLevel(1);
  ms1.setRT(10.0);
  ms1.setNativeID("scan=1");
  addEnvelope(ms1, 600.30, 2, 1e6);
  addEnvelope(ms1, 600.95, 3, 4e5);
  ms1.sortByPosition();
  exp.addSpectrum(ms1);

  MSSpectrum ms2;
  ms2.setMSLevel(2);
  ms2.setRT(10.5);
  Precursor p;
  p.setMZ(600.30);
  p.setCharge(2);
  p.setIsolationWindowLowerOffset(1.0);
  p.setIsolationWindowUpperOffset(1.0);
  ms2.setPrecursors({p});
  exp.addSpectrum(ms2);

  // without an isolation window: not annotated
  MSSpectrum no_window = ms2;
  Precursor q = p;
  q.setIsolationWindowLowerOffset(0.0);
  q.setIsolationWindowUpperOffset(0.0);
  no_window.setPrecursors({q});
  exp.addSpectrum(no_window);

  TEST_EQUAL(CoIsolationDetector::annotate(exp, 1, 4, 10.0, true, 0.3, 4), 1)
  TEST_EQUAL(exp[1].metaValueExists(CoIsolationDetector::META_MZ), true)
  const DoubleList mzs = exp[1].getMetaValue(CoIsolationDetector::META_MZ);
  const IntList charges = exp[1].getMetaValue(CoIsolationDetector::META_CHARGE);
  TEST_EQUAL(mzs.size(), 1)
  TEST_EQUAL(charges.size(), 1)
  ABORT_IF(mzs.size() != 1 || charges.size() != 1)
  TEST_REAL_SIMILAR(mzs[0], 600.95)
  TEST_EQUAL(charges[0], 3)
  TEST_EQUAL(exp[0].metaValueExists(CoIsolationDetector::META_MZ), false)
  TEST_EQUAL(exp[2].metaValueExists(CoIsolationDetector::META_MZ), false)

  // the tolerance in ppm of the selected precursor's m/z or in Th: isotope peaks 0.005 Th off their positions are
  // matched within 0.01 Th (17 ppm), not within 1 ppm
  PeakMap shifted;
  MSSpectrum ms1_shifted;
  ms1_shifted.setMSLevel(1);
  ms1_shifted.setRT(10.0);
  addEnvelope(ms1_shifted, 600.30, 2, 1e6);
  MSSpectrum envelope;
  addEnvelope(envelope, 600.95, 3, 4e5);
  for (Size k = 0; k < envelope.size(); ++k)
  {
    ms1_shifted.emplace_back(envelope[k].getMZ() + (k == 0 ? 0.0 : 0.005), envelope[k].getIntensity());
  }
  ms1_shifted.sortByPosition();
  shifted.addSpectrum(ms1_shifted);
  shifted.addSpectrum(ms2);
  auto found_at = [&shifted](double tolerance, bool ppm)
  {
    CoIsolationDetector::annotate(shifted, 1, 4, tolerance, ppm, 0.3, 4);
    const DoubleList mzs = shifted[1].getMetaValue(CoIsolationDetector::META_MZ);
    return std::any_of(mzs.begin(), mzs.end(), [](double mz) { return std::fabs(mz - 600.95) < 0.001; });
  };
  TEST_EQUAL(found_at(0.01, false), true)
  TEST_EQUAL(found_at(17.0, true), true)
  TEST_EQUAL(found_at(1.0, true), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
