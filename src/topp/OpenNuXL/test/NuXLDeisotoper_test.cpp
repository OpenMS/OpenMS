// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLDeisotoper.h"
///////////////////////////

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

using namespace OpenMS;
using namespace std;

START_TEST(NuXLDeisotoper, "$Id$")

START_SECTION([EXTRA] precursor mass constraint of deisotopeAndSingleCharge)
{
  // one doubly charged isotope cluster at m/z 200 (neutral monoisotopic mass 2 * (200 - proton) = 397.985).
  // Its first and last peak are also 1 u apart, so without the charge 2 interpretation it is read as charge 1.
  MSSpectrum base;
  Peak1D pk;
  pk.setIntensity(1.0);
  pk.setMZ(200.0);                                       base.push_back(pk);
  pk.setMZ(200.0 + 0.5 * Constants::C13C12_MASSDIFF_U);  base.push_back(pk);
  pk.setMZ(200.0 + 1.0 * Constants::C13C12_MASSDIFF_U);  base.push_back(pk);

  // deisotopes (10 ppm, charges 1..2, keep only deisotoped peaks, keep charges, annotate charge) and returns the charge
  // assigned to the peak at m/z 200, or 0 if that peak was removed
  auto charge_at_200 = [](MSSpectrum s)
  {
    NuXLDeisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, true, 2, 10, false, true);
    if (s.empty() || s[0].getMZ() > 200.001) return 0;
    for (const auto& array : s.getIntegerDataArrays())
    {
      if (array.getName() == "charge") return array[0];
    }
    return -1;
  };

  // no precursor: charge 2
  TEST_EQUAL(charge_at_200(base), 2)

  // precursor with unknown charge (0): no precursor mass constraint, same result (used to reject every cluster)
  MSSpectrum s = base;
  Precursor unknown_charge;
  unknown_charge.setMZ(200.0);
  unknown_charge.setCharge(0);
  s.setPrecursors({unknown_charge});
  TEST_EQUAL(charge_at_200(s), 2)

  // singly charged precursor at m/z 399.0 (neutral mass 397.993): the charge 2 cluster is slightly lighter and allowed.
  // Masses are in u, so the proton mass must be Constants::PROTON_MASS_U (with the kg constant 400 > 399 rejected it).
  s = base;
  Precursor precursor;
  precursor.setMZ(399.0);
  precursor.setCharge(1);
  s.setPrecursors({precursor});
  TEST_EQUAL(charge_at_200(s), 2)

  // a precursor lighter than the charge 2 cluster leaves only the charge 1 interpretation
  precursor.setMZ(390.0);
  s.setPrecursors({precursor});
  TEST_EQUAL(charge_at_200(s), 1)
}
END_SECTION

START_SECTION([EXTRA] first_isotope_check_max_mass of deisotopeAndSingleCharge)
{
  // a weak unrelated peak 1 u below a fragment ion (S30 spectrum 29492: m/z 445.2006 at 2.1 % and b3+U' at 446.2039, 3.8 %)
  MSSpectrum base;
  Peak1D pk;
  pk.setMZ(445.2006); pk.setIntensity(2.07); base.push_back(pk);
  pk.setMZ(446.2039); pk.setIntensity(3.78); base.push_back(pk);

  // OpenNuXL's settings (0.01 Da, charges 1..3, keep unassigned peaks, 2..10 isotopic peaks, decreasing model from the
  // second isotope on, add up intensities)
  auto deisotope = [](MSSpectrum s, double max_mass)
  {
    NuXLDeisotoper::deisotopeAndSingleCharge(s, 0.01, false, 1, 3, false, 2, 10, false, false, false, true, 2, true, false, false, -1e10, max_mass);
    return s;
  };

  // without the check, the stronger peak is taken as the first isotope of the weaker one and removed
  MSSpectrum s = deisotope(base, 0.0);
  TEST_EQUAL(s.size(), 1)
  TEST_REAL_SIMILAR(s[0].getMZ(), 445.2006)

  // with it, both peaks are kept with their own intensities
  s = deisotope(base, 1500.0);
  TEST_EQUAL(s.size(), 2)
  TEST_REAL_SIMILAR(s[1].getMZ(), 446.2039)
  TEST_REAL_SIMILAR(s[1].getIntensity(), 3.78)

  // a real isotope pattern (first isotope weaker) is still deisotoped
  base[1].setIntensity(0.5);
  s = deisotope(base, 1500.0);
  TEST_EQUAL(s.size(), 1)
  TEST_REAL_SIMILAR(s[0].getIntensity(), 2.57)

  // above the mass limit the old behaviour applies (here: a 2+ cluster of neutral mass 2 * (1000 - proton) = 1997.985)
  MSSpectrum heavy;
  pk.setMZ(1000.0); pk.setIntensity(1.0); heavy.push_back(pk);
  pk.setMZ(1000.0 + 0.5 * Constants::C13C12_MASSDIFF_U); pk.setIntensity(1.2); heavy.push_back(pk);
  s = deisotope(heavy, 1500.0);
  TEST_EQUAL(s.size(), 1)
}
END_SECTION

END_TEST
