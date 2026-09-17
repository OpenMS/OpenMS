// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/NUXL/NuXLDeisotoper.h>
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

END_TEST
