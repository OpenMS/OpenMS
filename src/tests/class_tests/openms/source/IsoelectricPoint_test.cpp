// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/IsoelectricPoint.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IsoelectricPoint, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static double computeCharge(const AASequence& seq, double pH, ProteomicsPkaScale scale))
{
  // Empty sequence should throw
  AASequence empty_seq;
  TEST_EXCEPTION(Exception::InvalidValue, IsoelectricPoint::computeCharge(empty_seq, 7.0));

  // At very low pH, charge should be positive (all groups protonated)
  AASequence seq = AASequence::fromString("ACDEFGHIKLMNPQRSTVWY");
  double charge_low = IsoelectricPoint::computeCharge(seq, 1.0);
  TEST_EQUAL(charge_low > 0, true);

  // At very high pH, charge should be negative (all groups deprotonated)
  double charge_high = IsoelectricPoint::computeCharge(seq, 14.0);
  TEST_EQUAL(charge_high < 0, true);

  // For a simple peptide "A" (Alanine): N-term pKa=9.69, C-term pKa=2.34 (Lehninger)
  // At pH 7: charge_nterm = 1/(1+10^(7-9.69)) = ~0.998, charge_cterm = -1/(1+10^(2.34-7)) = ~-1.0
  // Net charge ~ 0.0 (close to zero but slightly negative since pH > pI)
  AASequence ala = AASequence::fromString("A");
  double charge_ala_7 = IsoelectricPoint::computeCharge(ala, 7.0);
  TEST_EQUAL(charge_ala_7 > -0.1, true);
  TEST_EQUAL(charge_ala_7 < 0.1, true);

  // Test that different scales give different results
  double charge_lehninger = IsoelectricPoint::computeCharge(seq, 7.0, ProteomicsPkaScale::LEHNINGER);
  double charge_emboss = IsoelectricPoint::computeCharge(seq, 7.0, ProteomicsPkaScale::EMBOSS);
  TEST_NOT_EQUAL(charge_lehninger, charge_emboss);

  AASequence y = AASequence::fromString("Y");
  AASequence acetyl_y = AASequence::fromString(".(Acetyl)Y");
  TEST_EQUAL(IsoelectricPoint::computeCharge(acetyl_y, 1.0) < IsoelectricPoint::computeCharge(y, 1.0), true);

  AASequence amidated_y = AASequence::fromString("Y.(Amidated)");
  TEST_EQUAL(IsoelectricPoint::computeCharge(amidated_y, 14.0) > IsoelectricPoint::computeCharge(y, 14.0), true);

  AASequence sec = AASequence::fromString("U");
  TEST_EQUAL(IsoelectricPoint::computeCharge(sec, 7.0) < IsoelectricPoint::computeCharge(AASequence::fromString("C"), 7.0), true);

  const auto invalid_enum_value = static_cast<ProteomicsPkaScale>(static_cast<Int>(ProteomicsPkaScale::SIZE_OF_PROTEOMICS_PKA_SCALES));
  TEST_EXCEPTION(Exception::InvalidValue, IsoelectricPoint::computeCharge(ala, 7.0, invalid_enum_value));

  TEST_EXCEPTION(Exception::InvalidValue, IsoelectricPoint::computeCharge(AASequence::fromString("O"), 7.0));
}
END_SECTION

START_SECTION(static double computePI(const AASequence& seq, ProteomicsPkaScale scale, double tolerance))
{
  // Empty sequence should throw
  AASequence empty_seq;
  TEST_EXCEPTION(Exception::InvalidValue, IsoelectricPoint::computePI(empty_seq));

  // Alanine pI should be close to (9.69 + 2.34) / 2 = 6.015 (Lehninger)
  AASequence ala = AASequence::fromString("A");
  double pi_ala = IsoelectricPoint::computePI(ala);
  TEST_REAL_SIMILAR(pi_ala, 6.015);

  // Lysine (K) has an ionizable side chain (pKa = 10.53 Lehninger)
  // pI should be avg of the two highest pKa values: (9.69 + 10.53)/2 = 10.11
  AASequence lys = AASequence::fromString("K");
  double pi_lys = IsoelectricPoint::computePI(lys);
  TEST_REAL_SIMILAR(pi_lys, 10.11);

  // Aspartate (D) has an ionizable side chain (pKa = 3.65 Lehninger)
  // pI should be avg of the two lowest pKa values: (2.34 + 3.65)/2 = 2.995
  AASequence asp = AASequence::fromString("D");
  double pi_asp = IsoelectricPoint::computePI(asp);
  TEST_REAL_SIMILAR(pi_asp, 2.995);

  // At the pI, the charge should be approximately zero
  AASequence peptide = AASequence::fromString("PEPTIDE");
  double pi_peptide = IsoelectricPoint::computePI(peptide);
  double charge_at_pi = IsoelectricPoint::computeCharge(peptide, pi_peptide);
  // computePI() uses a configurable bisection tolerance, so the returned pH only needs to yield a small residual charge.
  TOLERANCE_ABSOLUTE(1e-3)
  TEST_REAL_SIMILAR(charge_at_pi, 0.0);
  TOLERANCE_ABSOLUTE(1e-5)

  // pI must be between 0 and 14
  TEST_EQUAL(pi_peptide > 0.0, true);
  TEST_EQUAL(pi_peptide < 14.0, true);

  // Test EMBOSS scale gives a different pI
  double pi_emboss = IsoelectricPoint::computePI(peptide, ProteomicsPkaScale::EMBOSS);
  double pi_lehninger = IsoelectricPoint::computePI(peptide, ProteomicsPkaScale::LEHNINGER);
  // Both should be reasonable but not identical
  TEST_EQUAL(pi_emboss > 0.0, true);
  TEST_EQUAL(pi_emboss < 14.0, true);
  TEST_NOT_EQUAL(pi_emboss, pi_lehninger);

  AASequence sec = AASequence::fromString("U");
  double pi_sec = IsoelectricPoint::computePI(sec);
  TEST_REAL_SIMILAR(pi_sec, 4.035);

  AASequence polybasic = AASequence::fromString(String(100, 'R'));
  double pi_polybasic = IsoelectricPoint::computePI(polybasic);
  TEST_REAL_SIMILAR(pi_polybasic, 14.0);

  TEST_EXCEPTION(Exception::InvalidValue, IsoelectricPoint::computePI(AASequence::fromString("O")));

  // Bjellqvist scale: pI(A) uses N-term pKa 7.59 (A-specific) and C-term pKa 3.55
  // Expected: (7.59 + 3.55) / 2 = 5.57 (vs Lehninger 6.015)
  double pi_ala_bjellqvist = IsoelectricPoint::computePI(ala, ProteomicsPkaScale::BJELLQVIST);
  TEST_REAL_SIMILAR(pi_ala_bjellqvist, 5.57);
  TEST_EQUAL(std::abs(pi_ala_bjellqvist - pi_ala) > 0.1, true);

  // Bjellqvist scale: pI(P) uses N-term pKa 8.36 (P-specific) and C-term pKa 3.55
  // Expected: (8.36 + 3.55) / 2 = 5.955
  AASequence pro = AASequence::fromString("P");
  double pi_pro_bjellqvist = IsoelectricPoint::computePI(pro, ProteomicsPkaScale::BJELLQVIST);
  TEST_REAL_SIMILAR(pi_pro_bjellqvist, 5.955);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
