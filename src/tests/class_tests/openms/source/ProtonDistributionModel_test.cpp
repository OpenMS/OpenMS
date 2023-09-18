// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(ProtonDistributionModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProtonDistributionModel* ptr = nullptr;
ProtonDistributionModel* nullPointer = nullptr;

START_SECTION(ProtonDistributionModel())
  ptr = new ProtonDistributionModel();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ProtonDistributionModel())
  delete ptr;
END_SECTION

ptr = new ProtonDistributionModel();

START_SECTION(ProtonDistributionModel(const ProtonDistributionModel& model))
  ProtonDistributionModel copy(*ptr);
  NOT_TESTABLE
END_SECTION

START_SECTION(ProtonDistributionModel& operator = (const ProtonDistributionModel& pdm))
  ProtonDistributionModel copy;
  copy = *ptr;
  NOT_TESTABLE
END_SECTION

START_SECTION(void getProtonDistribution(vector<double>& bb_charges, vector<double>& sc_charges, const AASequence& peptide, Int charge, Residue::ResidueType res_type = Residue::YIon))
  vector<double> bb_charges, sc_charges;
  double bb_tmp[] = {1.76496e-09, 2.9459e-13, 6.3724e-12, 2.96724e-13, 0.69332e-13, 6.56286e-13, 4.82365e-13, 3.51139e-13, 5.82514e-23, 1.35049e-12};
  AASequence peptide = AASequence::fromString("DFPIANGER");
  ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);
  for (Size i = 0; i <= peptide.size(); ++i)
  {
    TEST_REAL_SIMILAR(bb_charges[i], bb_tmp[i])
  }

  double sc_tmp[] = {2.7239e-23, 0, 0, 0, 0, 7.77547e-15, 0, 1.15343e-22, 1};
  for (Size i = 0; i != peptide.size(); ++i)
  {
    TEST_REAL_SIMILAR(sc_charges[i], sc_tmp[i])
  }

END_SECTION

START_SECTION((void setPeptideProtonDistribution(const std::vector< double > &bb_charge, const std::vector< double > &sc_charge)))
  vector<double> bb_charges, sc_charges;
  AASequence peptide = AASequence::fromString("DFPIANGER");
  ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);

  ptr->setPeptideProtonDistribution(bb_charges, sc_charges);
  NOT_TESTABLE
END_SECTION

START_SECTION((void getChargeStateIntensities(const AASequence &peptide, const AASequence &n_term_ion, const AASequence &c_term_ion, Int charge, Residue::ResidueType n_term_type, std::vector< double > &n_term_intensities, std::vector< double > &c_term_intensities, FragmentationType type)))
  vector<double> bb_charges, sc_charges;
  AASequence peptide = AASequence::fromString("DFPIANGER");
  ptr->getProtonDistribution(bb_charges, sc_charges, peptide, 1);

  // set the full proton distribution
  ptr->setPeptideProtonDistribution(bb_charges, sc_charges);

  AASequence pre1 = AASequence::fromString("DFP");
  AASequence suf1 = AASequence::fromString("IANGER");
  vector<double> pre_ints, suf_ints;
  ptr->getChargeStateIntensities(peptide, pre1, suf1, 1, Residue::YIon, pre_ints, suf_ints, ProtonDistributionModel::ChargeDirected);

  TEST_EQUAL(pre_ints.size(), 1)
  TEST_EQUAL(suf_ints.size(), 1)
  TEST_REAL_SIMILAR(pre_ints[0], 0.0);
  TEST_REAL_SIMILAR(suf_ints[0], 1.0);

  pre_ints.clear();
  suf_ints.clear();
  ptr->getChargeStateIntensities(peptide, pre1, suf1, 2, Residue::YIon, pre_ints, suf_ints, ProtonDistributionModel::ChargeDirected);
  TEST_EQUAL(pre_ints.size(), 2)
  TEST_EQUAL(suf_ints.size(), 2)
  TOLERANCE_ABSOLUTE(0.01)
  TEST_REAL_SIMILAR(pre_ints[0], 0.40526)
  TEST_REAL_SIMILAR(pre_ints[1], 0.0)
  TEST_REAL_SIMILAR(suf_ints[0], 0.4922)
  TEST_REAL_SIMILAR(suf_ints[1], 0.1025)

END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
