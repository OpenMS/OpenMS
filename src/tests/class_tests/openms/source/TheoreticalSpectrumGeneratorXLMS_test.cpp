// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <iostream>
#include <cmath>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <algorithm>


START_TEST(TheoreticalSpectrumGeneratorXLMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TheoreticalSpectrumGeneratorXLMS* ptr = nullptr;
TheoreticalSpectrumGeneratorXLMS* nullPointer = nullptr;

/// mostly copied from TheoreticalSpectrumGenerator_test
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
START_SECTION(TheoreticalSpectrumGeneratorXLMS())
  ptr = new TheoreticalSpectrumGeneratorXLMS();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(TheoreticalSpectrumGeneratorXLMS(const TheoreticalSpectrumGeneratorXLMS& source))
  TheoreticalSpectrumGeneratorXLMS copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())

  // getParameters() cannot see the members updateMembers_() caches, so also check
  // that a copy of a configured generator produces the same spectrum
  TheoreticalSpectrumGeneratorXLMS configured;
  Param p(configured.getParameters());
  p.setValue("add_b_ions", "false");
  p.setValue("a_intensity", 0.5);
  configured.setParameters(p);
  TheoreticalSpectrumGeneratorXLMS configured_copy(configured);

  AASequence seq = AASequence::fromString("IFSQVGK");
  PeakSpectrum from_source, from_copy;
  configured.getLinearIonSpectrum(from_source, seq, 3, true, 2);
  configured_copy.getLinearIonSpectrum(from_copy, seq, 3, true, 2);
  TEST_EQUAL(from_copy == from_source, true)
END_SECTION

START_SECTION(~TheoreticalSpectrumGeneratorXLMS())
  delete ptr;
END_SECTION

ptr = new TheoreticalSpectrumGeneratorXLMS();
AASequence peptide = AASequence::fromString("IFSQVGK");

START_SECTION(TheoreticalSpectrumGeneratorXLMS& operator = (const TheoreticalSpectrumGeneratorXLMS& tsg))
  TheoreticalSpectrumGeneratorXLMS copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())

  // same for assignment
  TheoreticalSpectrumGeneratorXLMS configured;
  Param p(configured.getParameters());
  p.setValue("add_b_ions", "false");
  p.setValue("a_intensity", 0.5);
  configured.setParameters(p);
  TheoreticalSpectrumGeneratorXLMS assigned;
  assigned = configured;

  PeakSpectrum from_source, from_assigned;
  configured.getLinearIonSpectrum(from_source, peptide, 3, true, 2);
  assigned.getLinearIonSpectrum(from_assigned, peptide, 3, true, 2);
  TEST_EQUAL(from_assigned == from_source, true)
END_SECTION
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

START_SECTION(virtual void getLinearIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, bool frag_alpha, int charge = 1, Size link_pos_2 = 0))
  PeakSpectrum spec;
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 2);
  TEST_EQUAL(spec.size(), 18)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {43.55185, 57.54930, 74.06004, 86.09642, 102.57077, 114.09134, 117.08605, 131.08351, 147.11280, 152.10497, 160.60207, 174.59953, 204.13426, 233.16484, 261.15975, 303.20268, 320.19686, 348.19178};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);
  TEST_EQUAL(spec.size(), 27)

  spec.clear(true);
  Param param(ptr->getParameters());
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);
  TEST_EQUAL(spec.size(), 54)


//  // test annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);

  // 6 ion types with 3 charges each are expected
  TEST_EQUAL(spec.size(), 30)

  set<std::string> ion_names;
  ion_names.insert("[alpha|ci$b1]");
  ion_names.insert("[alpha|ci$b2]");
  ion_names.insert("[alpha|ci$b2-H2O1]");
  ion_names.insert("[alpha|ci$b3]");
  ion_names.insert("[alpha|ci$b3-H2O1]");
  ion_names.insert("[alpha|ci$b3-H3N1]");
  ion_names.insert("[alpha|ci$x1]");
  ion_names.insert("[alpha|ci$x2]");
  ion_names.insert("[alpha|ci$x3]");
  ion_names.insert("[alpha|ci$x1-H3N1]");
  ion_names.insert("[alpha|ci$x2-H3N1]");
  ion_names.insert("[alpha|ci$x3-H3N1]");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all ion names have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // beta annotations
  spec.clear(true);
  ptr->getLinearIonSpectrum(spec, peptide, 3, false, 3);
  ion_names.clear();
  ion_names.insert("[beta|ci$b1]");
  ion_names.insert("[beta|ci$b2]");
  ion_names.insert("[beta|ci$b2-H2O1]");
  ion_names.insert("[beta|ci$b3]");
  ion_names.insert("[beta|ci$b3-H2O1]");
  ion_names.insert("[beta|ci$b3-H3N1]");
  ion_names.insert("[beta|ci$x1]");
  ion_names.insert("[beta|ci$x2]");
  ion_names.insert("[beta|ci$x3]");
  ion_names.insert("[beta|ci$x1-H3N1]");
  ion_names.insert("[beta|ci$x2-H3N1]");
  ion_names.insert("[beta|ci$x3-H3N1]");

  string_array = spec.getStringDataArrays().at(0);

  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum::IntegerDataArray charge_array = spec.getIntegerDataArrays().at(0);

  int charge_counts[3] = {0, 0, 0};
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 10)
  TEST_EQUAL(charge_counts[1], 10)
  TEST_EQUAL(charge_counts[2], 10)


  param = ptr->getParameters();
  param.setValue("add_losses", "false");
  ptr->setParameters(param);

  // the smallest examples, that make sense for cross-linking
  spec.clear(true);
  AASequence testseq = AASequence::fromString("HA");
  ptr->getLinearIonSpectrum(spec, testseq, 0, true, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getLinearIonSpectrum(spec, testseq, 1, true, 1);
  TEST_EQUAL(spec.size(), 1)

  // loop link
  spec.clear(true);
  testseq = AASequence::fromString("PEPTIDESAREWEIRD");
  ptr->getLinearIonSpectrum(spec, testseq, 1, true, 1, 14);
  TEST_EQUAL(spec.size(), 2)

  spec.clear(true);
  ptr->getLinearIonSpectrum(spec, testseq, 2, false, 1, 14);
  TEST_EQUAL(spec.size(), 3)

  // test isotopic peaks
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 1);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected
  TEST_EQUAL(spec.size(), 18)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2); //
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected, each with a second isotopic peak
  // + a few losses
  TEST_EQUAL(spec.size(), 48)


  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 3); // not supported yet, but it should at least run (with the maximal possible number of peaks)
  ptr->setParameters(param);
  ptr->getLinearIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected, each with a second isotopic peak
  // should be the same result as above for now
  TEST_EQUAL(spec.size(), 48)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  // for quick benchmarking of implementation chances
//  param = ptr->getParameters();
//  param.setValue("add_a_ions", "true");
//  param.setValue("add_b_ions", "true");
//  param.setValue("add_c_ions", "true");
//  param.setValue("add_x_ions", "true");
//  param.setValue("add_y_ions", "true");
//  param.setValue("add_z_ions", "true");
//  param.setValue("add_metainfo", "true");
//  param.setValue("add_losses", "true");
//  ptr->setParameters(param);
//  AASequence tmp_peptide = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
//  for (Size i = 0; i != 1e4; ++i)
//  {
//    PeakSpectrum spec;
//    ptr->getLinearIonSpectrum(spec, tmp_peptide, 9, true, 5);
//  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

END_SECTION

START_SECTION([EXTRA] getLinearIonSpectrum() places suffix-ion neutral-loss peaks correctly at charge >= 2 (CPP-042, issue 10148))
  // Regression test for CPP-042 (issue #10148): the x/y/z branch of addLinearPeaks_ passed the already
  // charge-divided position instead of the charged mono weight to addLinearIonLosses_, which subtracts
  // the loss and divides by the charge again. Suffix loss peaks at charge >= 2 therefore ended up at
  // (M/z - L)/z instead of (M - L)/z. Prefix (a/b/c) ions and cross-linked ions were not affected.
  TheoreticalSpectrumGeneratorXLMS tsg;
  Param p(tsg.getParameters());
  p.setValue("add_losses", "true");
  p.setValue("add_isotopes", "false");
  tsg.setParameters(p);

  // 1) charge independence: every charge-2 peak must have a charge-1 counterpart at 2 * mz - proton.
  //    Link on K (position 3); the suffix part IDESR carries H2O (D, E, S) and NH3 (R) losses.
  AASequence xl_peptide = AASequence::fromString("PEPKIDESR");
  PeakSpectrum spec_z1, spec_z2;
  tsg.getLinearIonSpectrum(spec_z1, xl_peptide, 3, true, 1);
  tsg.getLinearIonSpectrum(spec_z2, xl_peptide, 3, true, 2);
  // per charge: a1-a3, b1-b3, a2/a3/b2/b3-H2O, y1-y5, y1-NH3, y2-y5 -H2O and -NH3 = 24 peaks
  TEST_EQUAL(spec_z1.size(), 24)
  TEST_EQUAL(spec_z2.size(), 48)
  PeakSpectrum::IntegerDataArray charge_array_z2 = spec_z2.getIntegerDataArrays().at(0);
  ABORT_IF(charge_array_z2.size() != spec_z2.size())
  Size checked_z2 = 0;
  for (Size i = 0; i != spec_z2.size(); ++i)
  {
    if (charge_array_z2[i] != 2)
    {
      continue;
    }
    ++checked_z2;
    const double expected_z1 = 2.0 * spec_z2[i].getMZ() - Constants::PROTON_MASS_U;
    bool found = false;
    for (Size j = 0; j != spec_z1.size(); ++j)
    {
      if (std::fabs(spec_z1[j].getMZ() - expected_z1) <= 1e-5)
      {
        found = true;
        break;
      }
    }
    TEST_EQUAL(found, true)
  }
  TEST_EQUAL(checked_z2, 24)

  // 2) concrete value: y1 = S of "KS" (link on K, so only y1 is generated) at charge 2.
  //    M = 2 * 1.00727647 (protons) + 18.01056506 (y offset, H2O) + 87.03202916 (S internal) = 107.05714716
  //    y1(2+) = 53.52857, y1-H2O(2+) = (107.05714716 - 18.01056506) / 2 = 44.52329
  //    The unfixed code emitted the loss peak at (53.52857 - 18.01056506) / 2 = 17.75900 instead.
  AASequence ks = AASequence::fromString("KS");
  PeakSpectrum spec_ks;
  tsg.getLinearIonSpectrum(spec_ks, ks, 0, true, 2);
  TEST_EQUAL(spec_ks.size(), 4) // y1 and y1-H2O at charge 1 and 2
  PeakSpectrum::StringDataArray names_ks = spec_ks.getStringDataArrays().at(0);
  ABORT_IF(names_ks.size() != spec_ks.size())
  bool has_correct_loss = false;
  bool has_wrong_loss = false;
  for (Size i = 0; i != spec_ks.size(); ++i)
  {
    if (std::fabs(spec_ks[i].getMZ() - 44.5233) <= 1e-3)
    {
      has_correct_loss = true;
      TEST_EQUAL(names_ks[i] == "[alpha|ci$y1-H2O1]", true)
    }
    if (std::fabs(spec_ks[i].getMZ() - 17.7590) <= 1e-3)
    {
      has_wrong_loss = true;
    }
  }
  TEST_EQUAL(has_correct_loss, true)
  TEST_EQUAL(has_wrong_loss, false)
END_SECTION

START_SECTION(virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2 = 0))

  // reinitialize TSG to standard parameters
  Param param(ptr->getParameters());
  param.setValue("add_isotopes", "false");
  param.setValue("max_isotope", 2);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_losses", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);

  PeakSpectrum spec;
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);
  TEST_EQUAL(spec.size(), 17)

  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);
  TEST_EQUAL(spec.size(), 17)

  param.setValue("add_metainfo", "false");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);
  TEST_EQUAL(spec.size(), 39)

  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);
  TEST_EQUAL(spec.size(), 39)

  TOLERANCE_ABSOLUTE(0.001)

  param.setValue("add_losses", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);

  double result[] = {442.55421, 551.94577, 566.94214, 580.95645, 599.96494, 618.97210, 629.97925, 661.67042, 661.99842, 663.32768, 667.67394, 827.41502, 849.90957, 870.93103, 899.44378, 927.95451, 944.46524};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 4);
  TEST_EQUAL(spec.size(), 24)

  spec.clear(true);
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 4);
  TEST_EQUAL(spec.size(), 60)


  // test annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_losses", "true");
  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);

  // 6 ion types with 4 charges each are expected
  // + KLinked ions and precursors
  TEST_EQUAL(spec.size(), 75)

  set<std::string> ion_names;
  ion_names.insert("[alpha|xi$b4]");
  ion_names.insert("[alpha|xi$b5]");
  ion_names.insert("[alpha|xi$b6]");
  ion_names.insert("[alpha|xi$x4]");
  ion_names.insert("[alpha|xi$x5]");
  ion_names.insert("[alpha|xi$x6]");
  ion_names.insert("[Q-linked-beta]");
  ion_names.insert("[M+H]");
  ion_names.insert("[M+H]-H2O");
  ion_names.insert("[M+H]-NH3");
  ion_names.insert("[alpha|xi$x4-H3N1]");
  ion_names.insert("[alpha|xi$b4-H2O1]");
  ion_names.insert("[alpha|xi$b4-H3N1]");
  ion_names.insert("[alpha|xi$x5-H2O1]");
  ion_names.insert("[alpha|xi$x5-H3N1]");
  ion_names.insert("[alpha|xi$b5-H2O1]");
  ion_names.insert("[alpha|xi$b5-H3N1]");
  ion_names.insert("[alpha|xi$b6-H3N1]");
  ion_names.insert("[alpha|xi$b6-H2O1]");
  ion_names.insert("[alpha|xi$x6-H3N1]");
  ion_names.insert("[alpha|xi$x6-H2O1]");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all ion names have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // beta annotations
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, false, 2, 4);
  ion_names.clear();
  ion_names.insert("[beta|xi$b4]");
  ion_names.insert("[beta|xi$b5]");
  ion_names.insert("[beta|xi$b6]");
  ion_names.insert("[beta|xi$x4]");
  ion_names.insert("[beta|xi$x5]");
  ion_names.insert("[beta|xi$x6]");
  ion_names.insert("[Q-linked-alpha]");
  ion_names.insert("[M+H]");
  ion_names.insert("[M+H]-H2O");
  ion_names.insert("[M+H]-NH3");
  ion_names.insert("[beta|xi$b6-H2O1]");
  ion_names.insert("[beta|xi$b6-H3N1]");
  ion_names.insert("[beta|xi$x6-H2O1]");
  ion_names.insert("[beta|xi$x6-H3N1]");
  ion_names.insert("[beta|xi$x4-H3N1]");
  ion_names.insert("[beta|xi$b4-H2O1]");
  ion_names.insert("[beta|xi$b4-H3N1]");
  ion_names.insert("[beta|xi$x5-H2O1]");
  ion_names.insert("[beta|xi$x5-H3N1]");
  ion_names.insert("[beta|xi$b5-H2O1]");
  ion_names.insert("[beta|xi$b5-H3N1]");

  string_array = spec.getStringDataArrays().at(0);

  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum::IntegerDataArray charge_array = spec.getIntegerDataArrays().at(0);

  int charge_counts[5] = {0, 0, 0, 0, 0};
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 0)
  TEST_EQUAL(charge_counts[1], 18)
  TEST_EQUAL(charge_counts[2], 18)
  TEST_EQUAL(charge_counts[3], 21)
  TEST_EQUAL(charge_counts[4], 0)

  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  param.setValue("add_losses", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_k_linked_ions", "false");
  ptr->setParameters(param);

  // the smallest examples, that make sense for cross-linking
  spec.clear(true);
  AASequence testseq = AASequence::fromString("HA");
  ptr->getXLinkIonSpectrum(spec, testseq, 0, 2000.0, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, testseq, 1, 2000.0, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  // loop link
  spec.clear(true);
  testseq = AASequence::fromString("PEPTIDESAREWEIRD");
  ptr->getXLinkIonSpectrum(spec, testseq, 1, 2000.0, true, 1, 1, 14);
  TEST_EQUAL(spec.size(), 2)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, testseq, 2, 2000.0, false, 1, 1, 14);
  TEST_EQUAL(spec.size(), 3)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, testseq, 2, 2000.0, false, 1, 1, 13);
  TEST_EQUAL(spec.size(), 4)

  // test isotopic peaks
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 1);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected
  TEST_EQUAL(spec.size(), 24)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2); //
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 3); // not supported yet, but it should at least run (with the maximal possible number of peaks)
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // for quick benchmarking of implementation chances
//  param = ptr->getParameters();
//  param.setValue("add_a_ions", "true");
//  param.setValue("add_b_ions", "true");
//  param.setValue("add_c_ions", "true");
//  param.setValue("add_x_ions", "true");
//  param.setValue("add_y_ions", "true");
//  param.setValue("add_z_ions", "true");
//  param.setValue("add_metainfo", "true");
//  param.setValue("add_losses", "false");
//  ptr->setParameters(param);
//  AASequence tmp_peptide = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
//  for (Size i = 0; i != 1e3; ++i)
//  {
//    PeakSpectrum spec;
//    ptr->getXLinkIonSpectrum(spec, tmp_peptide, 9, 2000.0, false, 2, 5);
//  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
END_SECTION

START_SECTION(virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, int mincharge, int maxcharge))
  // reinitialize TSG to standard parameters
  Param param(ptr->getParameters());
  param.setValue("add_isotopes", "false");
  param.setValue("max_isotope", 2);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_losses", "false");
  param.setValue("add_metainfo", "false");
  param.setValue("add_precursor_peaks", "true");
  param.setValue("add_k_linked_ions", "true");
  ptr->setParameters(param);

  OPXLDataStructs::ProteinProteinCrossLink test_link;
  test_link.alpha = &peptide;
  AASequence beta = AASequence::fromString("TESTPEP");
  test_link.beta = &beta;
  test_link.cross_link_position = std::make_pair<SignedSize, SignedSize> (3, 4);
  test_link.cross_linker_mass = 150.0;

  PeakSpectrum spec;
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 3);
  TEST_EQUAL(spec.size(), 17)

  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 3);
  TEST_EQUAL(spec.size(), 17)

  param.setValue("add_metainfo", "false");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 3);
  TEST_EQUAL(spec.size(), 41)

  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 3);
  TEST_EQUAL(spec.size(), 41)

  TOLERANCE_ABSOLUTE(0.001)

  param.setValue("add_losses", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 3);

  // Example calculation for Residue-Linked Peptide (full peptide with one y/a fragmented residue cross-linked to it)
  // cross-link (with linker mass 150 Da):
  //  IFSQVGK
  //     |
  // TESTPEP

  // left over ion:
  //    yQa
  //     |
  // TESTPEP

  // alpha (M+2H)2+ = 389.72656
  // beta (M+2H)2+ = 380.67165
  // linker 2+ = 75
  //  = 845.39821 - 1 (remove 2/2 to reduce charges from +4 to +2)
  //  precursor mz with charge 2+ = 844.39821

  // IFS b3(2+) without charge protons = 173.59957
  // VGK x3(2+) without charge protons = 164.09466

  // 844.39821 - 173.59957 - 164.09466 =~ 506.70398 (with lazy proton masses)
  // corresponds to 6th ion: 506.71126

  double result[] = {338.14327, 447.53482, 462.53119, 476.54550, 495.55399, 506.71126, 514.56115, 525.56830, 557.25947, 557.58748, 563.26299, 670.79860, 693.29315, 714.31461, 742.82736, 771.33809, 787.84882};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 4);
  TEST_EQUAL(spec.size(), 24)

  spec.clear(true);
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 4);
  TEST_EQUAL(spec.size(), 60)


  // test annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_losses", "true");
  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 5);

  // 6 ion types with 4 charges each are expected
  // + KLinked ions and precursors
  TEST_EQUAL(spec.size(), 79)

  set<std::string> ion_names;
  ion_names.insert("[alpha|xi$b4]");
  ion_names.insert("[alpha|xi$b5]");
  ion_names.insert("[alpha|xi$b6]");
  ion_names.insert("[alpha|xi$x4]");
  ion_names.insert("[alpha|xi$x5]");
  ion_names.insert("[alpha|xi$x6]");
  ion_names.insert("[Q-linked-beta]");
  ion_names.insert("[M+H]");
  ion_names.insert("[M+H]-H2O");
  ion_names.insert("[M+H]-NH3");
  ion_names.insert("[alpha|xi$x4-H3N1]");
  ion_names.insert("[alpha|xi$x4-H2O1]");
  ion_names.insert("[alpha|xi$b4-H2O1]");
  ion_names.insert("[alpha|xi$b4-H3N1]");
  ion_names.insert("[alpha|xi$x5-H2O1]");
  ion_names.insert("[alpha|xi$x5-H3N1]");
  ion_names.insert("[alpha|xi$b5-H2O1]");
  ion_names.insert("[alpha|xi$b5-H3N1]");
  ion_names.insert("[alpha|xi$b6-H3N1]");
  ion_names.insert("[alpha|xi$b6-H2O1]");
  ion_names.insert("[alpha|xi$x6-H3N1]");
  ion_names.insert("[alpha|xi$x6-H2O1]");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all ion names have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    // TEST_EQUAL(name, "TESTSTRING")
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // beta annotations
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link, false, 2, 4);
  ion_names.clear();
  ion_names.insert("[beta|xi$b4]");
  ion_names.insert("[beta|xi$b5]");
  ion_names.insert("[beta|xi$b6]");
  ion_names.insert("[beta|xi$x3]");
  ion_names.insert("[beta|xi$x4]");
  ion_names.insert("[beta|xi$x5]");
  ion_names.insert("[beta|xi$x6]");
  ion_names.insert("[P-linked-alpha]");
  ion_names.insert("[M+H]");
  ion_names.insert("[M+H]-H2O");
  ion_names.insert("[M+H]-NH3");
  ion_names.insert("[beta|xi$x3-H3N1]");
  ion_names.insert("[beta|xi$x3-H2O1]");
  ion_names.insert("[beta|xi$b6-H2O1]");
  ion_names.insert("[beta|xi$b6-H3N1]");
  ion_names.insert("[beta|xi$x6-H2O1]");
  ion_names.insert("[beta|xi$x6-H3N1]");
  ion_names.insert("[beta|xi$x4-H3N1]");
  ion_names.insert("[beta|xi$x4-H2O1]");
  ion_names.insert("[beta|xi$b4-H2O1]");
  ion_names.insert("[beta|xi$b4-H3N1]");
  ion_names.insert("[beta|xi$x5-H2O1]");
  ion_names.insert("[beta|xi$x5-H3N1]");
  ion_names.insert("[beta|xi$b5-H2O1]");
  ion_names.insert("[beta|xi$b5-H3N1]");

  string_array = spec.getStringDataArrays().at(0);

  for (Size i = 0; i != spec.size(); ++i)
  {
    std::string name = string_array[i];
    // TEST_EQUAL(name, "TESTSTRING")
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum::IntegerDataArray charge_array = spec.getIntegerDataArrays().at(0);

  int charge_counts[5] = {0, 0, 0, 0, 0};
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 0)
  TEST_EQUAL(charge_counts[1], 19)
  TEST_EQUAL(charge_counts[2], 19)
  TEST_EQUAL(charge_counts[3], 22)
  TEST_EQUAL(charge_counts[4], 0)

  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  param.setValue("add_losses", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_k_linked_ions", "false");
  ptr->setParameters(param);

  // the smallest examples, that make sense for cross-linking
  spec.clear(true);
  AASequence testseq = AASequence::fromString("HA");

  OPXLDataStructs::ProteinProteinCrossLink test_link_short;
  test_link_short.alpha = &testseq;
  test_link_short.beta = &beta;
  test_link_short.cross_link_position = std::make_pair<SignedSize, SignedSize> (1, 4);
  test_link_short.cross_linker_mass = 150.0;

  ptr->getXLinkIonSpectrum(spec, test_link_short, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, test_link_short, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  // test isotopic peaks
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 1);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 5);
  // 6 ion types with 4 charges each are expected
  TEST_EQUAL(spec.size(), 24)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2); //
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 3); // not supported yet, but it should at least run (with the maximal possible number of peaks)
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, test_link, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

END_SECTION

START_SECTION([EXTRA] precursor isotope peaks are charge normalized)
{
  TheoreticalSpectrumGeneratorXLMS tsg;
  Param param = tsg.getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2);
  param.setValue("add_precursor_peaks", "true");
  param.setValue("add_losses", "false");
  param.setValue("add_metainfo", "true");
  tsg.setParameters(param);

  AASequence xl_peptide = AASequence::fromString("PEPTIDESAREWEIRD");
  PeakSpectrum spec;
  const double precursor_mass = 2000.0;
  const int charge = 3; // precursor peaks are added at the maximal charge
  tsg.getXLinkIonSpectrum(spec, xl_peptide, 3, precursor_mass, true, 2, charge);

  ABORT_IF(spec.getStringDataArrays().empty())
  const auto& names = spec.getStringDataArrays()[0];
  // m/z values of the peaks with the given name, sorted (mono- and isotope peak carry the same name)
  auto mzs_of = [&](const std::string& name)
  {
    std::vector<double> mzs;
    for (Size i = 0; i < spec.size(); ++i)
    {
      if (names[i] == name) mzs.push_back(spec[i].getMZ());
    }
    std::sort(mzs.begin(), mzs.end());
    return mzs;
  };

  TOLERANCE_ABSOLUTE(1e-6)
  // the precursor and its H2O and NH3 losses (added independently of add_losses), each with its first isotope peak,
  // which used to be placed at the charged mass plus a charge divided offset (~2003 instead of ~668)
  const std::vector<std::pair<std::string, double>> expected =
  {
    {"[M+H]", 0.0},
    {"[M+H]-H2O", EmpiricalFormula("H2O").getMonoWeight()},
    {"[M+H]-NH3", EmpiricalFormula("NH3").getMonoWeight()}
  };
  for (const auto& [name, loss] : expected)
  {
    const std::vector<double> mzs = mzs_of(name);
    TEST_EQUAL(mzs.size(), 2)
    ABORT_IF(mzs.size() != 2)
    const double mono_mz = (precursor_mass + charge * Constants::PROTON_MASS_U - loss) / charge;
    TEST_REAL_SIMILAR(mzs[0], mono_mz)
    TEST_REAL_SIMILAR(mzs[1], mono_mz + Constants::C13C12_MASSDIFF_U / charge)
  }
}
END_SECTION


delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
