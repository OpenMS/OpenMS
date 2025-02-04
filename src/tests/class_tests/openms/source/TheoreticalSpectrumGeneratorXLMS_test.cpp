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

  set<String> ion_names;
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
    String name = string_array[i];
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
    String name = string_array[i];
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

  set<String> ion_names;
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
    String name = string_array[i];
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
    String name = string_array[i];
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

  set<String> ion_names;
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
    String name = string_array[i];
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
    String name = string_array[i];
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

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
