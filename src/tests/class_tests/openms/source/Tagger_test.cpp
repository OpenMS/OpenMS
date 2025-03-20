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

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;

START_TEST(Tagger, "$Id$")

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags))

  TheoreticalSpectrumGenerator tsg;
  Param param = tsg.getParameters();
  param.setValue("add_metainfo", "false");
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "true");
  param.setValue("add_losses", "true");
  param.setValue("add_precursor_peaks", "true");
  tsg.setParameters(param);

  // spectrum with charges +1 and +2
  AASequence test_sequence = AASequence::fromString("PEPTIDETESTTHISTAGGER");
  PeakSpectrum spec;
  tsg.getSpectrum(spec, test_sequence, 1, 2);
  TEST_EQUAL(spec.size(), 357);

  std::vector<std::string> tags;

  // tagger searching only for charge +1
  Tagger tagger = Tagger(2, 10, 5, 1, 1);
  tagger.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 890);

  // first aa in prefixes is not recognized yet, unless as false positive
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)

  // last aa in suffixes is not recognized yet, unless as false positive
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), false)

  // tagger searching only for charge +2
  Tagger tagger2 = Tagger(2, 10, 5, 2, 2);
  tags.clear();
  tagger2.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 1006);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  // these are found as false positives with charge +2, in a +1 and +2 spectrum
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charges +1 and +2
  Tagger tagger3 = Tagger(2, 10, 5, 1, 2);
  tags.clear();
  tagger3.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 1094);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)


  // spectrum with charges +1 and +2
  AASequence test_sequence2 = AASequence::fromString("PEPTID(Oxidation)ETESTTHISTAGGER");
  PeakSpectrum spec2;
  tsg.getSpectrum(spec2, test_sequence2, 2, 2);
  TEST_EQUAL(spec2.size(), 180);

  tags.clear();
  tagger3.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 545);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)

  // not found due to modification
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), false)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charge +2 with fixed modification
  Tagger tagger4 = Tagger(2, 10, 5, 2, 2, ListUtils::create<String>("Oxidation (D)"));
  tags.clear();
  tagger4.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 667);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  // modified residue found again
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charge +2 with variable modification
  Tagger tagger5 = Tagger(2, 10, 5, 2, 2, StringList(), ListUtils::create<String>("Oxidation (D)"));
  tags.clear();
  tagger5.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 739);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  // modified residue found again
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // // runtime benchmark, research tags many times in the same spectrum
  // // takes currently about 90 sec
  // std::cout << std::endl;
  // for (int i = 0; i < 5000; i++)
  // {
  //   tags.clear();
  //   tagger3.getTag(spec, tags);
  // }

  // // write out found tags if necessary
  // for (const std::string& tag : tags)
  // {
  //   std::cout << "TEST TAG: " << tag << std::endl;
  // }

END_SECTION

END_TEST
