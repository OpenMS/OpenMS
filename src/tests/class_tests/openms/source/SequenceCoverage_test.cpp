// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Development Team $
// $Authors: Srikanth K N $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SequenceCoverage.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <OpenMS/CONCEPT/ClassTest.h>

#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(SequenceCoverage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(double getCoverage(const AASequence& protein,
                                 const std::vector<AASequence>& peptides))
{
  AASequence protein = AASequence::fromString("ACDEFGHIK");

  vector<AASequence> peptides;
  peptides.push_back(AASequence::fromString("ACD")); // covers 0–2
  peptides.push_back(AASequence::fromString("FGH")); // covers 4–6

  double coverage = SequenceCoverage::getCoverage(protein, peptides);

  // covered positions: A C D _ F G H _ _
  // total covered = 6 / 9 = 66.666...
  TEST_REAL_SIMILAR(coverage, 66.6667)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([edge case] empty peptide list)
{
  AASequence protein = AASequence::fromString("ACDEFG");

  vector<AASequence> peptides;

  double coverage = SequenceCoverage::getCoverage(protein, peptides);

  TEST_EQUAL(coverage, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([edge case] empty protein)
{
  AASequence protein = AASequence::fromString("");

  vector<AASequence> peptides;
  peptides.push_back(AASequence::fromString("ACD"));

  double coverage = SequenceCoverage::getCoverage(protein, peptides);

  TEST_EQUAL(coverage, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(overlapping peptides)
{
  AASequence protein = AASequence::fromString("ABCDE");

  vector<AASequence> peptides;
  peptides.push_back(AASequence::fromString("ABC")); // 0–2
  peptides.push_back(AASequence::fromString("BCD")); // 1–3

  double coverage = SequenceCoverage::getCoverage(protein, peptides);

  // covered positions: A B C D _
  // 4 / 5 = 80%
  TEST_REAL_SIMILAR(coverage, 80.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
