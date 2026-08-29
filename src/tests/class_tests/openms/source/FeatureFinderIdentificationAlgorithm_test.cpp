// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderIdentificationAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFinderIdentificationAlgorithm* ptr = 0;
FeatureFinderIdentificationAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderIdentificationAlgorithm())
{
  ptr = new FeatureFinderIdentificationAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFinderIdentificationAlgorithm())
{
  delete ptr;
}
END_SECTION


START_SECTION([EXTRA] run() rejects MS data without spectra)
{
  // Chromatogram extraction indexes the experiment unconditionally, so running on an
  // empty experiment used to read out of bounds instead of reporting the problem
  // (https://github.com/OpenMS/OpenMS/issues/9980).
  FeatureFinderIdentificationAlgorithm algo;

  ProteinIdentification protein;
  protein.setIdentifier("run0");
  vector<ProteinIdentification> proteins{protein};

  PeptideIdentification peptide;
  peptide.setIdentifier("run0");
  peptide.setRT(100.0);
  peptide.setMZ(500.0);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  peptide.setHits({hit});
  PeptideIdentificationList peptides{peptide};

  FeatureMap features;
  FeatureMap seeds;

  TEST_EQUAL(algo.getMSData().empty(), true)
  TEST_EXCEPTION(Exception::IllegalArgument,
                 algo.run(peptides, proteins, PeptideIdentificationList(), vector<ProteinIdentification>(),
                          features, seeds, ""))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



