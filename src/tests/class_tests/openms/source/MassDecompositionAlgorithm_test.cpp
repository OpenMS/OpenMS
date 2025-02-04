// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassDecompositionAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassDecompositionAlgorithm* ptr = nullptr;
MassDecompositionAlgorithm* nullPointer = nullptr;
START_SECTION(MassDecompositionAlgorithm())
{
  ptr = new MassDecompositionAlgorithm();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MassDecompositionAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((void getDecompositions(std::vector<MassDecomposition>& decomps, double weight)))
{
  vector<MassDecomposition> decomps;
  double mass = AASequence::fromString("DFPIANGER").getMonoWeight(Residue::Internal);
  cerr << mass << endl;

  MassDecompositionAlgorithm mda;
  Param p(mda.getParameters());
  p.setValue("tolerance", 0.0001);
  mda.setParameters(p);

  mda.getDecompositions(decomps, mass);
  TEST_EQUAL(decomps.size(), 842)

  p.setValue("tolerance", 0.001);
  mda.setParameters(p);
  decomps.clear();
  mda.getDecompositions(decomps, mass);
  TEST_EQUAL(decomps.size(), 911);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



