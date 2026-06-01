// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

START_TEST(ILPDCWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ILPDCWrapper* ptr = nullptr;
ILPDCWrapper* null_ptr = nullptr;
START_SECTION(ILPDCWrapper())
{
	ptr = new ILPDCWrapper();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~ILPDCWrapper())
{
	delete ptr;
}
END_SECTION


START_SECTION((double compute(const FeatureMap fm, PairsType &pairs, Size verbose_level) const))
{
  EmpiricalFormula ef("H1");
  Adduct a(+1, 1, ef.getMonoWeight(), "H1", 0.1, 0, "");
  MassExplainer::AdductsType potential_adducts_;
  potential_adducts_.push_back(a);
  MassExplainer me(potential_adducts_, 1, 3, 2, 0, 0);
  FeatureMap fm;
  ILPDCWrapper::PairsType pairs;

  ILPDCWrapper iw;
  iw.compute(fm, pairs, 1);

  // check that it runs without pairs (i.e. all clusters are singletons)
  TEST_EQUAL(pairs.size(), 0);

  // real data test


}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



