// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/WeightWrapper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(WeightWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

WeightWrapper* ptr = nullptr;
WeightWrapper* null_ptr = nullptr;
START_SECTION(WeightWrapper())
{
  ptr = new WeightWrapper();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~WeightWrapper())
{
  delete ptr;
}
END_SECTION

START_SECTION((WeightWrapper(const WEIGHTMODE weight_mode)))
{
  WeightWrapper ww(WeightWrapper::MONO);
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::MONO)
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeightMode(), WeightWrapper::AVERAGE)
}
END_SECTION

START_SECTION((WeightWrapper(const WeightWrapper &source)))
{
  WeightWrapper ww(WeightWrapper::AVERAGE);
  WeightWrapper ww2(ww);

  TEST_EQUAL(ww.getWeightMode(), ww2.getWeightMode())
}
END_SECTION

START_SECTION((void setWeightMode(const WEIGHTMODE mode)))
{
  WeightWrapper ww;
  TEST_EXCEPTION(Exception::IllegalArgument, ww.setWeightMode(WeightWrapper::SIZE_OF_WEIGHTMODE))
  ww.setWeightMode(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::AVERAGE)
}
END_SECTION

START_SECTION((WEIGHTMODE getWeightMode() const ))
{
  WeightWrapper ww;
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::MONO)
}
END_SECTION

START_SECTION((double getWeight(const AASequence &aa) const ))
{
  WeightWrapper ww;
  AASequence aa= AASequence::fromString("DFINAGER");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((double getWeight(const EmpiricalFormula &ef) const ))
{
  WeightWrapper ww;
  EmpiricalFormula aa("C12H544");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((double getWeight(const Residue &r, Residue::ResidueType res_type=Residue::Full) const ))
{
  WeightWrapper ww;
  Residue aa("L", "LEU", "L", EmpiricalFormula("C454H33"));
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



