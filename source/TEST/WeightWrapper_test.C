// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/WeightWrapper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(WeightWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

WeightWrapper* ptr = 0;
WeightWrapper* null_ptr = 0;
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

START_SECTION((DoubleReal getWeight(const AASequence &aa) const ))
{
  WeightWrapper ww;
  AASequence aa("DFINAGER");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((DoubleReal getWeight(const EmpiricalFormula &ef) const ))
{
  WeightWrapper ww;
  EmpiricalFormula aa("C12H544");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((DoubleReal getWeight(const Residue &r, Residue::ResidueType res_type=Residue::Full) const ))
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



