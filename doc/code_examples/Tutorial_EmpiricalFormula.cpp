// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [EmpiricalFormula]

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  EmpiricalFormula methanol("CH3OH"), water("H2O");

  // sum up empirical formula
  EmpiricalFormula sum = methanol + water;

  // get element from ElementDB
  const Element * carbon = ElementDB::getInstance()->getElement("Carbon");

  // output number of carbon atoms and average weight 
  cout << sum << " "
       << sum.getNumberOf(carbon) << " "
       << sum.getAverageWeight() << endl;

  // extract the isotope distribution
  IsotopeDistribution iso_dist = sum.getIsotopeDistribution(CoarseIsotopePatternGenerator(3));

  for (const auto& it : iso_dist)
  {
    cout << it.getMZ() << " " << it.getIntensity() << endl;
  }

  return 0;
} //end of main

//! [EmpiricalFormula]
