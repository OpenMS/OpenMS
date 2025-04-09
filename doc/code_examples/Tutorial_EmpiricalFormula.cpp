// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_EmpiricalFormula]

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  EmpiricalFormula methanol("CH3OH"), water("H2O");

  // sum up empirical formulae
  EmpiricalFormula sum = methanol + water;

  // get element from ElementDB
  const Element * carbon = ElementDB::getInstance()->getElement("Carbon");

  // output number of carbon atoms and average weight 
  cout << "Formula: " << sum 
       << "\n  average weight: " << sum.getAverageWeight() 
       << "\n  # of Carbons: " << sum.getNumberOf(carbon);

  // extract the isotope distribution
  IsotopeDistribution iso_dist = sum.getIsotopeDistribution(CoarseIsotopePatternGenerator(3));

  std::cout << "\n\nCoarse isotope distribution of " << sum << ": \n";
  for (const auto& it : iso_dist)
  {
    cout << "m/z: " << it.getMZ() << " abundance: " << it.getIntensity() << endl;
  }

} //end of main

//! [doxygen_snippet_EmpiricalFormula]
