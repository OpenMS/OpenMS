// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [Element]

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  const ElementDB * db = ElementDB::getInstance();

  // extract carbon element from ElementDB
  // .getResidue("C") would work as well
  Element carbon = *db->getElement("Carbon"); 

  // output name, symbol, monoisotopic weight and average weight
  cout << carbon.getName() << " "
       << carbon.getSymbol() << " "
       << carbon.getMonoWeight() << " "
       << carbon.getAverageWeight() << endl;

  return 0;
} //end of main

//! [Element]
