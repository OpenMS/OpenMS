// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [Enzyme]

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

#include <vector>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  ProteaseDigestion protease;

  // in this example, we don't produce peptides with missed cleavages
  protease.setMissedCleavages(0);

  // output the number of tryptic peptides (no cut before proline) 
  protease.setEnzyme("Trypsin");
  cout << protease.peptideCount(AASequence::fromString("ACKPDE")) << " "
       << protease.peptideCount(AASequence::fromString("ACRPDEKA"))
       << endl;

  // digest C-terminally amidated peptide 
  vector<AASequence> products;
  protease.digest(AASequence::fromString("ARCDRE.(Amidated)"), products);

  // output digestion products
  for (const AASequence& p : products)
  {
    cout << p.toString() << " ";
  }
  cout << endl;

  // allow many miss-cleavages
  protease.setMissedCleavages(10);
  protease.digest(AASequence::fromString("ARCDRE.(Amidated)"), products);

  // output digestion products
  for (const AASequence& p : products)
  {
    cout << p.toString() << " ";
  }
  cout << endl;

  // ... many more
  return 0;
}

//! [Enzyme]
