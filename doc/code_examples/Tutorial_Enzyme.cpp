// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_Enzyme]

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
  auto aa_seq = AASequence::fromString("ARCDRE.(Amidated)");
  protease.digest(aa_seq, products);

  // output digestion products
  std::cout << "digesting " << aa_seq.toString() << " into:\n";
  for (const AASequence& p : products)
  {
    cout << "-->  " << p.toString() << "\n";
  }
  cout << endl;

  // allow many miss-cleavages
  protease.setMissedCleavages(10);
  protease.digest(aa_seq, products);

  // output digestion products
  std::cout << "digesting " << aa_seq.toString() << " with 10 MCs into:\n";
  for (const AASequence& p : products)
  {
    cout << "-->  " << p.toString() << "\n";
  }
  cout << endl;

  // verify an infix of a protein is a digestion product:
  String peptide = "FFFRAAA";
  cout << "Is '" << peptide.prefix(4) << "' a valid digestion product of '" << peptide << "'? " 
       << std::boolalpha << protease.isValidProduct(peptide, 0, 4); // yes it is!

}

//! [doxygen_snippet_Enzyme]
