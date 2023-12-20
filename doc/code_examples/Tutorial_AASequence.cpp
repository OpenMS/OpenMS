// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [AASequence]

// This script calculates the mass-to-charge ratio of a 2+ charged b-ion and full peptide from a hardcoded sequence.

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // generate AASequence object from String
  const String s = "DEFIANGER";
  AASequence peptide1 = AASequence::fromString(s);

  // generate AASequence object from string literal
  AASequence peptide2 = AASequence::fromString("PEPTIDER");

  // extract prefix and suffix
  AASequence prefix(peptide1.getPrefix(2));
  AASequence suffix(peptide1.getSuffix(3));
  cout << peptide1.toString() << " "
       << prefix << " "
       << suffix << endl;
  
  // create chemically modified peptide
  AASequence peptide_meth_ox = AASequence::fromString("PEPTIDESEKUEM(Oxidation)CER");
  cout << peptide_meth_ox.toString() << " "
       << peptide_meth_ox.toUnmodifiedString()
       << endl;

  // mass of the full, uncharged peptide
  double peptide_mass_mono = peptide_meth_ox.getMonoWeight();
  cout << "Monoisotopic mass of the uncharged, full peptide: " << peptide_mass_mono << endl;

  double peptide_mass_avg = peptide_meth_ox.getAverageWeight();
  cout << "Average mass of the uncharged, full peptide: " << peptide_mass_avg << endl;

  // mass of the 2+ charged b-ion with the given sequence
  double ion_mass_2plus = peptide_meth_ox.getMonoWeight(Residue::BIon, 2);
  cout << "Mass of the doubly positively charged b-ion: " << ion_mass_2plus << endl;

  // mass-to-charge ratio (m/z) of the 2+ charged b-ion and full peptide with the given sequence
  cout << "Mass-to-charge of the doubly positively charged b-ion: " << peptide_meth_ox.getMZ(2, Residue::BIon) << endl;
  cout << "Mass-to-charge of the doubly positively charged peptide: " << peptide_meth_ox.getMZ(2) << endl;

  // ... many more
  return 0;
}

//! [AASequence]
