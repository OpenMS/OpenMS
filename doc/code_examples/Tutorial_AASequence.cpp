// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_AASequence]

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

  // ... or generate AASequence object from string literal
  AASequence peptide2 = AASequence::fromString("PEPTIDER");

  // extract prefix and suffix of the first/last AA residues
  AASequence prefix(peptide1.getPrefix(2)); // "PE"
  AASequence suffix(peptide1.getSuffix(3)); // "DER"
  cout << peptide1.toString() << " " << prefix << " " << suffix << endl;

  // create chemically modified peptide
  AASequence peptide_meth_ox = AASequence::fromString("PEPTIDESEKUEM(Oxidation)CER");
  cout << peptide_meth_ox.toString() << " --> unmodified: " << peptide_meth_ox.toUnmodifiedString() << endl;

  // mass of the full, uncharged peptide
  double peptide_mass_mono = peptide_meth_ox.getMonoWeight();
  cout << "Monoisotopic mass of the uncharged, full peptide: " << peptide_mass_mono << endl;

  double peptide_mass_avg = peptide_meth_ox.getAverageWeight();
  cout << "Average mass of the uncharged, full peptide: " << peptide_mass_avg << endl;

  // mass of the 2+ charged b-ion with the given sequence
  double ion_mass_b3_2plus = peptide_meth_ox.getPrefix(3).getMonoWeight(Residue::BIon, 2);
  cout << "Mass of the doubly positively charged b3-ion: " << ion_mass_b3_2plus << endl;

  // mass-to-charge ratio (m/z) of the 2+ charged b-ion and full peptide with the given sequence
  cout << "Mass-to-charge of the doubly positively charged b3-ion: " << peptide_meth_ox.getPrefix(3).getMZ(2, Residue::BIon) << endl;
  cout << "Mass-to-charge of the doubly positively charged peptide: " << peptide_meth_ox.getMZ(2) << endl;

  // count AA's to get a frequency table
  std::map<String, Size> aa_freq;
  peptide_meth_ox.getAAFrequencies(aa_freq);
  cout << "Number of Proline (P) residues in '" << peptide_meth_ox.toString() << "' is " << aa_freq['P'] << endl;


  return 0;
}

//! [doxygen_snippet_AASequence]
