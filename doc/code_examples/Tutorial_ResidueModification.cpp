//! [doxygen_snippet_ResidueModification]
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // construct a AASequence object, query a residue
  // and output some of its properties
  AASequence aas = AASequence::fromString("DECIANGER");
  cout << aas[2].getName() << " "
       << aas[2].getFormula().toString() << " "
       << aas[2].getModificationName() << " "
       << aas[2].getMonoWeight() << endl;
  
  // find a modification in ModificationsDB
  // and output some of its properties
  // getInstance() returns a pointer to a ModsDB instance
  const ResidueModification* mod = ModificationsDB::getInstance()->getModification("Carbamidomethyl (C)");
  cout << mod->getOrigin() << " "
       << mod->getFullId() << " "
       << mod->getDiffMonoMass() << " "
       << mod->getMonoMass() << endl;
  
  // set the modification on a residue of a peptide
  // and output some of its properties (the formula and mass have changed)
  // in this case ModificationsDB is used in the background
  // to relate the name of the mod to its attributes
  aas.setModification(2, "Carbamidomethyl (C)");
  cout << aas[2].getName() << " "
       << aas[2].getFormula().toString() << " "
       << aas[2].getModificationName() << " "
       << aas[2].getMonoWeight() << endl;

  return 0;
} //end of main

//! [doxygen_snippet_ResidueModification]
