// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_Residue]

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

int main()
{
  // get ResidueDB singleton
  ResidueDB const * res_db = ResidueDB::getInstance();

  // query Lysine
  Residue const * lys = res_db->getResidue("Lysine");

  cout << lys->getName() << " "
       << lys->getThreeLetterCode() << " "
       << lys->getOneLetterCode() << " "
       << lys->getFormula().toString() << " "
       << lys->getAverageWeight() << " "
       << lys->getMonoWeight() << endl;

  // one letter code query of Arginine
  Residue const * arg = res_db->getResidue('R');

  cout << arg->getName() << " "
       << arg->getFormula().toString() << " "
       << arg->getMonoWeight() << endl;

  // construct a AASequence object, query a residue 
  // and output some of its properties
  AASequence aas = AASequence::fromString("DEFIANGER");
  cout << aas[3].getName() << " "
       << aas[3].getFormula().toString() << " "
       << aas[3].getMonoWeight() << endl; 

  return 0;
} //end of main

//! [doxygen_snippet_Residue]
