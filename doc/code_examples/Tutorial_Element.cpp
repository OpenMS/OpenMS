// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

//! [doxygen_snippet_Element]

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <iostream>
#include <iomanip>

using namespace OpenMS;
using namespace std;

Int main()
{
  const ElementDB& db = *ElementDB::getInstance();

  // extract carbon element from ElementDB
  // .getResidue("C") would work as well
  const Element& carbon = *db.getElement("Carbon");

  // output name, symbol, monoisotopic weight and average weight
  cout << carbon.getName() << " " << carbon.getSymbol() << " " << carbon.getMonoWeight() << " " << carbon.getAverageWeight() << endl;


  if (db.hasElement("foo")) { std::cout << "worth a try..."; }

  // get all elements currently known; you can also get them by atomic number or symbols:
  const auto all_elements_name = db.getNames();
  const auto all_elements_AN = db.getAtomicNumbers();
  const auto all_elements_symbols = db.getSymbols();
  std::cout << "We currently know of: " << all_elements_name.size() << " elements (incl. isotopes)\n"
            << "                with: " << all_elements_AN.size() << " different atomic numbers (linking to the monoisotopic isotope)\n"
            << "                 and: " << all_elements_symbols.size() << " different symbols\n\n";

  std::cout << "\nLet's find all hydrogen isotopes:\n";
  for (const auto e : all_elements_name)
  {
    // all hydrogens have AN == 1
    if (e.second->getAtomicNumber() == 1)
    {
      std::cout << "  --> " << std::setw(30) << e.first 
                << "      Symbol: " << std::setw(5) << e.second->getSymbol() 
                << "          AN: " << std::setw(3) << e.second->getAtomicNumber()
                << " mono-weight: " << std::setw(14)<< e.second->getMonoWeight() << "\n";
    }
  }

  std::cout << "\nLets print all monoisotopic elements:\n";
  for (const auto e : all_elements_AN)
  {
      std::cout << std::setw(30) << e.first 
                << "      Symbol: " << std::setw(5) << e.second->getSymbol() 
                << "          AN: " << std::setw(3) << e.second->getAtomicNumber()
                << " mono-weight: " << std::setw(14)<< e.second->getMonoWeight() << "\n";
  }


} // end of main

//! [doxygen_snippet_Element]
