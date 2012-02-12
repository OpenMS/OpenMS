// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace ims;
using namespace std;

Weights createWeights()
{
  Map<char, DoubleReal> aa_to_weight;

  set<const Residue*> residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

  for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
  {
    aa_to_weight[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
  }

  // init mass decomposer
  IMSAlphabet alphabet;
  for (Map<char, DoubleReal>::ConstIterator it = aa_to_weight.begin(); it != aa_to_weight.end(); ++it)
  {
    alphabet.push_back(String(it->first), it->second);
  }

  // initializes weights
  Weights weights(alphabet.getMasses(), 0.01);

  // optimize alphabet by dividing by gcd
  weights.divideByGCD();

  return weights;
}

START_TEST(IntegerMassDecomposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IntegerMassDecomposer<>* ptr = 0;
IntegerMassDecomposer<>* null_ptr = 0;

START_SECTION((IntegerMassDecomposer(const Weights &alphabet_)))
{
  ptr = new IntegerMassDecomposer<>(createWeights());
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IntegerMassDecomposer())
{
	delete ptr;
}
END_SECTION

START_SECTION((bool exist(value_type mass)))
{
  // TODO
}
END_SECTION

START_SECTION((IntegerMassDecomposer< ValueType, DecompositionValueType >::decomposition_type getDecomposition(value_type mass)))
{
  // TODO
}
END_SECTION

START_SECTION((IntegerMassDecomposer< ValueType, DecompositionValueType >::decompositions_type getAllDecompositions(value_type mass)))
{
  // TODO
}
END_SECTION

START_SECTION((IntegerMassDecomposer< ValueType, DecompositionValueType >::decomposition_value_type getNumberOfDecompositions(value_type mass)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



