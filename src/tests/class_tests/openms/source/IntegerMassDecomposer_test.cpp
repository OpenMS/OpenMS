// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IntegerMassDecomposer.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <map>

using namespace OpenMS;
using namespace ims;
using namespace std;

Weights createWeights()
{
  std::map<char, double> aa_to_weight;

  set<const Residue*> residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

  for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
  {
    aa_to_weight[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
  }

  // init mass decomposer
  IMSAlphabet alphabet;
  for (std::map<char, double>::const_iterator it = aa_to_weight.begin(); it != aa_to_weight.end(); ++it)
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

IntegerMassDecomposer<>* ptr = nullptr;
IntegerMassDecomposer<>* null_ptr = nullptr;

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



