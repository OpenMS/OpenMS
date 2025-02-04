// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <map>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

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

START_TEST(RealMassDecomposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RealMassDecomposer* ptr = nullptr;
RealMassDecomposer* null_ptr = nullptr;

START_SECTION((RealMassDecomposer(const Weights &weights)))
{
  ptr = new RealMassDecomposer(createWeights());
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION


START_SECTION(~RealMassDecomposer())
{
	delete ptr;
}
END_SECTION


START_SECTION((decompositions_type getDecompositions(double mass, double error)))
{
  // TODO
}
END_SECTION

START_SECTION((decompositions_type getDecompositions(double mass, double error, const constraints_type &constraints)))
{
  // TODO
}
END_SECTION

START_SECTION((number_of_decompositions_type getNumberOfDecompositions(double mass, double error)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



