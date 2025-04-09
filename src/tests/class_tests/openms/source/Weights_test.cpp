// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <set>
///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/Weights.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_TEST(Weights, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


Weights* ptr = nullptr;
Weights* null_ptr = nullptr;

Weights::alphabet_masses_type masses;
masses.push_back(71.0456);
masses.push_back(180.0312);
masses.push_back(1.0186);
masses.push_back(4284.36894);
masses.push_back(255.0);


START_SECTION(Weights())
{
	ptr = new Weights();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~Weights())
{
	delete ptr;
}
END_SECTION

START_SECTION((Weights(const alphabet_masses_type &masses, alphabet_mass_type prec)))
{
  Weights::alphabet_mass_type precision = 0.01;

  ptr = new Weights(masses, precision);
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((Weights(const Weights &other)))
{
  Weights copy(*ptr);
  TEST_NOT_EQUAL(&copy, null_ptr)

  // test equality of copy and ptr
  ABORT_IF(ptr->size()!=copy.size())

  for(Size i = 0 ; i < ptr->size() ; ++i)
  {
    TEST_EQUAL(ptr->getAlphabetMass(i), copy.getAlphabetMass(i))
    TEST_EQUAL(ptr->getWeight(i), copy.getWeight(i))
    TEST_EQUAL((*ptr)[i], copy[i])
  }
}
END_SECTION

START_SECTION((Weights& operator=(const Weights &weights_)))
{
  Weights copy;
  copy = *ptr;
  TEST_NOT_EQUAL(&copy, null_ptr)

  // test equality of copy and ptr
  ABORT_IF(ptr->size()!=copy.size())

  for(Size i = 0 ; i < ptr->size() ; ++i)
  {
    TEST_EQUAL(ptr->getAlphabetMass(i), copy.getAlphabetMass(i))
    TEST_EQUAL(ptr->getWeight(i), copy.getWeight(i))
    TEST_EQUAL((*ptr)[i], copy[i])
  }
}
END_SECTION

START_SECTION((size_type size() const ))
{
  TEST_EQUAL(ptr->size(), 5)
  Weights w;
  TEST_EQUAL(w.size(),0)
}
END_SECTION

START_SECTION((weight_type getWeight(size_type i) const ))
{
  TEST_EQUAL(ptr->getWeight(0), 7105)
  TEST_EQUAL(ptr->getWeight(1), 18003)
  TEST_EQUAL(ptr->getWeight(2), 102)
  TEST_EQUAL(ptr->getWeight(3), 428437);
  TEST_EQUAL(ptr->getWeight(4), 25500)
}
END_SECTION

START_SECTION((void setPrecision(alphabet_mass_type precision_)))
{
  ptr->setPrecision(0.1);

  TEST_EQUAL(ptr->getWeight(0), 710)
  TEST_EQUAL(ptr->getWeight(1), 1800)
  TEST_EQUAL(ptr->getWeight(2), 10)
  TEST_EQUAL(ptr->getWeight(3), 42844)
  TEST_EQUAL(ptr->getWeight(4), 2550)

  ptr->setPrecision(1);

  TEST_EQUAL(ptr->getWeight(0), 71)
  TEST_EQUAL(ptr->getWeight(1), 180)
  TEST_EQUAL(ptr->getWeight(2), 1)
  TEST_EQUAL(ptr->getWeight(3), 4284)
  TEST_EQUAL(ptr->getWeight(4), 255)

  ptr->setPrecision(0.0001);

  TEST_EQUAL(ptr->getWeight(0), 710456)
  TEST_EQUAL(ptr->getWeight(1), 1800312)
  TEST_EQUAL(ptr->getWeight(2), 10186)
  TEST_EQUAL(ptr->getWeight(3), 42843689)
  TEST_EQUAL(ptr->getWeight(4), 2550000)

  ptr->setPrecision(0.01);

  TEST_EQUAL(ptr->getWeight(0), 7105)
  TEST_EQUAL(ptr->getWeight(1), 18003)
  TEST_EQUAL(ptr->getWeight(2), 102)
  TEST_EQUAL(ptr->getWeight(3), 428437)
  TEST_EQUAL(ptr->getWeight(4), 25500)
}
END_SECTION

START_SECTION((alphabet_mass_type getPrecision() const ))
{
  TEST_EQUAL(ptr->getPrecision(), 0.01)
  ptr->setPrecision(0.00025);
  TEST_EQUAL(ptr->getPrecision(), 0.00025)
  ptr->setPrecision(0.01);
  TEST_EQUAL(ptr->getPrecision(), 0.01)
}
END_SECTION

START_SECTION((weight_type operator[](size_type i) const ))
{
  TEST_EQUAL((*ptr)[0], 7105)
  TEST_EQUAL((*ptr)[1], 18003)
  TEST_EQUAL((*ptr)[2], 102)
  TEST_EQUAL((*ptr)[3], 428437)
  TEST_EQUAL((*ptr)[4], 25500)
}
END_SECTION

START_SECTION((weight_type back() const ))
{
  TEST_EQUAL((*ptr)[4], 25500);
  TEST_EQUAL(ptr->back(), 25500);
}
END_SECTION

START_SECTION((alphabet_mass_type getAlphabetMass(size_type i) const ))
{
  // compare with the masses it was created from
  ABORT_IF(ptr->size() != masses.size())

  for(Size i = 0 ; i < ptr->size() ; ++i)
  {
    TEST_EQUAL(ptr->getAlphabetMass(i), masses[i])
  }
}
END_SECTION

START_SECTION((alphabet_mass_type getParentMass(const std::vector< unsigned int > &decomposition) const ))
{
  vector<unsigned int> base_decomposition(5,0);

  for(Size i = 0 ; i < ptr->size() ; ++i)
  {
    vector<unsigned int> decomposition(base_decomposition);
    decomposition[i] = 1;
    TEST_REAL_SIMILAR(ptr->getParentMass(decomposition), masses[i])
    decomposition[i] = 2;
    TEST_REAL_SIMILAR(ptr->getParentMass(decomposition), 2*masses[i])
  }

  vector<unsigned int> wrong_decomposition(3,0);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ptr->getParentMass(wrong_decomposition),"The passed decomposition has the wrong size. Expected 5 but got 3.")
}
END_SECTION

START_SECTION((void swap(size_type index1, size_type index2)))
{
  Weights copy_to_swap(*ptr);

  copy_to_swap.swap(0,1);
  TEST_EQUAL(ptr->getAlphabetMass(0), copy_to_swap.getAlphabetMass(1))
  TEST_EQUAL(ptr->getAlphabetMass(1), copy_to_swap.getAlphabetMass(0))

  TEST_EQUAL(ptr->getWeight(0), copy_to_swap.getWeight(1))
  TEST_EQUAL(ptr->getWeight(1), copy_to_swap.getWeight(0))

  copy_to_swap.swap(1,3);

  TEST_EQUAL(ptr->getAlphabetMass(0), copy_to_swap.getAlphabetMass(3))
  TEST_EQUAL(ptr->getAlphabetMass(3), copy_to_swap.getAlphabetMass(1))

  TEST_EQUAL(ptr->getWeight(0), copy_to_swap.getWeight(3))
  TEST_EQUAL(ptr->getWeight(3), copy_to_swap.getWeight(1))

}
END_SECTION

START_SECTION((bool divideByGCD()))
{
  // we use the example from the documentation here to demonstrate that
  // it works
  // For example, given alphabet weights 3.0, 5.0, 8.0 with precision 0.1, the
  // integer weights would be 30, 50, 80. After calling this method, the new
  // weights are 3, 5, 8 with precision 1.0 (since the gcd of 30, 50, and 80
  //  is 10).

  Weights::alphabet_masses_type masses_local;
  masses_local.push_back(3.0);
  masses_local.push_back(5.0);
  masses_local.push_back(8.0);

  Weights weights_to_test_gcd(masses_local, 0.1);
  TEST_EQUAL(weights_to_test_gcd.divideByGCD(), true)
  TEST_EQUAL(weights_to_test_gcd[0],3)
  TEST_EQUAL(weights_to_test_gcd[1],5)
  TEST_EQUAL(weights_to_test_gcd[2],8)

  TEST_EQUAL(weights_to_test_gcd.getPrecision(), 1.0)

  // calling it again should not change anything
  TEST_EQUAL(weights_to_test_gcd.divideByGCD(), false)


  Weights::alphabet_masses_type prime_masses;
  prime_masses.push_back(1.13);
  prime_masses.push_back(1.67);
  prime_masses.push_back(2.41);

  Weights prime_weights(prime_masses, 0.01);
  // we cannot find a GCD here
  TEST_EQUAL(prime_weights.divideByGCD(), false)


  Weights::alphabet_masses_type not_enough_masses;
  not_enough_masses.push_back(40);

  Weights not_enough_entries_weights(not_enough_masses,0.01);
  // we cannot divide by GCD if we only have 1 entry
  TEST_EQUAL(not_enough_entries_weights.divideByGCD(), false)
}
END_SECTION

START_SECTION((alphabet_mass_type getMinRoundingError() const ))
{
  TEST_REAL_SIMILAR(ptr->getMinRoundingError(), -6.6655113114361e-06) // for 255.0 -> 25500
  ptr->setPrecision(10);
  TEST_REAL_SIMILAR(ptr->getMinRoundingError(), -1.0) // for 1.0186 -> 0
  ptr->setPrecision(0.01);
}
END_SECTION

START_SECTION((alphabet_mass_type getMaxRoundingError() const ))
{
  TEST_REAL_SIMILAR(ptr->getMaxRoundingError(), 0.00137443549970554)
  ptr->setPrecision(10);
  TEST_REAL_SIMILAR(ptr->getMaxRoundingError(), 0.0196078431372549)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



