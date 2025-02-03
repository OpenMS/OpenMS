// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Ribonucleotide, "$Id$")

/////////////////////////////////////////////////////////////

Ribonucleotide* e_ptr = nullptr;
Ribonucleotide* e_nullPointer = nullptr;
START_SECTION((Ribonucleotide()))
  e_ptr = new Ribonucleotide();
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((virtual ~Ribonucleotide()))
  delete e_ptr;
END_SECTION

Ribonucleotide test_ribo("test", "T");

START_SECTION(Ribonucleotide(const Ribonucleotide& ribonucleotide))
  Ribonucleotide copy(test_ribo);
  TEST_EQUAL(copy, test_ribo)
END_SECTION

START_SECTION(Ribonucleotide& operator=(const Ribonucleotide& ribonucleotide))
  Ribonucleotide copy;
  copy = test_ribo;
  TEST_EQUAL(copy, test_ribo);
END_SECTION

START_SECTION(bool operator==(const Ribonucleotide& ribonucleotide) const)
  TEST_EQUAL(Ribonucleotide(),
             Ribonucleotide("unknown ribonucleotide", ".", "", ".", EmpiricalFormula(), '.', 0.0, 0.0, Ribonucleotide::ANYWHERE, EmpiricalFormula("C5H10O5")));
END_SECTION

START_SECTION(String getName() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getName(), "unknown ribonucleotide");
END_SECTION

START_SECTION(String getCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getCode(), ".");
END_SECTION

START_SECTION(String getNewCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getNewCode(), "");
END_SECTION

START_SECTION(String getHTMLCode() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getHTMLCode(), ".");
END_SECTION

START_SECTION(EmpiricalFormula getFormula() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getFormula(), EmpiricalFormula());
END_SECTION

START_SECTION(char getOrigin() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getOrigin(), '.');
END_SECTION

START_SECTION(double getMonoMass() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getMonoMass(), 0.0);
END_SECTION

START_SECTION(double getAvgMass() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getAvgMass(), 0.0);
END_SECTION

START_SECTION(Ribonucleotide::TermSpecificity getTermSpecificity() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getTermSpecificity(), Ribonucleotide::ANYWHERE);
END_SECTION

START_SECTION(EmpiricalFormula getBaselossFormula() const)
  Ribonucleotide empty = Ribonucleotide();
  TEST_EQUAL(empty.getBaselossFormula(), EmpiricalFormula("C5H10O5"));
END_SECTION


START_SECTION(setName(string name))
  Ribonucleotide empty = Ribonucleotide();
  empty.setName("foo");
  TEST_EQUAL(empty.getName(), "foo");
END_SECTION

START_SECTION(setCode(string code))
  Ribonucleotide empty = Ribonucleotide();
  empty.setCode("x");
  TEST_EQUAL(empty.getCode(), "x");
END_SECTION

START_SECTION(setNewCode(String newCode))
  Ribonucleotide empty = Ribonucleotide();
  empty.setNewCode("y");
  TEST_EQUAL(empty.getNewCode(), "y");
END_SECTION

START_SECTION(setHTMLCode(String hmtlCode))
  Ribonucleotide empty = Ribonucleotide();
  empty.setHTMLCode("z");
  TEST_EQUAL(empty.getHTMLCode(), "z");
END_SECTION

START_SECTION(setFormula(EmpiricalFormula formula))
  Ribonucleotide empty = Ribonucleotide();
  empty.setFormula(EmpiricalFormula("H2O"));
  TEST_EQUAL(empty.getFormula(), EmpiricalFormula("H2O"));
END_SECTION

START_SECTION(setOrigin(char origin))
  Ribonucleotide empty = Ribonucleotide();
  empty.setOrigin('q');
  TEST_EQUAL(empty.getOrigin(), 'q');
END_SECTION

START_SECTION(setMonoMass(double mono_mass))
  Ribonucleotide empty = Ribonucleotide();
  empty.setMonoMass(2.0);
  TEST_EQUAL(empty.getMonoMass(), 2.0);
END_SECTION

START_SECTION(setAvgMass(double avg_mass))
  Ribonucleotide empty = Ribonucleotide();
  empty.setAvgMass(3.0);
  TEST_EQUAL(empty.getAvgMass(), 3.0);
END_SECTION

START_SECTION(setTermSpecificity(Ribonucleotide::TermSpecificity specificity))
  Ribonucleotide empty = Ribonucleotide();
  empty.setTermSpecificity(Ribonucleotide::FIVE_PRIME);
  TEST_EQUAL(empty.getTermSpecificity(), Ribonucleotide::FIVE_PRIME);
END_SECTION

START_SECTION(void setBaselossFormula(EmpiricalFormula formula))
  Ribonucleotide empty = Ribonucleotide();
  empty.setBaselossFormula(EmpiricalFormula("CO2"));
  TEST_EQUAL(empty.getBaselossFormula(), EmpiricalFormula("CO2"));
END_SECTION

START_SECTION(bool isModified() const)
  TEST_EQUAL(test_ribo.isModified(), true);
  test_ribo.setOrigin('T');
  TEST_EQUAL(test_ribo.isModified(), false);
  test_ribo.setCode("Tm");
  TEST_EQUAL(test_ribo.isModified(), true);
END_SECTION

END_TEST

