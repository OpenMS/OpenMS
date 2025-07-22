// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/SequenceBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <iostream>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SequenceBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SequenceBase<AASequence>* ptr = nullptr;
SequenceBase<AASequence>* nullPointer = nullptr;

SequenceBase<NASequence>* ptr_na = nullptr;
SequenceBase<NASequence>* nullPointer_na = nullptr;


START_SECTION(SequenceBase()=default)
{
  ptr = new AASequence("PEPTIDE");
  ptr_na = new NASequence("AUCG");
  TEST_NOT_EQUAL(ptr, nullptr);
  TEST_NOT_EQUAL(ptr_na, nullptr);
}
END_SECTION

START_SECTION(Derived& derived())
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(ptr->derived().toString(), seq.toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(ptr_na->derived().toString(), na_seq.toString());
}
END_SECTION

START_SECTION(const Derived& derived() const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(ptr->derived().toString(), seq.toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(ptr_na->derived().toString(), na_seq.toString());
}
END_SECTION

START_SECTION(bool empty() const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(ptr->derived().empty(), seq.empty());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(ptr_na->derived().empty(), na_seq.empty());
}
END_SECTION

START_SECTION(Size size() const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(ptr->derived().size(), seq.size());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(ptr_na->derived().size(), na_seq.size());
}
END_SECTION

START_SECTION(double getMonoWeight(Int charge) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_REAL_SIMILAR(seq.getMonoWeight(), ptr->derived().getMonoWeight());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_REAL_SIMILAR(na_seq.getMonoWeight(), ptr_na->derived().getMonoWeight());
}
END_SECTION

START_SECTION(double getAverageWeight(Int charge) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_REAL_SIMILAR(seq.getAverageWeight(), ptr->derived().getAverageWeight());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_REAL_SIMILAR(na_seq.getAverageWeight(), ptr_na->derived().getAverageWeight());
}
END_SECTION

START_SECTION(EmpiricalFormula getFormula(Int charge) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(seq.getFormula().toString(), ptr->derived().getFormula().toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(na_seq.getFormula().toString(), ptr_na->derived().getFormula().toString());
}
END_SECTION

START_SECTION(Derived getPrefix(Size length) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(seq.getPrefix(3).toString(), ptr->derived().getPrefix(3).toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(na_seq.getPrefix(2).toString(), ptr_na->derived().getPrefix(2).toString());
}
END_SECTION 

START_SECTION(Derived getSuffix(Size length) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(seq.getSuffix(3).toString(), ptr->derived().getSuffix(3).toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(na_seq.getSuffix(2).toString(), ptr_na->derived().getSuffix(2).toString());
}
END_SECTION

START_SECTION(Derived getSubsequence(Size begin, Size end) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(seq.getSubsequence(1, 4).toString(), ptr->derived().getSubsequence(1, 4).toString());
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(na_seq.getSubsequence(0, 2).toString(), ptr_na->derived().getSubsequence(0, 2).toString());
}
END_SECTION

START_SECTION(<typename Func> processSequence(Func func) const)
{
  AASequence seq = AASequence::fromString("PEPTIDE");
  TEST_EQUAL(ptr->processSequence([](const AASequence& aa) { return aa.hasCTerminalModification(); }), false);
  
  NASequence na_seq = NASequence::fromString("AUCG");
  TEST_EQUAL(ptr_na->processSequence([](const NASequence& na) { return na.hasFivePrimeMod(); }), false);
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST