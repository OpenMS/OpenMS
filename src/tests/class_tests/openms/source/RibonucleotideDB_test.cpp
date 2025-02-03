// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(RibonucleotideDB, "$Id$")

/////////////////////////////////////////////////////////////

RibonucleotideDB* ptr = nullptr;
RibonucleotideDB* null = nullptr;
START_SECTION(RibonucleotideDB* getInstance())
{
  ptr = RibonucleotideDB::getInstance();
  TEST_NOT_EQUAL(ptr, null);
}
END_SECTION

START_SECTION(virtual ~RibonucleotideDB())
  NOT_TESTABLE
END_SECTION

START_SECTION(void readFromJSON_(void const std::string& path))
  // Reading from the JSON gets tested as part of the constructor above.
  // We check the contents below in begin() and getRibonucleotide
  NOT_TESTABLE
END_SECTION

START_SECTION(void readFromFile_(void const std::string& path))
  // Reading from the TSV gets tested as part of the constructor above.
  // We check the contents below in getRibonucleotide and getRibonucleotideAlternatives
  NOT_TESTABLE
END_SECTION

START_SECTION(ConstIterator begin())
{
  //Loading of the JSON and TSV files gets tested during the 
  RibonucleotideDB::ConstIterator it = ptr->begin();
  TEST_STRING_EQUAL((*it)->getCode(), "io6A");
}
END_SECTION

START_SECTION(ConstIterator end())
{
  RibonucleotideDB::ConstIterator it = ptr->end();
  TEST_EQUAL(it != ptr->begin(), true);
}
END_SECTION

START_SECTION((const Ribonucleotide& getRibonucleotide(const String& code)))
{
  // These three load from the Modomics.json
  const Ribonucleotide * ribo = ptr->getRibonucleotide("Am");
  TEST_STRING_EQUAL(ribo->getCode(), "Am");
  TEST_STRING_EQUAL(ribo->getName(), "2'-O-methyladenosine");
  // This loads from Custom_RNA_modifications.tsv
  const Ribonucleotide * customRibo = ptr->getRibonucleotide("msU?");
  TEST_STRING_EQUAL(customRibo->getCode(), "msU?");
  TEST_EXCEPTION(Exception::ElementNotFound,
                 ptr->getRibonucleotide("bla"));
}
END_SECTION

START_SECTION( (pair<RibonucleotideDB::ConstRibonucleotidePtr, RibonucleotideDB::ConstRibonucleotidePtr> RibonucleotideDB::getRibonucleotideAlternatives(const std::string& code)))
{
  // THis also tests that loading from the TSV went well
  const pair<RibonucleotideDB::ConstRibonucleotidePtr, RibonucleotideDB::ConstRibonucleotidePtr> alts = ptr->getRibonucleotideAlternatives("msU?");
  TEST_STRING_EQUAL(alts.first->getCode(), "m5s2U");
  TEST_STRING_EQUAL(alts.second->getCode(), "s2Um");
}
END_SECTION

START_SECTION((const Ribonucleotide& getRibonucleotidePrefix(const String& seq)))
{
  const Ribonucleotide* ribo = ptr->getRibonucleotidePrefix("m1AmCGU");
  TEST_STRING_EQUAL(ribo->getCode(), "m1Am");
  TEST_EXCEPTION(Exception::ElementNotFound,
                 ptr->getRibonucleotidePrefix("blam1A"));
}
END_SECTION

START_SECTION(EmpiricalFormula getBaselossFormula())
{
  const Ribonucleotide* dna = ptr->getRibonucleotide("dT");
  TEST_EQUAL(EmpiricalFormula("C5H10O4") == dna->getBaselossFormula(), true);
  const Ribonucleotide* rnam = ptr->getRibonucleotide("Um");
  TEST_EQUAL(EmpiricalFormula("C6H12O5") == rnam->getBaselossFormula(), true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
