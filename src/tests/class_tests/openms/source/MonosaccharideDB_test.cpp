// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MonosaccharideDB.h>
///////////////////////////

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(MonosaccharideDB, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static const MonosaccharideDB* getInstance()))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  TEST_NOT_EQUAL(db, nullptr)
  // singleton: repeated calls return the same instance
  TEST_EQUAL(db == MonosaccharideDB::getInstance(), true)
}
END_SECTION

START_SECTION((bool hasSymbol(const std::string& symbol) const))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  // common ProForma glycan monosaccharide symbols
  TEST_EQUAL(db->hasSymbol("Hex"), true)
  TEST_EQUAL(db->hasSymbol("HexNAc"), true)
  TEST_EQUAL(db->hasSymbol("Fuc"), true)
  TEST_EQUAL(db->hasSymbol("Neu5Ac"), true)
  // not a known symbol
  TEST_EQUAL(db->hasSymbol("NotAGlycanSymbol"), false)
}
END_SECTION

START_SECTION((const Monosaccharide* getMonosaccharide(const std::string& symbol) const))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  TEST_NOT_EQUAL(db->getMonosaccharide("Hex"), nullptr)
  // an unknown symbol returns nullptr (no throw)
  TEST_EQUAL(db->getMonosaccharide("NotAGlycanSymbol") == nullptr, true)
}
END_SECTION

START_SECTION((const Monosaccharide& getMonosaccharideOrThrow(const std::string& symbol) const))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  // a known symbol returns a reference (no throw)
  const auto& m = db->getMonosaccharideOrThrow("HexNAc");
  TEST_EQUAL(&m == db->getMonosaccharide("HexNAc"), true)
  // an unknown symbol throws
  TEST_EXCEPTION(Exception::ElementNotFound, db->getMonosaccharideOrThrow("NotAGlycanSymbol"))
}
END_SECTION

START_SECTION((std::vector<std::string> getAllSymbols() const))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  std::vector<std::string> syms = db->getAllSymbols();
  TEST_EQUAL(syms.empty(), false)
  TEST_EQUAL(std::find(syms.begin(), syms.end(), std::string("Hex")) != syms.end(), true)
  TEST_EQUAL(std::find(syms.begin(), syms.end(), std::string("HexNAc")) != syms.end(), true)
}
END_SECTION

START_SECTION((Size getNumberOfMonosaccharides() const))
{
  const MonosaccharideDB* db = MonosaccharideDB::getInstance();
  TEST_EQUAL(db->getNumberOfMonosaccharides() > 0, true)
  // consistent with getAllSymbols()
  TEST_EQUAL(db->getNumberOfMonosaccharides(), db->getAllSymbols().size())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
