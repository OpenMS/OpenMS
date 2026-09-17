// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/GlycanStructure.h>
#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

START_TEST(GlycanStructure, "$Id$")

START_SECTION((Size addMonosaccharide(const ProForma::GlycanComposition::Monosaccharide&, std::optional<Size>, const std::string&)))
{
  GlycanStructure tree;
  TEST_EXCEPTION(Exception::InvalidParameter, tree.addMonosaccharide(std::string("Hex"), 0))
  TEST_EQUAL(tree.addMonosaccharide(std::string("HexNAc")), 0)
  TEST_EQUAL(tree.addMonosaccharide(std::string("Hex"), 0, "beta1-4"), 1)
  TEST_EQUAL(tree.addMonosaccharide(std::string("Fuc"), 0, "alpha1-6"), 2)
  TEST_FALSE(tree.getNodes()[0].parent.has_value())
  TEST_EQUAL(*tree.getNodes()[2].parent, 0)
  TEST_EQUAL(tree.getNodes()[1].linkage, "beta1-4")
  TEST_EXCEPTION(Exception::InvalidParameter, tree.addMonosaccharide(std::string("Hex")))
  TEST_EXCEPTION(Exception::InvalidParameter, tree.addMonosaccharide(std::string("Hex"), 3))
  TEST_EXCEPTION(Exception::ElementNotFound, tree.addMonosaccharide(std::string("Missing"), 0))
  ProForma::FormulaTag charged {"C6H10O5", 1};
  TEST_EXCEPTION(Exception::InvalidParameter, tree.addMonosaccharide(charged, 0))
  ProForma::FormulaTag negative {"C-1H10O5", std::nullopt};
  TEST_EXCEPTION(Exception::InvalidParameter, tree.addMonosaccharide(negative, 0))
  ProForma::FormulaTag formula {"C6H10O5", std::nullopt};
  TEST_EQUAL(tree.addMonosaccharide(formula, 1), 3)
  TEST_EQUAL(tree.getNodes().size(), 4)
  const auto copy = tree;
  TEST_EQUAL(copy.getComposition().components.size(), 4)
  TEST_EQUAL(std::get<std::string>(copy.getComposition().components[0].first), "HexNAc")
}
END_SECTION

END_TEST
