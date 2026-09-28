// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/GlycanStructure.h>
#include <OpenMS/CHEMISTRY/MonosaccharideDB.h>

namespace OpenMS
{
Size GlycanStructure::addMonosaccharide(const ProForma::GlycanComposition::Monosaccharide& monosaccharide,
                                        std::optional<Size> parent,
                                        const std::string& linkage)
{
  if ((nodes_.empty() && parent.has_value()) || (! nodes_.empty() && (! parent.has_value() || *parent >= nodes_.size())))
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A glycan must have one root and existing parents");
  }
  EmpiricalFormula formula;
  if (const auto* symbol = std::get_if<std::string>(&monosaccharide))
  {
    formula = EmpiricalFormula(MonosaccharideDB::getInstance()->getMonosaccharideOrThrow(*symbol).formula);
  }
  else
  {
    const auto& tag = std::get<ProForma::FormulaTag>(monosaccharide);
    if (tag.charge.value_or(0) != 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Glycan residues must be neutral");
    }
    formula = EmpiricalFormula(tag.formula_string);
  }
  if (formula.getCharge() != 0 || formula.getNumberOfAtoms() <= 0)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Glycan residues require a nonempty neutral formula");
  }
  for (const auto& atom : formula)
  {
    if (atom.second < 0) { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Negative atom counts in glycan residue"); }
  }
  nodes_.push_back({monosaccharide, parent, linkage});
  return nodes_.size() - 1;
}

const std::vector<GlycanStructure::Node>& GlycanStructure::getNodes() const
{ return nodes_; }

ProForma::GlycanComposition GlycanStructure::getComposition() const
{
  ProForma::GlycanComposition result;
  for (const auto& node : nodes_)
  {
    result.components.emplace_back(node.monosaccharide, 1);
  }
  return result;
}
} // namespace OpenMS
