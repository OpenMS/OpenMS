// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ProForma.h>

namespace OpenMS
{
/**
  @brief A rooted glycan tree, with the reducing-end monosaccharide at node zero.

  Nodes have stable, zero-based indices. Each new node is attached to an existing
  parent, so disconnected graphs and cycles cannot be constructed. Linkage labels
  (for example "beta1-4") are optional metadata; glycosidic fragment masses do not
  depend on them. Monosaccharide formulas describe residues (water excluded), as
  in ProForma. This is not a WURCS parser or a GNOme resolver.

  @ingroup Chemistry
*/
class OPENMS_DLLAPI GlycanStructure
{
public:
  /// A monosaccharide and its linkage towards the reducing end.
  struct OPENMS_DLLAPI Node
  {
    ProForma::GlycanComposition::Monosaccharide monosaccharide;
    std::optional<Size> parent;
    std::string linkage;
  };

  /**
    @brief Add the root or a child and return its stable index.
    @param[in] monosaccharide Named residue or neutral, positive empirical formula
    @param[in] parent No value for the root; an existing node index otherwise
    @param[in] linkage Optional linkage description
    @throws Exception::InvalidParameter for an invalid parent or charged/empty formula
    @throws Exception::ElementNotFound for an unknown monosaccharide
  */
  Size addMonosaccharide(const ProForma::GlycanComposition::Monosaccharide& monosaccharide,
                         std::optional<Size> parent = std::nullopt,
                         const std::string& linkage = "");

  /// Return the nodes in insertion order.
  const std::vector<Node>& getNodes() const;
  /// Return the composition, retaining named and formula components.
  ProForma::GlycanComposition getComposition() const;

private:
  std::vector<Node> nodes_;
};
} // namespace OpenMS
