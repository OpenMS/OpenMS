// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Per-product shuffling of RNA sequences for NASE's background search ('background:enable').

  Sequences are handled one digestion product at a time (the products of a digestion without missed cleavages, which
  tile the sequence). A shuffle permutes the movable residues of a product: every residue except the last one and
  except the cleavage residues of the enzyme (residues matching a cleavage pattern that does not match every canonical
  nucleotide, e.g. G for RNase T1). The shuffled product therefore digests like the original, and a sequence rewritten
  product by product digests into the rewritten products (and their missed-cleavage combinations).

  Shuffles are deterministic: they depend only on the seed and the product. rewrite() also gives the same product the
  same replacement every time, so that sequences sharing products (e.g. tRNA isodecoders) share their rewritten
  products, as the originals do. Replacements never repeat: every replacement is added to the forbidden keys.

  Products are compared by key(): the residues (with modifications) without terminal modifications.
*/
class RNABackgroundDatabase
{
public:
  /// @p min_length: products shorter than this are never rewritten (they are not searched)
  RNABackgroundDatabase(const DigestionEnzymeRNA& enzyme, Size min_length, std::uint64_t seed);

  /// The residues of @p oligo without terminal modifications, as a string
  static std::string key(const NASequence& oligo);

  /// Keys of the products (no missed cleavages) of @p sequence with at least min_length residues
  std::vector<std::string> productKeys(const std::string& sequence) const;

  /// Whether @p residue stays in place in a shuffle (a cleavage residue of the enzyme)
  bool isFixed(const Ribonucleotide& residue) const;

  /**
    @brief A shuffle of @p product whose key is not forbidden (by @p forbidden or by an earlier replacement of rewrite())

    Tries up to @p max_tries permutations of the movable residues, in an order that depends only on the seed and the
    product. Returns std::nullopt if none qualifies (e.g. fewer than two different movable residues).
  */
  std::optional<NASequence> shuffleProduct(const NASequence& product, const std::function<bool(const std::string&)>& forbidden,
                                           Size max_tries = 200) const;

  /**
    @brief @p sequence with each product of at least min_length residues for which @p select is true replaced by a
    shuffle (shuffleProduct(), avoiding @p forbidden)

    A product gets the same replacement in every sequence (the first one found). A product without a qualifying shuffle
    stays as it is and is counted in failed(). Terminal modifications of @p sequence are kept.

    With @p shuffle_short, products shorter than min_length (which are only searched as part of missed-cleavage
    products, e.g. across a modified G that RNase T1 does not cut) are shuffled as well, without a forbidden check.
  */
  std::string rewrite(const std::string& sequence, const std::function<bool(const std::string&)>& select,
                      const std::function<bool(const std::string&)>& forbidden, bool shuffle_short = false);

  /// Number of distinct products that rewrite() replaced
  Size rewritten() const { return rewritten_; }

  /// Number of distinct products that rewrite() selected but could not replace
  Size failed() const { return failed_; }

  /// Keys of all replacements made by rewrite()
  const std::set<std::string>& replacements() const { return used_; }

private:
  std::vector<NASequence> products_(const NASequence& sequence) const;

  RNaseDigestion digestor_;
  std::vector<std::regex> fixed_patterns_;
  Size min_length_;
  std::uint64_t seed_;
  std::map<std::string, std::optional<NASequence>> replacement_of_;
  std::map<std::string, std::optional<NASequence>> short_replacement_of_;
  std::set<std::string> used_;
  Size rewritten_ = 0;
  Size failed_ = 0;
};
} // namespace OpenMS
