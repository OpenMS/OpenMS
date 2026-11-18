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

  Shuffles are deterministic: they depend only on the seed and the product. rewrite() gives the same product the same
  replacement in every sequence, so that sequences sharing products (e.g. tRNA isodecoders) share their rewritten
  products, as the originals do, and the numbers of distinct original and rewritten products stay balanced.
  Replacements of products of at least min_length residues never repeat.

  rewrite() also checks the products that a search with missed cleavages makes from the rewritten sequences
  (combinations of consecutive products, including ones assembled from products shorter than min_length): if one is
  forbidden, one of its component products gets its next shuffle (everywhere), until none is forbidden.

  Products are compared by key(): the residues (with modifications) without terminal modifications.
*/
class RNABackgroundDatabase
{
public:
  /// @p min_length: products shorter than this are never rewritten (they are not searched)
  RNABackgroundDatabase(const DigestionEnzymeRNA& enzyme, Size min_length, std::uint64_t seed);

  /// The residues of @p oligo without terminal modifications, as a string
  static std::string key(const NASequence& oligo);

  /// Keys of the searchable products of @p sequence: products and combinations of up to @p missed_cleavages + 1
  /// consecutive products, with at least min_length residues
  std::vector<std::string> productKeys(const std::string& sequence, Size missed_cleavages = 0) const;

  /// Whether @p residue stays in place in a shuffle (a cleavage residue of the enzyme)
  bool isFixed(const Ribonucleotide& residue) const;

  /**
    @brief A shuffle of @p product whose key is not forbidden (by @p forbidden or by an earlier replacement of rewrite())

    Draws up to @p max_tries permutations of the movable residues, in an order that depends only on the seed and the
    product, and returns the (@p skip + 1)-th that qualifies. Returns std::nullopt if there is none (e.g. fewer than two
    different movable residues).
  */
  std::optional<NASequence> shuffleProduct(const NASequence& product, const std::function<bool(const std::string&)>& forbidden,
                                           Size max_tries = 200, Size skip = 0) const;

  /**
    @brief @p sequences with products replaced by shuffles

    - Products of at least min_length residues for which @p select is true are replaced by a shuffle that is not
      @p forbidden (shuffleProduct()). A product without a qualifying shuffle stays as it is (counted in failed()).
    - With @p shuffle_short, products shorter than min_length are shuffled as well (they are searched only as part of
      missed-cleavage products, e.g. across a modified G that RNase T1 does not cut).
    - Then the combinations of up to @p missed_cleavages + 1 consecutive products of every rewritten sequence with at
      least min_length residues are checked against @p forbidden. The component products of forbidden ones (also
      products that were not selected) then try other shuffles (up to 20 further qualifying ones; a product shorter
      than min_length may also stay as it is), and keep the one that leaves the fewest forbidden combinations around
      its occurrences, if it improves on the current one. Each product keeps one replacement in all sequences. This
      repeats until nothing improves; forbidden combinations left are counted in failedCombinations().

    A product keeps its replacement across calls. Terminal modifications of the sequences are kept.
  */
  std::vector<std::string> rewrite(const std::vector<std::string>& sequences, const std::function<bool(const std::string&)>& select,
                                   const std::function<bool(const std::string&)>& forbidden, bool shuffle_short = false,
                                   Size missed_cleavages = 0);

  /// rewrite() of a single sequence
  std::string rewrite(const std::string& sequence, const std::function<bool(const std::string&)>& select,
                      const std::function<bool(const std::string&)>& forbidden, bool shuffle_short = false, Size missed_cleavages = 0)
  {
    return rewrite(std::vector<std::string>{sequence}, select, forbidden, shuffle_short, missed_cleavages).front();
  }

  /// Number of distinct products of at least min_length residues that rewrite() replaced
  Size rewritten() const { return rewritten_; }

  /// Number of distinct products of at least min_length residues that rewrite() had to replace but could not
  Size failed() const { return failed_; }

  /// Number of forbidden missed-cleavage combinations left after the last rewrite()
  Size failedCombinations() const { return failed_combinations_; }

  /// Keys of all replacements of products of at least min_length residues made by rewrite()
  const std::set<std::string>& replacements() const { return used_; }

private:
  std::vector<NASequence> products_(const NASequence& sequence) const;

  /// The state of one distinct product in rewrite()
  struct Choice_
  {
    NASequence product;
    bool movable = false;      ///< at least two movable residues
    bool decided = false;      ///< select/shuffle_short evaluated
    bool active = false;       ///< to be replaced
    bool exhausted = false;    ///< no further shuffle qualifies
    Size skip = 0;             ///< which qualifying shuffle (see shuffleProduct())
    std::optional<NASequence> replacement;
  };

  /// (Re)compute the replacement of @p choice from its 'skip'
  void assign_(Choice_& choice, const std::function<bool(const std::string&)>& forbidden);

  RNaseDigestion digestor_;
  std::vector<std::regex> fixed_patterns_;
  Size min_length_;
  std::uint64_t seed_;
  std::map<std::string, Choice_> choices_;
  std::set<std::string> used_;
  Size rewritten_ = 0;
  Size failed_ = 0;
  Size failed_combinations_ = 0;
};
} // namespace OpenMS
