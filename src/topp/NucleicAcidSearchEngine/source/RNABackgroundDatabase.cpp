// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <RNABackgroundDatabase.h>

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{
namespace
{
  /// FNV-1a over the seed and a string: a platform-independent hash
  std::uint64_t hash_(std::uint64_t seed, const std::string& text)
  {
    std::uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 8; ++i)
    {
      h ^= (seed >> (8 * i)) & 0xFF;
      h *= 1099511628211ULL;
    }
    for (unsigned char c : text)
    {
      h ^= c;
      h *= 1099511628211ULL;
    }
    return h;
  }

  /// SplitMix64: a small generator with the same output on every platform (unlike the std distributions)
  struct SplitMix64_
  {
    std::uint64_t state;
    std::uint64_t operator()()
    {
      std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
      return z ^ (z >> 31);
    }
  };
} // namespace

RNABackgroundDatabase::RNABackgroundDatabase(const DigestionEnzymeRNA& enzyme, Size min_length, std::uint64_t seed) :
  min_length_(min_length), seed_(seed)
{
  digestor_.setEnzyme(&enzyme);
  digestor_.setMissedCleavages(0);
  StringList patterns, before;
  StringUtils::split(enzyme.getCutsAfterRegEx(), ',', patterns);
  StringUtils::split(enzyme.getCutsBeforeRegEx(), ',', before);
  patterns.insert(patterns.end(), before.begin(), before.end());
  for (const std::string& pattern : patterns)
  {
    if (pattern.empty()) continue;
    std::regex regex(pattern);
    // a pattern that matches every canonical nucleotide (e.g. ".*") does not single out cleavage residues
    bool matches_all = true;
    for (const char* code : {"A", "C", "G", "U"})
    {
      if (!std::regex_search(code, regex)) matches_all = false;
    }
    if (!matches_all) fixed_patterns_.push_back(regex);
  }
}

std::string RNABackgroundDatabase::key(const NASequence& oligo)
{
  NASequence copy = oligo;
  copy.setFivePrimeMod(nullptr);
  copy.setThreePrimeMod(nullptr);
  return copy.toString();
}

bool RNABackgroundDatabase::isFixed(const Ribonucleotide& residue) const
{
  for (const auto& regex : fixed_patterns_)
  {
    if (std::regex_search(residue.getCode(), regex)) return true;
  }
  return false;
}

std::vector<NASequence> RNABackgroundDatabase::products_(const NASequence& sequence) const
{
  std::vector<NASequence> products;
  digestor_.digest(sequence, products);
  return products;
}

std::vector<std::string> RNABackgroundDatabase::productKeys(const std::string& sequence) const
{
  std::vector<std::string> keys;
  for (const NASequence& product : products_(NASequence::fromString(sequence)))
  {
    if (product.size() >= min_length_) keys.push_back(key(product));
  }
  return keys;
}

std::optional<NASequence> RNABackgroundDatabase::shuffleProduct(const NASequence& product,
                                                                const std::function<bool(const std::string&)>& forbidden,
                                                                Size max_tries) const
{
  std::vector<Size> movable;
  for (Size i = 0; i + 1 < product.size(); ++i)
  {
    if (!isFixed(*product[i])) movable.push_back(i);
  }
  if (movable.size() < 2) return std::nullopt;
  const std::string original = key(product);
  SplitMix64_ rng{hash_(seed_, original)};
  for (Size attempt = 0; attempt < max_tries; ++attempt)
  {
    std::vector<const Ribonucleotide*> residues;
    for (Size i : movable) residues.push_back(product[i]);
    for (Size i = residues.size() - 1; i > 0; --i) // Fisher-Yates
    {
      std::swap(residues[i], residues[rng() % (i + 1)]);
    }
    NASequence shuffled = product;
    for (Size i = 0; i < movable.size(); ++i) shuffled[movable[i]] = residues[i];
    const std::string shuffled_key = key(shuffled);
    if ((shuffled_key != original) && !forbidden(shuffled_key) && !used_.count(shuffled_key)) return shuffled;
  }
  return std::nullopt;
}

std::string RNABackgroundDatabase::rewrite(const std::string& sequence, const std::function<bool(const std::string&)>& select,
                                           const std::function<bool(const std::string&)>& forbidden, bool shuffle_short)
{
  const NASequence parent = NASequence::fromString(sequence);
  std::vector<const Ribonucleotide*> residues;
  residues.reserve(parent.size());
  for (const NASequence& product : products_(parent))
  {
    const std::string product_key = key(product);
    const NASequence* replacement = nullptr;
    if ((product.size() >= min_length_) && select(product_key))
    {
      auto pos = replacement_of_.find(product_key);
      if (pos == replacement_of_.end())
      {
        std::optional<NASequence> shuffled = shuffleProduct(product, forbidden);
        if (shuffled)
        {
          used_.insert(key(*shuffled));
          ++rewritten_;
        }
        else
        {
          ++failed_;
        }
        pos = replacement_of_.emplace(product_key, std::move(shuffled)).first;
      }
      if (pos->second) replacement = &*pos->second;
    }
    else if (shuffle_short && (product.size() < min_length_))
    {
      auto pos = short_replacement_of_.find(product_key);
      if (pos == short_replacement_of_.end())
      {
        pos = short_replacement_of_.emplace(product_key, shuffleProduct(product, [](const std::string&) { return false; })).first;
      }
      if (pos->second) replacement = &*pos->second;
    }
    const NASequence& source = replacement ? *replacement : product;
    for (Size i = 0; i < source.size(); ++i) residues.push_back(source[i]);
  }
  return NASequence(residues, parent.getFivePrimeMod(), parent.getThreePrimeMod()).toString();
}
} // namespace OpenMS
