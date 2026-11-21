// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <RNABackgroundDatabase.h>

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <tuple>

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

std::vector<std::string> RNABackgroundDatabase::productKeys(const std::string& sequence, Size missed_cleavages) const
{
  const std::vector<NASequence> products = products_(NASequence::fromString(sequence));
  std::vector<std::string> keys;
  for (Size i = 0; i < products.size(); ++i)
  {
    std::string joined;
    Size length = 0;
    for (Size j = i; (j < products.size()) && (j <= i + missed_cleavages); ++j)
    {
      joined += key(products[j]);
      length += products[j].size();
      if (length >= min_length_) keys.push_back(joined);
    }
  }
  return keys;
}

std::optional<NASequence> RNABackgroundDatabase::shuffleProduct(const NASequence& product,
                                                                const std::function<bool(const std::string&)>& forbidden,
                                                                Size max_tries, Size skip) const
{
  std::vector<Size> movable;
  for (Size i = 0; i + 1 < product.size(); ++i)
  {
    if (!isFixed(*product[i])) movable.push_back(i);
  }
  if (movable.size() < 2) return std::nullopt;
  const std::string original = key(product);
  SplitMix64_ rng{hash_(seed_, original)};
  std::set<std::string> seen; // a skipped shuffle does not count twice
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
    if ((shuffled_key == original) || forbidden(shuffled_key) || used_.count(shuffled_key) || !seen.insert(shuffled_key).second) continue;
    if (skip == 0) return shuffled;
    --skip;
  }
  return std::nullopt;
}

void RNABackgroundDatabase::assign_(Choice_& choice, const std::function<bool(const std::string&)>& forbidden)
{
  const bool searched_alone = choice.product.size() >= min_length_;
  if (searched_alone && choice.replacement) used_.erase(key(*choice.replacement));
  std::optional<NASequence> shuffled =
    shuffleProduct(choice.product, searched_alone ? forbidden : [](const std::string&) { return false; }, 200, choice.skip);
  if (shuffled)
  {
    choice.replacement = std::move(shuffled);
  }
  else
  {
    choice.exhausted = true; // keep the previous replacement (if any)
  }
  if (searched_alone && choice.replacement) used_.insert(key(*choice.replacement));
}

std::vector<std::string> RNABackgroundDatabase::rewrite(const std::vector<std::string>& sequences,
                                                        const std::function<bool(const std::string&)>& select,
                                                        const std::function<bool(const std::string&)>& forbidden, bool shuffle_short,
                                                        Size missed_cleavages)
{
  struct Parsed
  {
    NASequence parent;
    std::vector<std::string> keys;
  };
  std::vector<Parsed> parsed;
  parsed.reserve(sequences.size());
  for (const std::string& sequence : sequences)
  {
    Parsed p{NASequence::fromString(sequence), {}};
    for (NASequence& product : products_(p.parent))
    {
      std::string product_key = key(product);
      auto pos = choices_.find(product_key);
      if (pos == choices_.end())
      {
        Choice_ choice;
        Size movable = 0;
        for (Size i = 0; i + 1 < product.size(); ++i) movable += !isFixed(*product[i]);
        choice.movable = movable >= 2;
        choice.product = std::move(product);
        choices_.emplace(product_key, std::move(choice));
      }
      p.keys.push_back(std::move(product_key));
    }
    parsed.push_back(std::move(p));
  }

  // first choice of every product, in the order of first appearance
  for (const Parsed& p : parsed)
  {
    for (const std::string& product_key : p.keys)
    {
      Choice_& choice = choices_[product_key];
      if (choice.decided) continue;
      choice.decided = true;
      choice.active = (choice.product.size() >= min_length_) ? select(product_key) : shuffle_short;
      if (choice.active) assign_(choice, forbidden);
    }
  }

  // forbidden missed-cleavage combinations: local search over the replacements of their component products
  std::map<std::string, std::vector<std::pair<Size, Size>>> occurrences; // product -> (sequence, position)
  for (Size s = 0; s < parsed.size(); ++s)
  {
    for (Size i = 0; i < parsed[s].keys.size(); ++i) occurrences[parsed[s].keys[i]].emplace_back(s, i);
  }
  auto text = [&](const std::string& product_key) -> std::string {
    const Choice_& choice = choices_[product_key];
    return choice.replacement ? key(*choice.replacement) : product_key;
  };
  // forbidden combinations (of at least two products) that start at i..j in sequence s
  auto forbidden_from = [&](Size s, Size i, Size j, std::set<std::pair<Size, Size>>* found) {
    const auto& keys = parsed[s].keys;
    Size count = 0;
    std::string joined = text(keys[i]);
    Size length = choices_[keys[i]].product.size();
    for (Size k = i + 1; (k < keys.size()) && (k <= j); ++k)
    {
      joined += text(keys[k]);
      length += choices_[keys[k]].product.size();
      if ((length >= min_length_) && forbidden(joined))
      {
        if (!found || found->emplace(i, k).second) ++count;
      }
    }
    return count;
  };
  // forbidden combinations (sequence, first product, last product) that contain any occurrence of a product
  using Combination = std::tuple<Size, Size, Size>;
  auto local_forbidden = [&](const std::string& product_key) {
    std::set<Combination> found;
    for (const auto& [s, position] : occurrences[product_key])
    {
      const Size first = (position >= missed_cleavages) ? position - missed_cleavages : 0;
      for (Size i = first; i <= position; ++i)
      {
        // combinations i..k with k >= position (those that contain the product)
        const auto& keys = parsed[s].keys;
        std::string joined;
        Size length = 0;
        for (Size k = i; (k < keys.size()) && (k <= i + missed_cleavages); ++k)
        {
          joined += text(keys[k]);
          length += choices_[keys[k]].product.size();
          if ((k > i) && (k >= position) && (length >= min_length_) && forbidden(joined)) found.emplace(s, i, k);
        }
      }
    }
    return found;
  };
  auto total_forbidden = [&]() {
    Size count = 0;
    for (Size s = 0; s < parsed.size(); ++s)
    {
      for (Size i = 0; i < parsed[s].keys.size(); ++i) count += forbidden_from(s, i, i + missed_cleavages, nullptr);
    }
    return count;
  };
  const Size candidates_per_product = 20;
  std::set<std::string> settled; // products whose replacement no candidate improves (until something around them changes)
  for (Size round = 0; round < 50; ++round)
  {
    // the components of forbidden combinations
    std::set<std::string> components;
    for (Size s = 0; s < parsed.size(); ++s)
    {
      const auto& keys = parsed[s].keys;
      for (Size i = 0; i < keys.size(); ++i)
      {
        std::set<std::pair<Size, Size>> found;
        if (forbidden_from(s, i, i + missed_cleavages, &found) == 0) continue;
        for (const auto& [start, end] : found)
        {
          for (Size k = start; k <= end; ++k)
          {
            if (choices_[keys[k]].movable && !settled.count(keys[k])) components.insert(keys[k]);
          }
        }
      }
    }
    bool improved = false;
    for (const std::string& product_key : components)
    {
      Choice_& choice = choices_[product_key];
      const bool searched_alone = choice.product.size() >= min_length_;
      const std::set<Combination> before = local_forbidden(product_key);
      if (before.empty()) continue;
      // candidates: further shuffles; a product not searched on its own may also stay as it is. A candidate is only taken
      // if it removes forbidden combinations without making any new one (trading collisions could move one onto an
      // abundant real molecule)
      std::optional<NASequence> current = choice.replacement;
      if (searched_alone && current) used_.erase(key(*current));
      std::optional<NASequence> best = current;
      Size best_count = before.size();
      auto try_candidate = [&](const std::optional<NASequence>& candidate) {
        choice.replacement = candidate;
        const std::set<Combination> after = local_forbidden(product_key);
        if ((after.size() < best_count) && std::includes(before.begin(), before.end(), after.begin(), after.end()))
        {
          best_count = after.size();
          best = candidate;
        }
      };
      auto none = [](const std::string&) { return false; };
      for (Size skip = 0; skip < candidates_per_product; ++skip)
      {
        std::optional<NASequence> candidate = shuffleProduct(choice.product, searched_alone ? forbidden : none, 200, skip);
        if (!candidate) break;
        if (current && (key(*candidate) == key(*current))) continue;
        try_candidate(candidate);
      }
      if (!searched_alone && current) try_candidate(std::nullopt);
      choice.replacement = best;
      if (searched_alone && choice.replacement) used_.insert(key(*choice.replacement));
      if (best_count < before.size())
      {
        improved = true;
        choice.active = choice.active || bool(choice.replacement);
      }
      else
      {
        settled.insert(product_key);
      }
    }
    if (!improved) break;
  }
  failed_combinations_ = total_forbidden();

  rewritten_ = 0;
  failed_ = 0;
  for (const auto& [product_key, choice] : choices_)
  {
    if (!choice.active || (choice.product.size() < min_length_)) continue;
    if (choice.replacement) ++rewritten_;
    else ++failed_;
  }

  std::vector<std::string> result;
  result.reserve(parsed.size());
  for (const Parsed& p : parsed)
  {
    std::vector<const Ribonucleotide*> residues;
    residues.reserve(p.parent.size());
    for (const std::string& product_key : p.keys)
    {
      const Choice_& choice = choices_[product_key];
      const NASequence& source = choice.replacement ? *choice.replacement : choice.product;
      for (Size i = 0; i < source.size(); ++i) residues.push_back(source[i]);
    }
    result.push_back(NASequence(residues, p.parent.getFivePrimeMod(), p.parent.getThreePrimeMod()).toString());
  }
  return result;
}
} // namespace OpenMS
