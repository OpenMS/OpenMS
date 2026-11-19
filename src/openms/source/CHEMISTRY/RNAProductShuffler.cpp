// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RNAProductShuffler.h>

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

RNAProductShuffler::RNAProductShuffler(const DigestionEnzymeRNA& enzyme, Size min_length, std::uint64_t seed) :
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

std::string RNAProductShuffler::key(const NASequence& oligo)
{
  NASequence copy = oligo;
  copy.setFivePrimeMod(nullptr);
  copy.setThreePrimeMod(nullptr);
  return copy.toString();
}

bool RNAProductShuffler::isFixed(const Ribonucleotide& residue) const
{
  for (const auto& regex : fixed_patterns_)
  {
    if (std::regex_search(residue.getCode(), regex)) return true;
  }
  return false;
}

std::vector<NASequence> RNAProductShuffler::products_(const NASequence& sequence) const
{
  std::vector<NASequence> products;
  digestor_.digest(sequence, products);
  return products;
}

std::vector<std::string> RNAProductShuffler::productKeys(const std::string& sequence, Size missed_cleavages) const
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

std::optional<NASequence> RNAProductShuffler::shuffleProduct(const NASequence& product,
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

std::vector<std::string> RNAProductShuffler::failedProducts() const
{
  std::vector<std::string> failed;
  for (const auto& [product_key, choice] : choices_)
  {
    if (choice.active && !choice.replacement && (choice.product.size() >= min_length_)) failed.push_back(product_key);
  }
  return failed;
}

void RNAProductShuffler::assign_(Choice_& choice, const std::function<bool(const std::string&)>& forbidden)
{
  const bool searched_alone = choice.product.size() >= min_length_;
  if (searched_alone && choice.replacement) used_.erase(key(*choice.replacement));
  std::optional<NASequence> shuffled =
    shuffleProduct(choice.product, searched_alone ? forbidden : [](const std::string&) { return false; }, max_tries_, choice.skip);
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

std::vector<std::string> RNAProductShuffler::rewrite(const std::vector<std::string>& sequences,
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

  // Index all searchable windows, including shared originals. Different originals
  // must not acquire the same rewritten key; repeated originals may share one.
  struct Window
  {
    Size sequence, first, last;
    std::string original, rewritten;
  };
  std::vector<Window> windows;
  std::map<std::string, std::vector<Size>> affected_windows;
  std::map<std::string, std::map<std::string, Size>> owners; // rewritten -> original -> occurrence count
  auto text = [&](const std::string& product_key) -> std::string {
    const Choice_& choice = choices_.at(product_key);
    return choice.replacement ? key(*choice.replacement) : product_key;
  };
  for (Size sequence = 0; sequence < parsed.size(); ++sequence)
  {
    const auto& keys = parsed[sequence].keys;
    for (Size first = 0; first < keys.size(); ++first)
    {
      std::string original, rewritten;
      Size length = 0;
      std::set<std::string> components;
      for (Size last = first; last < keys.size() && last - first <= missed_cleavages; ++last)
      {
        original += keys[last];
        rewritten += text(keys[last]);
        length += choices_.at(keys[last]).product.size();
        components.insert(keys[last]);
        if (length < min_length_) { continue; }
        for (const auto& component : components)
        {
          affected_windows[component].push_back(windows.size());
        }
        windows.push_back({sequence, first, last, original, rewritten});
        ++owners[rewritten][original];
      }
    }
  }
  // Update only windows which contain the changed product. The occurrence counts
  // preserve ownership when a shared original occurs in several parent sequences.
  auto update_windows = [&](const std::string& product_key) {
    for (Size index : affected_windows.at(product_key))
    {
      auto& window = windows[index];
      std::string rewritten;
      for (Size position = window.first; position <= window.last; ++position)
      {
        rewritten += text(parsed[window.sequence].keys[position]);
      }
      if (rewritten == window.rewritten) { continue; }
      auto previous = owners.find(window.rewritten);
      auto owner = previous->second.find(window.original);
      if (--owner->second == 0) { previous->second.erase(owner); }
      if (previous->second.empty()) { owners.erase(previous); }
      window.rewritten = std::move(rewritten);
      ++owners[window.rewritten][window.original];
    }
  };
  auto forbidden_window = [&](const Window& window) { return window.first != window.last && forbidden(window.rewritten); };
  auto duplicate_window = [&](const Window& window) { return owners.at(window.rewritten).size() > 1; };
  // Include the rewritten key and conflict kind, so a candidate cannot exchange
  // a forbidden match for a duplicate, or move a collision onto another molecule.
  using Conflict = std::tuple<Size, bool, std::string>;
  auto local_conflicts = [&](const std::string& product_key) {
    std::set<Conflict> conflicts;
    for (Size index : affected_windows.at(product_key))
    {
      const auto& window = windows[index];
      if (forbidden_window(window)) { conflicts.emplace(index, false, window.rewritten); }
      if (duplicate_window(window)) { conflicts.emplace(index, true, window.rewritten); }
    }
    return conflicts;
  };
  const Size candidates_per_product = 20;
  std::set<std::string> settled;
  for (Size round = 0; round < 50; ++round)
  {
    std::set<std::string> components;
    for (const auto& window : windows)
    {
      if (! forbidden_window(window) && ! duplicate_window(window)) { continue; }
      for (Size position = window.first; position <= window.last; ++position)
      {
        const auto& product_key = parsed[window.sequence].keys[position];
        if (choices_.at(product_key).movable && ! settled.count(product_key)) { components.insert(product_key); }
      }
    }
    bool improved = false;
    for (const auto& product_key : components)
    {
      auto& choice = choices_.at(product_key);
      const bool searched_alone = choice.product.size() >= min_length_;
      const auto before = local_conflicts(product_key);
      if (before.empty()) { continue; }
      const auto current = choice.replacement;
      if (searched_alone && current) { used_.erase(key(*current)); }
      auto best = current;
      Size best_count = before.size();
      auto try_candidate = [&](const std::optional<NASequence>& candidate) {
        choice.replacement = candidate;
        update_windows(product_key);
        const auto after = local_conflicts(product_key);
        if (after.size() < best_count && std::includes(before.begin(), before.end(), after.begin(), after.end()))
        {
          best_count = after.size();
          best = candidate;
        }
      };
      auto none = [](const std::string&) { return false; };
      for (Size skip = 0; skip < candidates_per_product; ++skip)
      {
        auto candidate = shuffleProduct(choice.product, searched_alone ? forbidden : none, max_tries_, skip);
        if (! candidate) { break; }
        if (current && key(*candidate) == key(*current)) { continue; }
        try_candidate(candidate);
      }
      if (! searched_alone && current) { try_candidate(std::nullopt); }
      choice.replacement = best;
      update_windows(product_key);
      if (searched_alone && choice.replacement) { used_.insert(key(*choice.replacement)); }
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
    if (! improved) { break; }
    // Changes to shared products can make a previously blocked repair possible.
    settled.clear();
  }
  failed_combinations_
    = std::count_if(windows.begin(), windows.end(), [&](const Window& window) { return forbidden_window(window) || duplicate_window(window); });

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
