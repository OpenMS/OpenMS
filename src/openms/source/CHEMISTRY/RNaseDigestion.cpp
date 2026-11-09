// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>

#include <algorithm>
#include <functional>
#include <set>

using namespace std;

namespace OpenMS
{
  void RNaseDigestion::setEnzyme(const DigestionEnzyme* enzyme)
  {
    EnzymaticDigestion::setEnzyme(enzyme);

    const DigestionEnzymeRNA* rnase =
        dynamic_cast<const DigestionEnzymeRNA*>(enzyme_);
    std::string five_prime_code = rnase->getFivePrimeGain();
    if (five_prime_code == "p")
    {
      five_prime_code = "5'-p";
    }
    std::string three_prime_code = rnase->getThreePrimeGain();
    if (three_prime_code == "p")
    {
      three_prime_code = "3'-p";
    } 
    else if (three_prime_code == "c")
    {
      three_prime_code = "3'-c";
    }
    else if (three_prime_code != "")
    {
      three_prime_code = '['+three_prime_code+']';
    }
    

    static const RibonucleotideDB* ribo_db = RibonucleotideDB::getInstance();
    five_prime_gain_ = five_prime_code.empty() ?
                           nullptr :
                           ribo_db->getRibonucleotide(five_prime_code);
    three_prime_gain_ = three_prime_code.empty() ?
                            nullptr :
                            ribo_db->getRibonucleotide(three_prime_code);

    cuts_after_regexes_.clear();
    cuts_before_regexes_.clear();

    StringList CAregexes, CBregexes;
    StringUtils::split(rnase->getCutsAfterRegEx(), ',', CAregexes);
    StringUtils::split(rnase->getCutsBeforeRegEx(), ',', CBregexes);
    for (auto it = std::begin(CAregexes); it != std::end(CAregexes); ++it)
    {
      cuts_after_regexes_.emplace_back(*it);
    }
    for (auto it = std::begin(CBregexes); it != std::end(CBregexes); ++it)
    {
      cuts_before_regexes_.emplace_back(*it);
    }
  }


  void RNaseDigestion::setEnzyme(const std::string& enzyme_name)
  {
    setEnzyme(RNaseDB::getInstance()->getEnzyme(enzyme_name));
  }


  vector<pair<Size, Size>> RNaseDigestion::getFragmentPositions_(
      const NASequence& rna, Size min_length, Size max_length) const
  {
    if (min_length == 0)
      {
        min_length = 1;
      }
    if ((max_length == 0) || (max_length > rna.size()))
    {
      max_length = rna.size();
    }

    vector<pair<Size, Size>> result;
    // no fragment can reach the minimum length (also keeps the unsigned "rna.size() - min_length" below from wrapping)
    if (rna.size() < min_length)
    {
      return result;
    }
    // no length fits (also keeps the unsigned "max_length - min_length + 1" below from wrapping)
    if (max_length < min_length)
    {
      return result;
    }
    if (enzyme_->getName() == NoCleavage) // no cleavage
    {
      Size length = rna.size();
      if ((length >= min_length) && (length <= max_length))
      {
        result.emplace_back(0, length);
      }
    }
    else if (enzyme_->getName() == UnspecificCleavage) // unspecific cleavage
    {
      result.reserve(rna.size() * (max_length - min_length + 1));
      for (Size i = 0; i <= rna.size() - min_length; ++i)
      {
        const Size right = std::min(i + max_length, rna.size());
        for (Size j = i + min_length; j <= right; ++j)
        {
          result.emplace_back(i, j - i);
        }
      }
    }
    else // proper enzyme cleavage
    {
      vector<Size> fragment_pos(1, 0);
      for (Size i = 1; i < rna.size(); ++i)
      {
        bool is_match = true;
        // can't match if we don't have enough bases before or after
        if (i < cuts_after_regexes_.size() || rna.size() - i < cuts_before_regexes_.size())
        {
          is_match = false;
        }
        for (auto it = cuts_after_regexes_.begin(); it != cuts_after_regexes_.end() && is_match; ++it) // Check if the cuts_after_regexes all match
        {
          if (!it->search(rna[i - cuts_after_regexes_.size() + (it - cuts_after_regexes_.begin())]->getCode()))
          {
            is_match = false;
          }
        }
        for (auto it = cuts_before_regexes_.begin(); it != cuts_before_regexes_.end() && is_match; ++it) // Check if the cuts_before_regexes all match
        {
          if (!it->search(rna[i + (it - cuts_before_regexes_.begin())]->getCode()))
          {
            is_match = false;
          }
        }
        if (is_match)
        {
          fragment_pos.push_back(i);
        }
      }
      fragment_pos.push_back(rna.size());

      // "fragment_pos" has at least two elements (zero and "rna.size()"):
      for (Size start_it = 0; start_it < fragment_pos.size() - 1; ++start_it)
      {
        Size start_pos = fragment_pos[start_it];
        for (Size offset = 0; offset <= missed_cleavages_; ++offset)
        {
          Size end_it = start_it + offset + 1;
          if (end_it >= fragment_pos.size())
          {
            break;
          }
          Size end_pos = fragment_pos[end_it];

          Size length = end_pos - start_pos;
          if ((length >= min_length) && (length <= max_length))
          {
            result.emplace_back(start_pos, length);
          }
        }
      }
    }

    return result;
  }

  void RNaseDigestion::digest(const NASequence& rna, vector<NASequence>& output,
                              Size min_length, Size max_length) const
  {
    output.clear();
    vector<DigestionProduct> products;
    digest(rna, products, min_length, max_length);
    output.reserve(products.size());
    for (const auto& product : products)
    {
      output.push_back(product.fragment);
    }
  }


  void RNaseDigestion::digest(const NASequence& rna,
                              vector<DigestionProduct>& output,
                              Size min_length, Size max_length) const
  {
    output.clear();
    if (rna.empty())
    {
      return;
    }

    vector<pair<Size, Size>> positions = getFragmentPositions_(rna, min_length,
                                                               max_length);
    output.reserve(positions.size());
    for (const auto& pos : positions)
    {
      NASequence fragment = rna.getSubsequence(pos.first, pos.second);
      applyTerminalGains_(fragment, pos, rna.size());
      output.push_back({fragment, pos});
    }
  }


  vector<pair<Size, Size>> RNaseDigestion::getFragmentPositions(
    const NASequence& rna, Size min_length, Size max_length) const
  {
    return getFragmentPositions_(rna, min_length, max_length);
  }


  RNaseDigestion::CleavageSensitiveModGroups
  RNaseDigestion::inferCleavageSensitiveMods(
    const set<ConstRibonucleotidePtr>& variable_modifications) const
  {
    CleavageSensitiveModGroups groups;

    for (ConstRibonucleotidePtr mod : variable_modifications)
    {
      if (mod->getTermSpecificity() != Ribonucleotide::ANYWHERE)
      {
        continue;
      }

      std::string origin_code(1, mod->getOrigin());
      const std::string& modified_code = mod->getCode();

      for (const auto& pattern : cuts_after_regexes_)
      {
        bool origin_matches = pattern.search(origin_code);
        bool modified_matches = pattern.search(modified_code);
        if (origin_matches && !modified_matches)
        {
          groups.cuts_after_sensitive.insert(mod);
          break;
        }
      }

      for (const auto& pattern : cuts_before_regexes_)
      {
        bool origin_matches = pattern.search(origin_code);
        bool modified_matches = pattern.search(modified_code);
        if (origin_matches && !modified_matches)
        {
          groups.cuts_before_sensitive.insert(mod);
          break;
        }
      }
    }

    return groups;
  }


  void RNaseDigestion::digestWithCleavageSensitiveMods(
    const NASequence& rna,
    const CleavageSensitiveModGroups& cleavage_sensitive_mods,
    Size max_sensitive_mods_per_fragment,
    vector<DigestionProduct>& output,
    Size min_length,
    Size max_length) const
  {
    output.clear();
    if (rna.empty())
    {
      return;
    }
    if (cleavage_sensitive_mods.combined().empty())
    {
      digest(rna, output, min_length, max_length);
      return;
    }

    RNaseDigestion no_missed = *this;
    no_missed.setMissedCleavages(0);
    const auto pieces = no_missed.getFragmentPositions_(rna, 1, 0);
    vector<Size> cuts{0};
    for (const auto& piece : pieces)
    {
      cuts.push_back(piece.first + piece.second);
    }

    vector<vector<ConstRibonucleotidePtr>> mods_by_pos(rna.size());
    const auto sensitive = cleavage_sensitive_mods.combined();
    for (Size i = 0; i < rna.size(); ++i)
    {
      if (rna[i]->isModified())
      {
        continue;
      }
      const auto& code = rna[i]->getCode();
      for (auto* mod : sensitive)
      {
        if (mod->getTermSpecificity() == Ribonucleotide::ANYWHERE &&
            code.size() == 1 && code[0] == mod->getOrigin())
        {
          mods_by_pos[i].push_back(mod);
        }
      }
    }

    NASequence parent = rna;
    // Test actual cleavage, including multi-residue enzyme patterns. Sensitive
    // modifications can remove cuts, but cannot create new cuts by definition.
    const auto is_cut = [&](Size position)
    {
      if (position == 0 || position == parent.size())
      {
        return true;
      }
      if (position < cuts_after_regexes_.size() ||
          parent.size() - position < cuts_before_regexes_.size()) return false;
      for (Size i = 0; i < cuts_after_regexes_.size(); ++i)
      {
        if (!cuts_after_regexes_[i].search(parent[position - cuts_after_regexes_.size() + i]->getCode()))
        {
          return false;
        }
      }
      for (Size i = 0; i < cuts_before_regexes_.size(); ++i)
      {
        if (!cuts_before_regexes_[i].search(parent[position + i]->getCode()))
        {
          return false;
        }
      }
      return true;
    };
    // A single substituted residue can block at most this many original cuts.
    const Size pattern_width = std::max<Size>(1, cuts_after_regexes_.size() + cuts_before_regexes_.size());
    const Size extra_cuts = std::min(rna.size(), max_sensitive_mods_per_fragment) * pattern_width;
    for (Size left = 0; left + 1 < cuts.size(); ++left)
    {
      for (Size right = left + 1; right < cuts.size(); ++right)
      {
        if (right - left - 1 > missed_cleavages_ &&
            right - left - 1 - missed_cleavages_ > extra_cuts)
        {
          break;
        }
        const Size start = cuts[left], end = cuts[right], length = end - start;
        if (length < std::max<Size>(1, min_length))
        {
          continue;
        }
        if (max_length != 0 && length > max_length)
        {
          break;
        }
        std::function<void(Size, Size)> enumerate = [&](Size next, Size used)
        {
          // A blocked boundary cannot be restored by another blocking mod.
          if (!is_cut(start) || !is_cut(end))
          {
            return;
          }
          Size missed = 0;
          for (Size cut = left + 1; cut < right; ++cut)
          {
            missed += is_cut(cuts[cut]);
          }
          if (missed <= missed_cleavages_)
          {
            NASequence fragment = parent.getSubsequence(start, length);
            applyTerminalGains_(fragment, {start, length}, rna.size());
            output.push_back({fragment, {start, length}});
          }
          if (used >= max_sensitive_mods_per_fragment)
          {
            return;
          }
          for (Size position = next; position < end; ++position)
          {
            const auto* original = parent[position];
            for (auto* mod : mods_by_pos[position])
            {
              parent.set(position, mod);
              enumerate(position + 1, used + 1);
            }
            parent.set(position, original);
          }
        };
        enumerate(start, 0);
      }
    }
  }


  void RNaseDigestion::applyTerminalGains_(NASequence& fragment,
                                           const pair<Size, Size>& pos,
                                           Size parent_size) const
  {
    if ((pos.first > 0) && (five_prime_gain_ != nullptr))
    {
      fragment.setFivePrimeMod(five_prime_gain_);
    }
    if ((pos.first + pos.second < parent_size) &&
        (three_prime_gain_ != nullptr))
    {
      fragment.setThreePrimeMod(three_prime_gain_);
    }
  }


  void RNaseDigestion::DigestedOligo::addEvidence(const IdentificationData::DatabaseSequence& entry, const NASequence& rna,
                                                 const pair<Size, Size>& position)
  {
    using ID = IdentificationData;
    if (sequence_evidence.empty())
      target_decoy = entry.target_decoy;
    else if (target_decoy == ID::TargetDecoy::UNKNOWN || entry.target_decoy == ID::TargetDecoy::UNKNOWN)
      target_decoy = ID::TargetDecoy::UNKNOWN;
    else if (target_decoy != entry.target_decoy)
      target_decoy = ID::TargetDecoy::BOTH;
    const Size end = position.first + position.second;
    ID::SequenceEvidence evidence;
    evidence.database = entry.database;
    evidence.accession = entry.accession;
    evidence.start = static_cast<UInt32>(position.first);
    evidence.end = static_cast<UInt32>(end - 1);
    evidence.before = position.first ? rna[position.first - 1]->getCode()[0] : '[';
    evidence.after = end < rna.size() ? rna[end]->getCode()[0] : ']';
    if (std::find(sequence_evidence.begin(), sequence_evidence.end(), evidence) == sequence_evidence.end())
      sequence_evidence.push_back(std::move(evidence));
  }


  std::vector<RNaseDigestion::DigestedOligo> RNaseDigestion::digest(const IdentificationData::Run& run, Size min_length, Size max_length) const
  {
    using ID = IdentificationData;
    if (run.getMoleculeKind() != ID::MoleculeKind::OLIGONUCLEOTIDE)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RNA digestion requires an oligonucleotide run");
    std::map<NASequence, DigestedOligo> candidates;
    if (run.getDatabaseSequences())
      for (const auto& entry : *run.getDatabaseSequences())
      {
        NASequence rna = NASequence::fromString(entry.sequence);
        for (const auto& pos : getFragmentPositions_(rna, min_length, max_length))
        {
          NASequence fragment = rna.getSubsequence(pos.first, pos.second);
          applyTerminalGains_(fragment, pos, rna.size());
          auto [it, inserted] = candidates.try_emplace(fragment);
          if (inserted) it->second.sequence = fragment;
          it->second.addEvidence(entry, rna, pos);
        }
      }
    std::vector<DigestedOligo> result;
    result.reserve(candidates.size());
    for (auto& [sequence, candidate] : candidates)
      result.push_back(std::move(candidate));
    return result;
  }

} // namespace OpenMS
