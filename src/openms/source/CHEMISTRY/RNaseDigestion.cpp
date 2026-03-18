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
    String five_prime_code = rnase->getFivePrimeGain();
    if (five_prime_code == "p")
    {
      five_prime_code = "5'-p";
    }
    String three_prime_code = rnase->getThreePrimeGain();
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
    

    static RibonucleotideDB* ribo_db = RibonucleotideDB::getInstance();
    five_prime_gain_ = five_prime_code.empty() ?
                           nullptr :
                           ribo_db->getRibonucleotide(five_prime_code);
    three_prime_gain_ = three_prime_code.empty() ?
                            nullptr :
                            ribo_db->getRibonucleotide(three_prime_code);

    cuts_after_regexes_.clear();
    cuts_before_regexes_.clear();

    StringList CAregexes, CBregexes;
    rnase->getCutsAfterRegEx().split(',', CAregexes);
    rnase->getCutsBeforeRegEx().split(',', CBregexes);
    for (auto it = std::begin(CAregexes); it != std::end(CAregexes); ++it)
    {
      cuts_after_regexes_.emplace_back(*it);
    }
    for (auto it = std::begin(CBregexes); it != std::end(CBregexes); ++it)
    {
      cuts_before_regexes_.emplace_back(*it);
    }
  }


  void RNaseDigestion::setEnzyme(const String& enzyme_name)
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
          if (!boost::regex_search(rna[i - cuts_after_regexes_.size() + (it - cuts_after_regexes_.begin())]->getCode(), *it))
          {
            is_match = false;
          }
        }
        for (auto it = cuts_before_regexes_.begin(); it != cuts_before_regexes_.end() && is_match; ++it) // Check if the cuts_before_regexes all match
        {
          if (!boost::regex_search(rna[i + (it - cuts_before_regexes_.begin())]->getCode(), *it))
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

      String origin_code(1, mod->getOrigin());
      const String& modified_code = mod->getCode();

      for (const auto& pattern : cuts_after_regexes_)
      {
        bool origin_matches = boost::regex_search(origin_code, pattern);
        bool modified_matches = boost::regex_search(modified_code, pattern);
        if (origin_matches && !modified_matches)
        {
          groups.cuts_after_sensitive.insert(mod);
          break;
        }
      }

      for (const auto& pattern : cuts_before_regexes_)
      {
        bool origin_matches = boost::regex_search(origin_code, pattern);
        bool modified_matches = boost::regex_search(modified_code, pattern);
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

    vector<pair<Size, Size>> base_positions =
      getFragmentPositions_(rna, 1, 0);

    vector<vector<ConstRibonucleotidePtr>> before_mods_by_pos(rna.size());
    vector<vector<ConstRibonucleotidePtr>> after_mods_by_pos(rna.size());
    for (Size i = 0; i < rna.size(); ++i)
    {
      const auto& residue = rna[i];
      if (residue->isModified())
      {
        continue;
      }

      const String& code = residue->getCode();
      if (code.size() != 1)
      {
        continue;
      }

      for (ConstRibonucleotidePtr mod : cleavage_sensitive_mods.cuts_before_sensitive)
      {
        if ((mod->getTermSpecificity() == Ribonucleotide::ANYWHERE) &&
            (code[0] == mod->getOrigin()))
        {
          before_mods_by_pos[i].push_back(mod);
        }
      }
      for (ConstRibonucleotidePtr mod : cleavage_sensitive_mods.cuts_after_sensitive)
      {
        if ((mod->getTermSpecificity() == Ribonucleotide::ANYWHERE) &&
            (code[0] == mod->getOrigin()))
        {
          after_mods_by_pos[i].push_back(mod);
        }
      }
    }

    RNaseDigestion digestor_nomissed = *this;
    digestor_nomissed.setMissedCleavages(0);
    vector<pair<Size, Size>> cut_fragments =
      digestor_nomissed.getFragmentPositions_(rna, 1, 0);
    set<Size> cut_points_set = {0, rna.size()};
    for (const auto& pos : cut_fragments)
    {
      cut_points_set.insert(pos.first);
      cut_points_set.insert(pos.first + pos.second);
    }
    vector<Size> cut_points(cut_points_set.begin(), cut_points_set.end());

    set<String> emitted;
    std::function<void(NASequence&, Size, Size, Size)> recurse =
      [&](NASequence& current_parent, Size start, Size end, Size used_mods)
    {
      const Size length = end - start;
      if ((length >= std::max<Size>(min_length, 1)) &&
          ((max_length == 0) || (length <= max_length)))
      {
        NASequence fragment = current_parent.getSubsequence(start, length);
        applyTerminalGains_(fragment, make_pair(start, length), rna.size());
        String key = String(start);
        key += ":";
        key += String(end);
        key += ":";
        key += fragment.toString();
        if (emitted.insert(key).second)
        {
          output.push_back({fragment, make_pair(start, length)});
        }
      }

      if (used_mods >= max_sensitive_mods_per_fragment)
      {
        return;
      }

      if (start > 0)
      {
        auto cut_it = lower_bound(cut_points.begin(), cut_points.end(), start);
        if (cut_it != cut_points.begin())
        {
          const Size left_cut = *(cut_it - 1);
          if (left_cut < start)
          {
            ConstRibonucleotidePtr original = current_parent[start];
            if (!original->isModified())
            {
              for (ConstRibonucleotidePtr mod : before_mods_by_pos[start])
              {
                current_parent.set(start, mod);
                recurse(current_parent, left_cut, end, used_mods + 1);
              }
              current_parent.set(start, original);
            }
          }
        }
      }

      if ((end > 0) && (end < rna.size()))
      {
        auto cut_it = upper_bound(cut_points.begin(), cut_points.end(), end);
        if (cut_it != cut_points.end())
        {
          const Size right_cut = *cut_it;
          if (right_cut > end)
          {
            const Size boundary_pos = end - 1;
            ConstRibonucleotidePtr original = current_parent[boundary_pos];
            if (!original->isModified())
            {
              for (ConstRibonucleotidePtr mod : after_mods_by_pos[boundary_pos])
              {
                current_parent.set(boundary_pos, mod);
                recurse(current_parent, start, right_cut, used_mods + 1);
              }
              current_parent.set(boundary_pos, original);
            }
          }
        }
      }
    };

    for (const auto& pos : base_positions)
    {
      NASequence current_parent = rna;
      recurse(current_parent, pos.first, pos.first + pos.second, 0);
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


  void RNaseDigestion::digest(IdentificationData& id_data, Size min_length,
                              Size max_length) const
  {
    for (IdentificationData::ParentSequenceRef parent_ref = id_data.getParentSequences().begin();
         parent_ref != id_data.getParentSequences().end(); ++parent_ref)
    {
      if (parent_ref->molecule_type != IdentificationData::MoleculeType::RNA)
      {
        continue;
      }

      NASequence rna = NASequence::fromString(parent_ref->sequence);
      vector<pair<Size, Size>> positions =
          getFragmentPositions_(rna, min_length, max_length);

      for (const auto& pos : positions)
      {
        NASequence fragment = rna.getSubsequence(pos.first, pos.second);
        applyTerminalGains_(fragment, pos, rna.size());
        IdentificationData::IdentifiedOligo oligo(fragment);
        Size end_pos = pos.first + pos.second; // past-the-end position!
        IdentificationData::ParentMatch match(pos.first, end_pos - 1);
        match.left_neighbor = ((pos.first > 0) ?
                               rna[pos.first - 1]->getCode() :
                               IdentificationData::ParentMatch::LEFT_TERMINUS);
        match.right_neighbor = ((end_pos < rna.size()) ?
                                rna[end_pos]->getCode() :
                                IdentificationData::ParentMatch::RIGHT_TERMINUS);
        oligo.parent_matches[parent_ref].insert(match);
        id_data.registerIdentifiedOligo(oligo);
      }
    }
  }

} // namespace OpenMS
