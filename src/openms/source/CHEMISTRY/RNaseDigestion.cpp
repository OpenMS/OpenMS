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
    if (rna.empty())
      return;

    vector<pair<Size, Size>> positions = getFragmentPositions_(rna, min_length,
                                                               max_length);

    for (const auto& pos : positions)
    {
      NASequence fragment = rna.getSubsequence(pos.first, pos.second);
      if (pos.first > 0)
      {
        fragment.setFivePrimeMod(five_prime_gain_);
      }
      if (pos.first + pos.second < rna.size())
      {
        fragment.setThreePrimeMod(three_prime_gain_);
      }
      output.push_back(fragment);
    }
  }


  std::vector<RNaseDigestion::DigestedOligo> RNaseDigestion::digest(const IdentificationData::Run& run, Size min_length, Size max_length) const
  {
    using ID = IdentificationData;
    if (run.getMoleculeKind() != ID::MoleculeKind::OLIGONUCLEOTIDE)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RNA digestion requires an oligonucleotide run");
    std::map<NASequence, DigestedOligo> candidates;
    if (run.getDatabaseSequences())
      for (const auto& parent : *run.getDatabaseSequences())
      {
        NASequence rna = NASequence::fromString(parent.sequence);
        for (const auto& pos : getFragmentPositions_(rna, min_length, max_length))
        {
          NASequence fragment = rna.getSubsequence(pos.first, pos.second);
          if (pos.first > 0) fragment.setFivePrimeMod(five_prime_gain_);
          const Size end = pos.first + pos.second;
          if (end < rna.size()) fragment.setThreePrimeMod(three_prime_gain_);
          auto [it, inserted] = candidates.try_emplace(fragment);
          auto& candidate = it->second;
          if (inserted)
          {
            candidate.sequence = fragment;
            candidate.target_decoy = parent.target_decoy;
          }
          else if (candidate.target_decoy == ID::TargetDecoy::UNKNOWN || parent.target_decoy == ID::TargetDecoy::UNKNOWN)
            candidate.target_decoy = ID::TargetDecoy::UNKNOWN;
          else if (candidate.target_decoy != parent.target_decoy)
            candidate.target_decoy = ID::TargetDecoy::BOTH;
          ID::SequenceEvidence evidence;
          evidence.database = parent.database;
          evidence.accession = parent.accession;
          evidence.start = pos.first;
          evidence.end = end - 1;
          evidence.before = std::string(1, pos.first ? rna[pos.first - 1]->getCode()[0] : '[');
          evidence.after = std::string(1, end < rna.size() ? rna[end]->getCode()[0] : ']');
          if (std::find(candidate.sequence_evidence.begin(), candidate.sequence_evidence.end(), evidence) == candidate.sequence_evidence.end())
            candidate.sequence_evidence.push_back(std::move(evidence));
        }
      }
    std::vector<DigestedOligo> result;
    result.reserve(candidates.size());
    for (auto& [sequence, candidate] : candidates)
      result.push_back(std::move(candidate));
    return result;
  }

} // namespace OpenMS
