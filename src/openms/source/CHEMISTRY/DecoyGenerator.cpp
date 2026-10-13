// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Sam Wein $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/RegularExpression.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <chrono>
#include <algorithm>
#include <cmath>

using namespace OpenMS;

namespace
{
  NASequence transformRNA(const NASequence& rna, const std::string& rnase,
                         bool shuffle, UInt64 seed, int max_attempts,
                         double identity_threshold, bool keep_five, bool keep_three)
  {
    if (max_attempts < 0 || !std::isfinite(identity_threshold) || identity_threshold < 0 || identity_threshold > 1)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RNA decoy generation requires a non-negative budget and an identity target in [0,1].");
    }
    if (rnase == "unspecific cleavage")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Unspecific RNase products overlap and cannot define a parent-preserving decoy.");
    }
    RNaseDigestion digestor;
    digestor.setEnzyme(rnase);
    digestor.setMissedCleavages(0);
    const auto positions = digestor.getFragmentPositions(rna);
    Size covered = 0;
    for (const auto& [start, length] : positions)
    {
      if (start != covered || start > rna.size() || length == 0 || length > rna.size() - start)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "RNase products must form a contiguous, non-overlapping partition of the RNA parent.");
      }
      covered += length;
    }
    if (covered != rna.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RNase products do not cover the RNA parent.");
    }
    std::vector<bool> fixed(rna.size(), false);
    if (!rna.empty())
    {
      fixed.front() = keep_five;
      fixed.back() = keep_three;
    }
    if (rnase != "no cleavage" && !rna.empty())
    {
      const auto* enzyme = RNaseDB::getInstance()->getEnzyme(rnase);
      StringList after, before;
      StringUtils::split(enzyme->getCutsAfterRegEx(), ',', after);
      StringUtils::split(enzyme->getCutsBeforeRegEx(), ',', before);
      for (const auto& [start, length] : positions)
      {
        if (start == 0) continue;
        // Preserve the complete sequence context responsible for each actual cut.
        for (Size i = 0; i < after.size(); ++i)
        {
          if (i + start >= after.size() && !after[i].empty() && after[i] != ".*") fixed[start - after.size() + i] = true;
        }
        for (Size i = 0; i < before.size() && start + i < rna.size(); ++i)
        {
          if (!before[i].empty() && before[i] != ".*") fixed[start + i] = true;
        }
      }
      const auto preserve_terminal_context = [&](const StringList& patterns, bool at_five)
      {
        if (patterns.empty() || patterns.size() > rna.size() ||
            std::all_of(patterns.begin(), patterns.end(), [](const auto& pattern) { return pattern.empty() || pattern == ".*"; })) return;
        const Size start = at_five ? 0 : rna.size() - patterns.size();
        for (Size i = 0; i < patterns.size(); ++i)
        {
          if (!RegularExpression(patterns[i]).search(rna[start + i]->getCode())) return;
        }
        for (Size i = 0; i < patterns.size(); ++i)
        {
          if (!patterns[i].empty() && patterns[i] != ".*") fixed[start + i] = true;
        }
      };
      // A product used independently may end at an enzyme cut in its original parent.
      preserve_terminal_context(before, true);
      preserve_terminal_context(after, false);
    }
    NASequence result = rna;
    for (const auto& [start, length] : positions)
    {
      const auto fragment = rna.getSubsequence(start, length);
      std::vector<Size> movable;
      std::vector<const Ribonucleotide*> residues;
      for (Size i = 0; i < length; ++i)
      {
        if (!fixed[start + i] && fragment[i]->getTermSpecificity() == Ribonucleotide::ANYWHERE)
        {
          movable.push_back(i);
          residues.push_back(fragment[i]);
        }
      }
      NASequence best = fragment;
      double best_identity = 1.0;
      const int budget = shuffle ? max_attempts : 1;
      Math::RandomShuffler product_shuffler(seed);
      for (int attempt = 0; attempt < budget && movable.size() > 1; ++attempt)
      {
        auto permutation = residues;
        if (shuffle) product_shuffler.portable_random_shuffle(permutation.begin(), permutation.end());
        else std::reverse(permutation.begin(), permutation.end());
        NASequence candidate = fragment;
        for (Size i = 0; i < movable.size(); ++i) candidate[movable[i]] = permutation[i];
        // Context motifs can otherwise create new cuts inside a product.
        const auto cuts = digestor.getFragmentPositions(candidate);
        if (cuts.size() != 1 || cuts[0] != std::make_pair(Size(0), length)) continue;
        const double identity = fragment.getSequenceIdentity(candidate);
        if (!shuffle || identity < best_identity)
        {
          best = candidate;
          best_identity = identity;
        }
        if (!shuffle || best_identity <= identity_threshold) break;
      }
      if (shuffle && max_attempts > 0 && best_identity > identity_threshold)
      {
        OPENMS_LOG_WARN << "RNA decoy product '" << fragment.toString() << "' did not reach the best-effort identity target "
                        << identity_threshold << " within " << max_attempts << " attempts (best " << best_identity
                        << "); composition and modifications are retained.\n";
      }
      for (Size i = 0; i < length; ++i) result[start + i] = best[i];
    }
    if (digestor.getFragmentPositions(result) != positions)
    {
      OPENMS_LOG_WARN << "RNA decoy transformation would change RNase cleavage positions; retaining the original parent.\n";
      return rna;
    }
    return result;
  }
}

NASequence DecoyGenerator::reversePeptides(const NASequence& rna, const std::string& rnase, bool keep_five, bool keep_three) const
{
  return transformRNA(rna, rnase, false, 0, 1, 0.0, keep_five, keep_three);
}

NASequence DecoyGenerator::shufflePeptides(const NASequence& rna, const std::string& rnase, int max_attempts,
                                         double identity_threshold, bool keep_five, bool keep_three)
{
  return transformRNA(rna, rnase, true, seed_, max_attempts, identity_threshold, keep_five, keep_three);
}

DecoyGenerator::DecoyGenerator()
{
  const UInt64 seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  setSeed(seed);
}

void DecoyGenerator::setSeed(UInt64 seed)
{
  seed_ = seed;
  shuffler_.seed(seed);
}

AASequence DecoyGenerator::reverseProtein(const AASequence& protein) const
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.")
  std::string s = protein.toUnmodifiedString();
  std::reverse(s.begin(), s.end());
  return AASequence::fromString(s);
}

AASequence DecoyGenerator::reversePeptides(const AASequence& protein, const std::string& protease) const
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.")
  std::vector<AASequence> peptides;
  ProteaseDigestion ed;
  ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
  ed.setEnzyme(protease);
  ed.setSpecificity(EnzymaticDigestion::SPEC_FULL);
  ed.digest(protein, peptides);    
  std::string pseudo_reversed;
  for (int i = 0; i < static_cast<int>(peptides.size()) - 1; ++i)
  {
    std::string s = peptides[i].toUnmodifiedString();
    auto last = --s.end(); // don't reverse enzymatic cutting site
    std::reverse(s.begin(), last);
    pseudo_reversed += s;
  }
  // the last peptide of a protein is not an enzymatic cutting site so we do a full reverse
  std::string s = peptides[peptides.size() - 1 ].toUnmodifiedString();
  std::reverse(s.begin(), s.end());
  pseudo_reversed += s;
  return AASequence::fromString(pseudo_reversed);
}

// generate decoy protein sequences
std::vector<AASequence> DecoyGenerator::shuffle(const AASequence& protein, const std::string& protease, int decoy_factor)
{
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.");
  
  ProteaseDigestion digestor;
  digestor.setEnzyme(protease);
  digestor.setMissedCleavages(0);  // for decoy generation disable missed cleavages
  digestor.setSpecificity(EnzymaticDigestion::SPEC_FULL);

  std::vector<AASequence> output;
  digestor.digest(protein, output);

  // generate decoy_factor number of complete decoy proteins
  std::vector<AASequence> decoy_proteins;
  for (int variant = 0; variant < decoy_factor; ++variant)
  {
    std::string decoy_sequence;
    for (const auto & aas : output)
    {
      if (aas.size() <= 2)
      {
        decoy_sequence += aas.toUnmodifiedString();
        continue;
      }

      // Important: create DecoyGenerator instance per peptide with same seed
      // Otherwise same peptides end up creating different decoys -> much more decoys than targets
      // But: we add variant to seed to get different decoys in multiple decoy generation
      DecoyGenerator dg;
      dg.setSeed(4711 + variant); // + variant to get different decoys in multiple decoy generation
      decoy_sequence += dg.shufflePeptides(aas, protease).toUnmodifiedString();
    }
    decoy_proteins.push_back(AASequence::fromString(decoy_sequence));
  }
  
  return decoy_proteins;
}

AASequence DecoyGenerator::shufflePeptides(
        const AASequence& protein,
        const std::string& protease,
        const int max_attempts)
{  
  OPENMS_PRECONDITION(!protein.isModified(), "Decoy generation only supports unmodified proteins.");

  std::vector<AASequence> peptides;
  ProteaseDigestion ed;
  ed.setMissedCleavages(0); // important as we want to reverse between all cutting sites
  ed.setEnzyme(protease);
  ed.setSpecificity(EnzymaticDigestion::SPEC_FULL);
  ed.digest(protein, peptides);    
  std::string protein_shuffled;
  for (int i = 0; i < static_cast<int>(peptides.size()) - 1; ++i)
  {
    const std::string peptide_string = peptides[i].toUnmodifiedString();

    // add from cache if available
    bool cached(false);
    #pragma omp critical (td_cache_)
    {
      auto it = td_cache_.find(peptide_string);
      if (it != td_cache_.end())
      {
        protein_shuffled += it->second; // add if cached
        cached = true;
      }
    }
    if (cached) continue;

    std::string peptide_string_shuffled = peptide_string;
    auto last = --peptide_string_shuffled.end();
    double lowest_identity(1.0);
    std::string lowest_identity_string(peptide_string_shuffled);
    for (int i = 0; i < max_attempts; ++i) // try to find sequence with low identity
    {
      shuffler_.portable_random_shuffle(std::begin(peptide_string_shuffled), last);

      double identity = SequenceIdentity_(peptide_string_shuffled, peptide_string);
      if (identity < lowest_identity)
      {
        lowest_identity = identity;
        lowest_identity_string = peptide_string_shuffled;

        if (identity <= (1.0/peptide_string_shuffled.size() + 1e-6)) 
        {
          break; // found perfect shuffle (only 1 (=cutting site) of all AAs match)
        }
      }
    }
    protein_shuffled += lowest_identity_string;
    #pragma omp critical (td_cache_)
    {
      td_cache_[peptide_string] = lowest_identity_string;
    }
  }
  // the last peptide of a protein is not an enzymatic cutting site so we do a full shuffle
  const std::string peptide_string = peptides[peptides.size() - 1 ].toUnmodifiedString();
  bool cached(false);
  #pragma omp critical (td_cache_)
  {
    auto it = td_cache_.find(peptide_string);
    if (it != td_cache_.end())
    {
      protein_shuffled += it->second; // add if cached
      cached = true;
    }
  }
  if (cached) return AASequence::fromString(protein_shuffled);

  std::string peptide_string_shuffled = peptide_string;
  double lowest_identity(1.0);
  std::string lowest_identity_string(peptide_string_shuffled);
  for (int i = 0; i < max_attempts; ++i) // try to find sequence with low identity
  {
    shuffler_.portable_random_shuffle(std::begin(peptide_string_shuffled), std::end(peptide_string_shuffled));
    double identity = SequenceIdentity_(peptide_string_shuffled, peptide_string);
    if (identity < lowest_identity)
    {
      lowest_identity = identity;
      lowest_identity_string = peptide_string_shuffled;
      if (identity == 0)
      {
        break; // found best shuffle
      }
    }
  }
  protein_shuffled += lowest_identity_string;
  #pragma omp critical (td_cache_)
  {
    td_cache_[peptide_string] = lowest_identity_string;
  }
  return AASequence::fromString(protein_shuffled);
}

// static
double DecoyGenerator::SequenceIdentity_(const std::string& decoy, const std::string& target)
{
  int match = 0;
  for (Size i = 0; i < target.size(); ++i)
  {
    if (target[i] == decoy[i]) { ++match; }
  }
  double identity = (double) match / target.size();

  // also compare against reverse
  match = 0;
  for (int i = (int)target.size() - 1; i >= 0; --i)
  {
    int j = (int)target.size() - 1 - i;
    if (target[j] == decoy[i]) { ++match; }
  }
  double rev_identity = (double) match / target.size();
   
  return std::max(identity, rev_identity);
}
