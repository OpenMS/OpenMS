// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, kg290 $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <chrono>
#include <algorithm>
#include <map>
#include <random>

using namespace OpenMS;

DecoyGenerator::DecoyGenerator()
{
  const UInt64 seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  shuffler_.seed(seed);
}

void DecoyGenerator::setSeed(UInt64 seed)
{
  shuffler_.seed(seed);
}

void DecoyGenerator::startDeBruijn(Size k, UInt64 seed, const std::string& keep_residues)
{
  OPENMS_PRECONDITION(k > 0, "The de Bruijn k-mer size must be positive.")
  debruijn_k_ = k;
  debruijn_seed_ = seed;
  debruijn_keep_residues_ = keep_residues;
  debruijn_edge_counts_.clear();
  debruijn_residue_counts_.clear();
  debruijn_edge_labels_.clear();
  debruijn_ready_ = false;
}

void DecoyGenerator::addProteinToDeBruijn(const AASequence& protein)
{
  OPENMS_PRECONDITION(debruijn_k_ > 0 && !debruijn_ready_, "Call startDeBruijn() before adding target proteins and finalize only after all proteins are added.")
  OPENMS_PRECONDITION(!protein.isModified(), "De Bruijn decoy generation only supports unmodified proteins.")

  const std::string sequence = protein.toUnmodifiedString();
  const std::string padded = std::string(debruijn_k_, '-') + sequence;
  for (Size i = 0; i < sequence.size(); ++i)
  {
    ++debruijn_edge_counts_[padded.substr(i, debruijn_k_ + 1)];
    ++debruijn_residue_counts_[sequence[i]];
  }
}

void DecoyGenerator::finalizeDeBruijn()
{
  OPENMS_PRECONDITION(debruijn_k_ > 0 && !debruijn_ready_, "Call startDeBruijn() once before finalizing de Bruijn decoys.")

  // Stable edge order keeps seeded decoys reproducible across standard-library hash implementations.
  std::vector<std::pair<std::string, Size>> edges(debruijn_edge_counts_.begin(), debruijn_edge_counts_.end());
  std::sort(edges.begin(), edges.end());

  std::map<char, Size> mutable_residue_counts;
  for (const auto& [residue, count] : debruijn_residue_counts_)
  {
    mutable_residue_counts[residue] = count;
  }
  for (const auto& [edge, count] : edges)
  {
    if (debruijn_keep_residues_.find(edge.back()) != std::string::npos)
    {
      mutable_residue_counts[edge.back()] -= count;
      debruijn_edge_labels_[edge] = edge.back();
    }
  }

  std::vector<char> alphabet;
  std::vector<double> weights;
  for (const auto& [residue, count] : mutable_residue_counts)
  {
    if (count > 0)
    {
      alphabet.push_back(residue);
      weights.push_back(static_cast<double>(count));
    }
  }
  std::mt19937_64 rng(debruijn_seed_);
  if (alphabet.empty())
  {
    OPENMS_PRECONDITION(debruijn_edge_labels_.size() == edges.size(), "No replacement residues are available for the de Bruijn decoy.")
  }
  else
  {
    std::discrete_distribution<Size> residue_distribution(weights.begin(), weights.end());
    for (const auto& [edge, count] : edges)
    {
      if (debruijn_edge_labels_.contains(edge)) continue;
      debruijn_edge_labels_[edge] = alphabet[residue_distribution(rng)];
    }
  }
  debruijn_ready_ = true;
}

AASequence DecoyGenerator::deBruijn(const AASequence& protein) const
{
  OPENMS_PRECONDITION(debruijn_ready_, "Call finalizeDeBruijn() before generating a de Bruijn decoy.")
  OPENMS_PRECONDITION(!protein.isModified(), "De Bruijn decoy generation only supports unmodified proteins.")

  const std::string sequence = protein.toUnmodifiedString();
  const std::string padded = std::string(debruijn_k_, '-') + sequence;
  std::string decoy;
  decoy.reserve(sequence.size());
  for (Size i = 0; i < sequence.size(); ++i)
  {
    const auto it = debruijn_edge_labels_.find(padded.substr(i, debruijn_k_ + 1));
    OPENMS_PRECONDITION(it != debruijn_edge_labels_.end(), "The protein was not included when preparing the de Bruijn decoy database.")
    decoy.push_back(it->second);
  }
  return AASequence::fromString(decoy);
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


