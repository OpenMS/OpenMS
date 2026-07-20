// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <algorithm>
#include <string>
#include <unordered_map>

namespace
{
  constexpr float CHARGE_SCALE = 0.1f;
  constexpr float NCE_SCALE = 0.01f;

  size_t encodedLength_(const std::string& peptide, const OpenMS::ML::PeptDeepInputConfig& config)
  {
    // Note: AASequence::fromString is used to get the true length minus the modification brackets
    OpenMS::AASequence seq = OpenMS::AASequence::fromString(peptide);
    return seq.size() + (config.add_terminal_tokens ? 2 : 0);
  }

  // Generates the 109-element tensor based on atomic composition (Periodic Table index)
  std::vector<float> createModVector(int C, int H, int N, int O, int P, int S) {
      std::vector<float> mod(109, 0.0f);
      if (H != 0) mod[0]  = static_cast<float>(H); // Hydrogen (Atomic #1)
      if (C != 0) mod[5]  = static_cast<float>(C); // Carbon   (Atomic #6)
      if (N != 0) mod[6]  = static_cast<float>(N); // Nitrogen (Atomic #7)
      if (O != 0) mod[7]  = static_cast<float>(O); // Oxygen   (Atomic #8)
      if (P != 0) mod[14] = static_cast<float>(P); // Phospho  (Atomic #15)
      if (S != 0) mod[15] = static_cast<float>(S); // Sulfur   (Atomic #16)
      return mod;
  }

  // Map OpenMS string names to the AlphaPeptDeep elemental tensors
  const std::unordered_map<std::string, std::vector<float>> PEPTDEEP_MOD_DICT = {
      {"Oxidation", createModVector(0, 0, 0, 1, 0, 0)},          // +1 Oxygen
      {"Phospho", createModVector(0, 1, 0, 3, 1, 0)},            // +1 H, +3 O, +1 P
      {"Acetyl", createModVector(2, 2, 0, 1, 0, 0)},             // +2 C, +2 H, +1 O
      {"Carbamidomethyl", createModVector(2, 3, 1, 1, 0, 0)},    // +2 C, +3 H, +1 N, +1 O
      {"Deamidation", createModVector(0, -1, -1, 1, 0, 0)},      // +1 O, -1 H, -1 N
      {"Amidation", createModVector(0, 1, 1, -1, 0, 0)}          // +1 H, +1 N, -1 O
  };
}

namespace OpenMS
{
  namespace ML
  {
    PeptDeepInputBatch PeptDeepInputBuilder::buildPeptideBatch(
      const std::vector<std::string>& peptides,
      const PeptDeepInputConfig& config)
    {
      if (peptides.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
      }

      size_t longest_encoded = 0;
      for (const std::string& peptide : peptides)
      {
        validatePeptide(peptide);
        longest_encoded = std::max(longest_encoded, encodedLength_(peptide, config));
      }

      const size_t sequence_length = config.fixed_sequence_length == 0 ? longest_encoded : config.fixed_sequence_length;

      PeptDeepInputBatch batch;
      batch.batch_size = peptides.size();
      batch.sequence_length = sequence_length;
      batch.aa_indices.reserve(batch.batch_size * batch.sequence_length);

      // Zero-fill the mod_x tensor upfront for the whole batch
      batch.mod_x.assign(batch.batch_size * batch.sequence_length * PEPTDEEP_MOD_ELEMENTS, 0.0f);

      for (size_t batch_idx = 0; batch_idx < peptides.size(); ++batch_idx)
      {
        OpenMS::AASequence seq = OpenMS::AASequence::fromString(peptides[batch_idx]);
        size_t written = 0;

        if (config.add_terminal_tokens)
        {
          batch.aa_indices.push_back(0);
          ++written;
        }

        for (size_t i = 0; i < seq.size(); ++i)
        {
          // 1. Tokenize the un-modified base amino acid
          char aa = seq[i].getOneLetterCode()[0];
          batch.aa_indices.push_back(getAAIndex(aa));

          // 2. Extract Modifications
          if (seq[i].isModified())
          {
            std::string mod_name = seq[i].getModificationName();
            size_t tensor_offset = (batch_idx * sequence_length * PEPTDEEP_MOD_ELEMENTS) + (written * PEPTDEEP_MOD_ELEMENTS);

            //  Mapping: Look up the modification in our dictionary
            auto it = PEPTDEEP_MOD_DICT.find(mod_name);
            if (it != PEPTDEEP_MOD_DICT.end())
            {
                // Overwrite the zeros with the 109-element chemical features
                const std::vector<float>& mod_features = it->second;
                for (size_t elem = 0; elem < PEPTDEEP_MOD_ELEMENTS; ++elem)
                {
                    batch.mod_x[tensor_offset + elem] = mod_features[elem];
                }
            }
            else
            {
                // Throw a clean exception if an unsupported modification is passed in
                throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unsupported PeptDeep modification: " + mod_name);
            }
          }

          ++written;
        }

        if (config.add_terminal_tokens)
        {
          batch.aa_indices.push_back(0);
          ++written;
        }

        while (written < batch.sequence_length)
        {
          batch.aa_indices.push_back(0);
          ++written;
        }
      }

      return batch;
    }

    PeptDeepInputBatch PeptDeepInputBuilder::buildChargedBatch(
      const std::vector<std::string>& peptides,
      const std::vector<float>& charges,
      const PeptDeepInputConfig& config)
    {
      const size_t batch_size = peptides.size();
      if (charges.size() != batch_size)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide and charge input vectors must have the same size.");
      }

      PeptDeepInputBatch batch = buildPeptideBatch(peptides, config);
      batch.charges.reserve(batch_size);

      for (size_t i = 0; i < batch_size; ++i)
      {
        batch.charges.push_back(charges[i] * CHARGE_SCALE);
      }

      return batch;
    }

    PeptDeepInputBatch PeptDeepInputBuilder::buildInstrumentBatch(
      const std::vector<std::string>& peptides,
      const std::vector<float>& charges,
      const std::vector<float>& nces,
      const std::vector<int64_t>& instrument_indices,
      const PeptDeepInputConfig& config)
    {
      const size_t batch_size = peptides.size();
      if (charges.size() != batch_size || nces.size() != batch_size || instrument_indices.size() != batch_size)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide, charge, NCE, and instrument input vectors must have the same size.");
      }

      PeptDeepInputBatch batch = buildChargedBatch(peptides, charges, config);
      batch.nces.reserve(batch_size);
      batch.instrument_indices.reserve(batch_size);

      for (size_t i = 0; i < batch_size; ++i)
      {
        batch.nces.push_back(nces[i] * NCE_SCALE);
        batch.instrument_indices.push_back(instrument_indices[i]);
      }

      return batch;
    }
  }
}