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
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <algorithm>
#include <string>
#include <vector>

namespace
{
  constexpr float CHARGE_SCALE = 0.1f;
  constexpr float NCE_SCALE = 0.01f;

  size_t encodedLength_(const OpenMS::AASequence& seq, const OpenMS::ML::PeptDeepInputConfig& config)
  {
    return seq.size() + (config.add_terminal_tokens ? 2 : 0);
  }
}

namespace OpenMS
{
  namespace ML
  {
    PeptDeepInputBatch PeptDeepInputBuilder::buildPeptideBatch(
      const std::vector<OpenMS::AASequence>& peptides,
      const PeptDeepInputConfig& config)
    {
      if (peptides.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
      }

      size_t longest_encoded = 0;
      for (const auto& seq : peptides)
      {
        longest_encoded = std::max(longest_encoded, encodedLength_(seq, config));
      }

      if (config.fixed_sequence_length > 0 && longest_encoded > config.fixed_sequence_length)
      {
         throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide length exceeds fixed sequence length.");
      }

      const size_t sequence_length = config.fixed_sequence_length == 0 ? longest_encoded : config.fixed_sequence_length;

      PeptDeepInputBatch batch;
      batch.batch_size = peptides.size();
      batch.sequence_length = sequence_length;
      batch.aa_indices.reserve(batch.batch_size * batch.sequence_length);

      batch.mod_x.assign(batch.batch_size * batch.sequence_length * PEPTDEEP_MOD_ELEMENTS, 0.0f);

      for (size_t batch_idx = 0; batch_idx < peptides.size(); ++batch_idx)
      {
        const OpenMS::AASequence& seq = peptides[batch_idx];
        size_t written = 0;

        // Apply modifications by placing them at their specific AlphaPeptDeep index
        auto apply_modification = [&](const OpenMS::ResidueModification* mod, size_t site_index) {
          if (mod != nullptr)
          {
              OpenMS::EmpiricalFormula diff_formula = mod->getDiffFormula();
              size_t tensor_offset = (batch_idx * sequence_length * PEPTDEEP_MOD_ELEMENTS) + (site_index * PEPTDEEP_MOD_ELEMENTS);

              for (auto const& element_count : diff_formula)
              {
                  std::string symbol = element_count.first->getSymbol();
                  int count = element_count.second;

                  // Normalize OpenMS isotope format "(13)C" to AlphaPeptDeep format "13C"
                  if (symbol.front() == '(') {
                      size_t close_paren = symbol.find(')');
                      if (close_paren != std::string::npos) {
                          symbol = symbol.substr(1, close_paren - 1) + symbol.substr(close_paren + 1);
                      }
                  }

                  auto it = std::find(ALPHAPEPTDEEP_MOD_ELEMENTS.begin(), ALPHAPEPTDEEP_MOD_ELEMENTS.end(), symbol);
                  int tensor_elem_index = static_cast<int>(PEPTDEEP_MOD_ELEMENTS) - 1; // Default to "Other"

                  if (it != ALPHAPEPTDEEP_MOD_ELEMENTS.end())
                  {
                      tensor_elem_index = std::distance(ALPHAPEPTDEEP_MOD_ELEMENTS.begin(), it);
                      if (tensor_elem_index >= static_cast<int>(PEPTDEEP_MOD_ELEMENTS))
                      {
                          tensor_elem_index = static_cast<int>(PEPTDEEP_MOD_ELEMENTS) - 1;
                      }
                  }

                  batch.mod_x[tensor_offset + tensor_elem_index] += static_cast<float>(count);
              }
          }
        };

        if (config.add_terminal_tokens)
        {
          batch.aa_indices.push_back(0);
          if (seq.hasNTerminalModification())
          {
              apply_modification(seq.getNTerminalModification(), written);
          }
          ++written;
        }
        else if (seq.hasNTerminalModification())
        {
          apply_modification(seq.getNTerminalModification(), written);
        }

        for (size_t i = 0; i < seq.size(); ++i)
        {
          char aa = seq[i].getOneLetterCode()[0];
          batch.aa_indices.push_back(getAAIndex(aa));

          if (seq[i].isModified())
          {
              apply_modification(seq[i].getModification(), written);
          }
          ++written;
        }

        if (config.add_terminal_tokens)
        {
          batch.aa_indices.push_back(0);
          if (seq.hasCTerminalModification())
          {
              apply_modification(seq.getCTerminalModification(), written);
          }
          ++written;
        }
        else if (seq.hasCTerminalModification())
        {
           // Explicitly guard against empty sequence underflow before writing offset
           if (written > 0)
           {
               apply_modification(seq.getCTerminalModification(), written - 1);
           }
        }

        while (written < batch.sequence_length)
        {
          batch.aa_indices.push_back(0);
          ++written;
        }
      }

      return batch;
    }

    PeptDeepInputBatch PeptDeepInputBuilder::buildPrecursorBatch(
      const std::vector<OpenMS::AASequence>& peptides,
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

    PeptDeepInputBatch PeptDeepInputBuilder::buildProductMetaBatch(
      const std::vector<OpenMS::AASequence>& peptides,
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

      PeptDeepInputBatch batch = buildPrecursorBatch(peptides, charges, config);
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