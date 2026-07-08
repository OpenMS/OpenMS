// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PeptDeepInput.h>
#include <OpenMS/ML/PeptDeepUtils.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <string>

namespace
{
  constexpr float CHARGE_SCALE = 0.1f;
  constexpr float NCE_SCALE = 0.01f;

  size_t encodedLength_(const std::string& peptide, const OpenMS::ML::PeptDeepInputConfig& config)
  {
    return peptide.size() + (config.add_terminal_tokens ? 2 : 0);
  }
}

namespace OpenMS
{
  namespace ML
  {
    PeptDeepInputBatch PeptDeepInputBuilder::buildUnmodifiedPeptideBatch(
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
      if (sequence_length > PEPTDEEP_MAX_SEQUENCE_LENGTH)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Encoded peptide sequence exceeds the maximum allowed PeptDeep input length.",
          std::to_string(sequence_length));
      }
      if (longest_encoded > sequence_length)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Configured PeptDeep sequence length is shorter than an encoded peptide.",
          std::to_string(sequence_length));
      }

      PeptDeepInputBatch batch;
      batch.batch_size = peptides.size();
      batch.sequence_length = sequence_length;
      batch.aa_indices.reserve(batch.batch_size * batch.sequence_length);

      for (const std::string& peptide : peptides)
      {
        size_t written = 0;

        if (config.add_terminal_tokens)
        {
          batch.aa_indices.push_back(0);
          ++written;
        }

        for (const char aa : peptide)
        {
          batch.aa_indices.push_back(getAAIndex(aa));
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

      batch.mod_x.assign(batch.batch_size * batch.sequence_length * PEPTDEEP_MOD_ELEMENTS, 0.0f);
      return batch;
    }

    PeptDeepInputBatch PeptDeepInputBuilder::buildUnmodifiedChargedBatch(
      const std::vector<std::string>& peptides,
      const std::vector<float>& charges,
      const PeptDeepInputConfig& config)
    {
      const size_t batch_size = peptides.size();
      if (charges.size() != batch_size)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Peptide and charge input vectors must have the same size.");
      }

      PeptDeepInputBatch batch = buildUnmodifiedPeptideBatch(peptides, config);
      batch.charges.reserve(batch_size);

      for (size_t i = 0; i < batch_size; ++i)
      {
        batch.charges.push_back(charges[i] * CHARGE_SCALE);
      }

      return batch;
    }

    PeptDeepInputBatch PeptDeepInputBuilder::buildUnmodifiedInstrumentBatch(
      const std::vector<std::string>& peptides,
      const std::vector<float>& charges,
      const std::vector<float>& nces,
      const std::vector<int64_t>& instrument_indices,
      const PeptDeepInputConfig& config)
    {
      const size_t batch_size = peptides.size();
      if (charges.size() != batch_size || nces.size() != batch_size || instrument_indices.size() != batch_size)
      {
        throw Exception::IllegalArgument(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Peptide, charge, NCE, and instrument input vectors must have the same size.");
      }

      PeptDeepInputBatch batch = buildUnmodifiedChargedBatch(peptides, charges, config);
      batch.nces.reserve(batch_size);
      batch.instrument_indices.reserve(batch_size);

      for (size_t i = 0; i < batch_size; ++i)
      {
        batch.nces.push_back(nces[i] * NCE_SCALE);
        batch.instrument_indices.push_back(instrument_indices[i]);
      }

      return batch;
    }
  } // namespace ML
} // namespace OpenMS
