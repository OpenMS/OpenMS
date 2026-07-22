// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS
{
  namespace ML
  {
    /// @brief Sequence-padding policy for PeptDeep ONNX inputs.
    struct OPENMS_DLLAPI PeptDeepInputConfig
    {
      /// Add PeptDeep terminal padding tokens around each peptide.
      bool add_terminal_tokens = true;

      /// If non-zero, pad every batch to this fixed sequence length.
      /// If zero, pad only to the longest encoded peptide in the batch.
      size_t fixed_sequence_length = 0;
    };

    /// @brief Flat tensor buffers and dimensions for a PeptDeep ONNX batch.
    struct OPENMS_DLLAPI PeptDeepInputBatch
    {
      size_t batch_size = 0;
      size_t sequence_length = 0;

      std::vector<int64_t> aa_indices;
      std::vector<float> mod_x;

      std::vector<float> charges;
      std::vector<float> nces;
      std::vector<int64_t> instrument_indices;
    };

    /// @brief Shared featurization for PeptDeep RT, CCS, and MS2 ONNX predictors.
    ///
    /// Dynamically handles both unmodified and modified peptides, parsing modifications
    /// natively into a 109-element elemental feature tensor using empirical formulas.
    class OPENMS_DLLAPI PeptDeepInputBuilder
    {
    public:
      /// @brief Builds a baseline tensor batch containing amino acid indices and modification features.
      /// @param peptides A vector of peptide sequences (supports OpenMS AASequence bracket notation).
      /// @param config Configuration dictating padding and terminal token policies.
      /// @return A PeptDeepInputBatch populated with sequence lengths, aa_indices, and mod_x.
      /// @throws Exception::IllegalArgument If the batch is empty or a peptide exceeds a fixed sequence length.
      static PeptDeepInputBatch buildPeptideBatch(
        const std::vector<std::string>& peptides,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Builds a tensor batch incorporating precursor charge states (primarily for CCS and RT).
      /// @param peptides A vector of peptide sequences.
      /// @param charges A vector of precursor charges corresponding to the peptides.
      /// @param config Configuration dictating padding and terminal token policies.
      /// @return A PeptDeepInputBatch populated with baseline features and scaled charges.
      /// @throws Exception::IllegalArgument If the size of the peptides and charges vectors do not match.
      static PeptDeepInputBatch buildChargedBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Builds a full tensor batch including Normalized Collision Energies and instrument details (for MS2).
      /// @param peptides A vector of peptide sequences.
      /// @param charges A vector of precursor charges corresponding to the peptides.
      /// @param nces A vector of Normalized Collision Energies (NCE).
      /// @param instrument_indices A vector of instrument identifier indices.
      /// @param config Configuration dictating padding and terminal token policies.
      /// @return A fully populated PeptDeepInputBatch ready for MS2 ONNX inference.
      /// @throws Exception::IllegalArgument If any of the input vectors differ in size.
      static PeptDeepInputBatch buildInstrumentBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());
    };
  } // namespace ML
} // namespace OpenMS