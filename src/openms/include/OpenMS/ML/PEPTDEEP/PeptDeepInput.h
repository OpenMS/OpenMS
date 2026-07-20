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
    /// Dynamically handles both unmodified and modified peptides.
    class OPENMS_DLLAPI PeptDeepInputBuilder
    {
    public:
      static PeptDeepInputBatch buildPeptideBatch(
        const std::vector<std::string>& peptides,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      static PeptDeepInputBatch buildChargedBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      static PeptDeepInputBatch buildInstrumentBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());
    };
  } // namespace ML
} // namespace OpenMS
