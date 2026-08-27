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
  class AASequence;

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
    /// The buildUnmodified*() overloads take plain sequence strings and emit a
    /// zero-filled mod_x, i.e. they describe unmodified peptides. The buildModified*()
    /// overloads take AASequence and encode each modification's elemental composition
    /// into mod_x via PeptDeepModEncoder, which is what the models expect for modified
    /// peptides.
    class OPENMS_DLLAPI PeptDeepInputBuilder
    {
    public:
      /// @brief Tokenizes a batch of unmodified peptides into a padded aa_indices/mod_x tensor pair.
      /// @param peptides A vector of raw uppercase peptide strings. Must be non-empty.
      /// @param config Padding policy (terminal tokens, fixed sequence length).
      /// @return A PeptDeepInputBatch with aa_indices and a zero-filled mod_x populated; charges/nces/instrument_indices left empty.
      /// @throws Exception::IllegalArgument if peptides is empty or contains an invalid sequence (see validatePeptide()).
      /// @throws Exception::InvalidValue if config.fixed_sequence_length is shorter than the longest encoded peptide.
      static PeptDeepInputBatch buildUnmodifiedPeptideBatch(
        const std::vector<std::string>& peptides,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Like buildUnmodifiedPeptideBatch(), additionally populating the (scaled) charges field.
      /// @param peptides A vector of raw uppercase peptide strings.
      /// @param charges Precursor charge states, one per peptide. Must match peptides size.
      /// @param config Padding policy (terminal tokens, fixed sequence length).
      /// @return A PeptDeepInputBatch with aa_indices, mod_x, and charges populated.
      /// @throws Exception::IllegalArgument if charges.size() != peptides.size(), or via buildUnmodifiedPeptideBatch().
      static PeptDeepInputBatch buildUnmodifiedChargedBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Like buildUnmodifiedChargedBatch(), additionally populating the (scaled) nces and instrument_indices fields.
      /// @param peptides A vector of raw uppercase peptide strings.
      /// @param charges Precursor charge states, one per peptide. Must match peptides size.
      /// @param nces Normalized collision energies, one per peptide. Must match peptides size.
      /// @param instrument_indices Categorical instrument identifiers, one per peptide. Must match peptides size.
      /// @param config Padding policy (terminal tokens, fixed sequence length).
      /// @return A PeptDeepInputBatch with all fields populated.
      /// @throws Exception::IllegalArgument if any of charges/nces/instrument_indices does not match peptides.size(), or via buildUnmodifiedChargedBatch().
      static PeptDeepInputBatch buildUnmodifiedInstrumentBatch(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Tokenizes a batch of possibly modified peptides, encoding modifications into mod_x.
      /// @param peptides Peptides whose modifications are taken from the AASequence itself.
      /// @param config Padding policy (terminal tokens, fixed sequence length).
      /// @return A batch with aa_indices and a modification-aware mod_x; charges/nces/instrument_indices left empty.
      /// @throws Exception::IllegalArgument if peptides is empty or an unmodified sequence is invalid.
      /// @throws Exception::InvalidValue if config.fixed_sequence_length is shorter than the longest encoded peptide.
      static PeptDeepInputBatch buildModifiedPeptideBatch(
        const std::vector<AASequence>& peptides,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Like buildModifiedPeptideBatch(), additionally populating the (scaled) charges field.
      /// @throws Exception::IllegalArgument if charges.size() != peptides.size().
      static PeptDeepInputBatch buildModifiedChargedBatch(
        const std::vector<AASequence>& peptides,
        const std::vector<float>& charges,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());

      /// @brief Like buildModifiedChargedBatch(), additionally populating nces and instrument_indices.
      /// @throws Exception::IllegalArgument if any input vector size differs from peptides.size().
      static PeptDeepInputBatch buildModifiedInstrumentBatch(
        const std::vector<AASequence>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices,
        const PeptDeepInputConfig& config = PeptDeepInputConfig());
    };
  } // namespace ML
} // namespace OpenMS
