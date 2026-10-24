// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <OpenMS/config.h>
#include <OpenMS/ML/ONNX/ONNXPredictorBase.h>

namespace OpenMS {

/// @brief Inference engine for PeptDeep MS2 Fragment Intensity prediction.
///
/// This class parses and tokenizes raw peptide string sequences into structural integer arrays
/// and evaluates them through the ONNX Runtime backend to predict fragment intensities.
class OPENMS_DLLAPI PeptDeepMS2Inference {
public:
    /// @brief Constructor initializes the ONNX session with the specified model graph.
    /// @param model_path Absolute or relative path to the PeptDeep MS2 ONNX model file.
    /// @param threads The model runs (of at most @p batch_size peptides of one length each) that run at once, each on one thread (default 4).
    /// @param batch_size Maximum number of peptides to process in a single ONNX run (default 500).
    /// @throws Exception::BaseException if the model file cannot be found or loaded by the runtime.
    explicit PeptDeepMS2Inference(const std::string& model_path, int threads = 4, size_t batch_size = 500);

    /// @brief Destructor
    ~PeptDeepMS2Inference();

    /// @brief Predicts MS2 fragment intensities for a batch of peptide sequences.
    /// @param peptides A vector of peptide strings. Supports OpenMS AASequence modification notation (e.g., "PEPTIDEK", "M(Oxidation)PEP").
    /// @param charges A vector of precursor charge states (must match peptides size).
    /// @param nces A vector of normalized collision energies (must match peptides size).
    /// @param instrument_indices One peptdeep instrument index per peptide. Derive these with
    ///        ML::instrumentIndex() rather than writing literals: the encoding is positional and
    ///        belongs to peptdeep (ML::ALPHAPEPTDEEP_INSTRUMENTS), where 0 is QE and 1 is Lumos.
    ///        Indices outside [0, ML::PEPTDEEP_MAX_INSTRUMENT_NUM) are rejected, because the
    ///        model's OneHot would otherwise answer them with an all-off row.
    /// @return A vector of flattened fragment intensity arrays, one for each peptide.
    ///         Each cleavage position contains eight contiguous AlphaPeptDeep channels in this order:
    ///         [b_z1, b_z2, y_z1, y_z2, b_modloss_z1, b_modloss_z2, y_modloss_z1, y_modloss_z2].
    /// @throws Exception::IllegalArgument if peptides is empty, size constraints fail, or a sequence is chemically invalid.
    std::vector<std::vector<float>> predictMS2(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices);

    /// The values per cleavage position of a prediction of predictMS2(), in this order (Channel)
    static constexpr size_t ION_CHANNELS = 8;
    /// The channels of a cleavage position: [b_z1, b_z2, y_z1, y_z2, b_modloss_z1, b_modloss_z2, y_modloss_z1, y_modloss_z2]
    enum Channel : size_t { B_Z1 = 0, B_Z2 = 1, Y_Z1 = 2, Y_Z2 = 3 };
    /// ionIndex() of an ion the prediction does not have
    static constexpr size_t NO_ION = static_cast<size_t>(-1);

    /**
      @brief The index into a prediction of predictMS2() of the b or y ion of a peptide

      Cleavage position @c pos (0-based) carries b(pos + 1) and y(length - 1 - pos), with @p channels values each.

      @param length The peptide's length (residues)
      @param type 'b' or 'y'
      @param ordinal The ion's ordinal, 1 to @p length - 1
      @param charge 1 or 2
      @param channels The values per cleavage position (ION_CHANNELS, unless a model has another layout)
      @return The index, NO_ION for another ion type, ordinal or charge
    */
    static size_t ionIndex(size_t length, char type, size_t ordinal, int charge, size_t channels = ION_CHANNELS);

private:
    ONNXPredictorBase model_;
    size_t batch_size_;
};

} // namespace OpenMS