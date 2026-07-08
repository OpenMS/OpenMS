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
    /// @throws Exception::BaseException if the model file cannot be found or loaded by the runtime.
    PeptDeepMS2Inference(const std::string& model_path);

    /// @brief Destructor
    ~PeptDeepMS2Inference();

    /// @brief Predicts MS2 fragment intensities for a batch of peptide sequences.
    /// @param peptides A vector of raw uppercase peptide strings.
    /// @param charges A vector of precursor charge states (must match peptides size).
    /// @param nces A vector of normalized collision energies (must match peptides size).
    /// @param instrument_indices A vector of categorical integers representing MS instruments (e.g., 0=Lumos, 1=QE, 2=timsTOF, 3=Sciex).
    /// @return A vector of flattened fragment intensity arrays, one for each peptide.
    ///         Native shape ordering is contiguous by fragment and ion type [b_1, y_1, b_2, y_2...].
    /// @throws std::invalid_argument if the sequence string is empty, too long, or contains unmapped characters.
    std::vector<std::vector<float>> predictMS2(
        const std::vector<std::string>& peptides,
        const std::vector<float>& charges,
        const std::vector<float>& nces,
        const std::vector<int64_t>& instrument_indices);

private:
    ONNXPredictorBase model_;
};

} // namespace OpenMS
