// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ML/ONNXPredictorBase.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /// @brief Inference engine for PeptDeep collisional cross section prediction.
  class OPENMS_DLLAPI PeptDeepCCSInference
  {
  public:
    /// @brief Constructor initializes the ONNX session with the specified model graph.
    /// @param model_path Absolute or relative path to the PeptDeep CCS ONNX model file.
    explicit PeptDeepCCSInference(const std::string& model_path);

    /// @brief Destructor
    ~PeptDeepCCSInference();

    /// @brief Predicts CCS values for a batch of unmodified peptide sequences.
    /// @param peptides A vector of raw uppercase peptide strings.
    /// @param charges A vector of precursor charge states. Must match peptides size.
    /// @return A vector of predicted CCS values corresponding to the input peptides.
    std::vector<float> predictCCS(
      const std::vector<std::string>& peptides,
      const std::vector<float>& charges);

  private:
    ONNXPredictorBase model_;
  };
} // namespace OpenMS
