// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ML/ONNX/ONNXPredictorBase.h>
#include <string>
#include <vector>

namespace OpenMS
{
    class OPENMS_DLLAPI PeptDeepRTInference
    {
    public:
        /**
         * @brief Constructor initializes the ONNX environment and loads the model
         * @param model_path Absolute path to peptdeep_rt_dynamic.onnx
         * @param intra_op_threads Number of ONNX execution threads (default 4).
         * @param batch_size Maximum number of peptides to process in a single ONNX run (default 500).
         */
        explicit PeptDeepRTInference(const std::string& model_path, int intra_op_threads = 4, size_t batch_size = 500);

        /**
         * @brief Destructor
         */
        ~PeptDeepRTInference();

        /**
         * @brief Predicts Retention Times for a list of peptide sequences
         * @param peptides A vector of raw uppercase, unmodified peptide strings (e.g., "PEPTIDEK")
         * @return A vector of predicted RT values corresponding to the input peptides
         * @throws Exception::IllegalArgument if peptides is empty or contains an invalid/modified sequence.
         */
        std::vector<float> predictRT(const std::vector<std::string>& peptides);

    private:
        ONNXPredictorBase model_;
        size_t batch_size_;
    };
} // namespace OpenMS
