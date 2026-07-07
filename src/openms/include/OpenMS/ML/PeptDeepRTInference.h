// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ML/ONNXPredictorBase.h>
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
         */
        explicit PeptDeepRTInference(const std::string& model_path);

        /**
         * @brief Destructor
         */
        ~PeptDeepRTInference();

        /**
         * @brief Predicts Retention Times for a list of peptide sequences
         * @param peptides A vector of peptide strings (e.g., "PEPTIDEK")
         * @return A vector of predicted RT values corresponding to the input peptides
         */
        std::vector<float> predictRT(const std::vector<std::string>& peptides);

    private:
        ONNXPredictorBase model_;
    };
} // namespace OpenMS
