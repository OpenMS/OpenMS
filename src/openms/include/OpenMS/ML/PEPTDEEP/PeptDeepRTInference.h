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
    /**
      @brief Predicts peptide retention times with the PeptDeep RT model (@c peptdeep_rt_dynamic.onnx)

      Runs the model through ONNXPredictorBase and returns one predicted retention time per peptide
      sequence, processing at most @c batch_size peptides per model run. Available in builds with
      @c WITH_ONNX.
    */
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
         * @brief Predicts Retention Times for a list of peptide sequences.
         * @param peptides A vector of peptide strings. Supports OpenMS AASequence modification notation (e.g., "PEPTIDEK", "M(Oxidation)PEP").
         * @return A vector of predicted RT values corresponding to the input peptides.
         * @throws Exception::IllegalArgument if peptides is empty, size constraints fail, or a sequence is chemically invalid.
         */
        std::vector<float> predictRT(const std::vector<std::string>& peptides);

    private:
        ONNXPredictorBase model_;
        size_t batch_size_;
    };
} // namespace OpenMS