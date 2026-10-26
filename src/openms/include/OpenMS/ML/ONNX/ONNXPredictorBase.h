// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <memory>
#include <vector>

// Forward declaration to avoid exposing ONNX headers globally
namespace Ort {
    struct Session;
    struct SessionOptions;
    struct MemoryInfo;
}

namespace OpenMS
{
    /**
      @brief Owns the ONNX Runtime session of one model file, for the predictors built on it

      Loads the model with all graph optimizations, without ONNX Runtime's memory pattern (the inputs of
      the models change shape from call to call), and gives a predictor (PeptDeepRTInference,
      PeptDeepCCSInference, PeptDeepMS2Inference) the session and the CPU memory info it needs to build
      input tensors. The ONNX Runtime types are only forward-declared here. Available in builds with
      @c WITH_ONNX.

      A predictor splits its input into model runs (calls of Ort::Session::Run(), one per batch of
      peptides of one length) and runs them with forEachBatch(): up to the constructor's @p threads runs
      at once, each on one thread (the session has one intra-op thread). On CPU, the small models of PeptDeep scale much
      better over their calls than within one call (2.2x faster with 16 threads), and the result of a
      call does not depend on the number of its threads.
    */
    class OPENMS_DLLAPI ONNXPredictorBase
    {
    public:
        /// @brief Constructor initializes the generic ONNX session safely across platforms.
        /// @param model_path Path to the ONNX model file.
        /// @param threads The model runs forEachBatch() runs at once (each on one thread); at least 1.
        explicit ONNXPredictorBase(const std::string& model_path, int threads = 4);

        virtual ~ONNXPredictorBase();

        ONNXPredictorBase(const ONNXPredictorBase&) = delete;
        ONNXPredictorBase& operator=(const ONNXPredictorBase&) = delete;

        ONNXPredictorBase(ONNXPredictorBase&&) noexcept;
        ONNXPredictorBase& operator=(ONNXPredictorBase&&) noexcept;

        /// @brief Access the wrapped ONNX Runtime session.
        Ort::Session& session();
        const Ort::Session& session() const;

        /// @brief Access CPU memory info used to construct input tensors.
        Ort::MemoryInfo& memoryInfo();
        const Ort::MemoryInfo& memoryInfo() const;

        /// @brief Return model input names in ONNX graph order.
        std::vector<std::string> getInputNames() const;

        /// @brief Return model output names in ONNX graph order.
        std::vector<std::string> getOutputNames() const;

        /// @brief Return the declared input shape for an input index.
        std::vector<int64_t> getInputShape(size_t input_index) const;

        /// @brief Return the number of model inputs.
        size_t getInputCount() const;

        /// @brief Return the number of model outputs.
        size_t getOutputCount() const;

        /**
          @brief Runs @p run on each batch of @p groups, up to the constructor's @p threads batches at once (with OpenMP)

          The batches of a group are its consecutive slices of at most @p batch_size indices (the last one shorter);
          @p run(indices, n) gets one of them, e.g. to run the session on those peptides (Ort::Session::Run() may be
          called concurrently) and store their predictions. Once all have run, the first exception of a batch, if any,
          is rethrown.
        */
        void forEachBatch(const std::vector<std::vector<size_t>>& groups, size_t batch_size,
                          const std::function<void(const size_t* indices, size_t n)>& run) const;

    private:
        int threads_ = 1;
        std::unique_ptr<Ort::SessionOptions> session_options_;
        std::unique_ptr<Ort::Session> session_;
        std::unique_ptr<Ort::MemoryInfo> memory_info_;
    };
} // namespace OpenMS
