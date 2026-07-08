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
    class OPENMS_DLLAPI ONNXPredictorBase
    {
    public:
        /// @brief Constructor initializes the generic ONNX session safely across platforms.
        /// @param model_path Path to the ONNX model file.
        ONNXPredictorBase(const std::string& model_path);

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

    private:
        std::unique_ptr<Ort::SessionOptions> session_options_;
        std::unique_ptr<Ort::Session> session_;
        std::unique_ptr<Ort::MemoryInfo> memory_info_;
    };
} // namespace OpenMS
