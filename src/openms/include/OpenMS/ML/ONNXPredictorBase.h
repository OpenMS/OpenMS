// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <string>
#include <memory>

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

    protected:
        std::unique_ptr<Ort::SessionOptions> session_options_;
        std::unique_ptr<Ort::Session> session_;
        std::unique_ptr<Ort::MemoryInfo> memory_info_;
    };
} // namespace OpenMS