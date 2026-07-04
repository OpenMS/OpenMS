// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/ONNXPredictorBase.h>
#include <OpenMS/ML/ONNXEnvironment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <onnxruntime_cxx_api.h>

#ifdef _WIN32
#include <windows.h>
#endif

namespace OpenMS
{
    ONNXPredictorBase::ONNXPredictorBase(const std::string& model_path)
    {
        try {
            session_options_ = std::make_unique<Ort::SessionOptions>();

            // Removed thread limit (per PR review) to allow ONNX default parallelization
            session_options_->SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

            memory_info_ = std::make_unique<Ort::MemoryInfo>(
                Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)
            );

#ifdef _WIN32
            if (model_path.empty()) {
                session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), L"", *session_options_);
            } else {
                int size_needed = MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.length(), NULL, 0);
                std::wstring w_model_path(size_needed, 0);
                MultiByteToWideChar(CP_UTF8, 0, model_path.c_str(), (int)model_path.length(), &w_model_path[0], size_needed);
                session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), w_model_path.c_str(), *session_options_);
            }
#else
            session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), *session_options_);
#endif
        } catch (const Ort::Exception& e) {
            throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, model_path);
        }
    }

    ONNXPredictorBase::~ONNXPredictorBase() = default;

} // namespace OpenMS