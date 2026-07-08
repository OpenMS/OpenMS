// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/ONNX/ONNXPredictorBase.h>
#include "ONNXEnvironment.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>

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

    ONNXPredictorBase::ONNXPredictorBase(ONNXPredictorBase&&) noexcept = default;

    ONNXPredictorBase& ONNXPredictorBase::operator=(ONNXPredictorBase&&) noexcept = default;

    ONNXPredictorBase::~ONNXPredictorBase() = default;

    Ort::Session& ONNXPredictorBase::session()
    {
        return *session_;
    }

    const Ort::Session& ONNXPredictorBase::session() const
    {
        return *session_;
    }

    Ort::MemoryInfo& ONNXPredictorBase::memoryInfo()
    {
        return *memory_info_;
    }

    const Ort::MemoryInfo& ONNXPredictorBase::memoryInfo() const
    {
        return *memory_info_;
    }

    std::vector<std::string> ONNXPredictorBase::getInputNames() const
    {
        Ort::AllocatorWithDefaultOptions allocator;
        std::vector<std::string> names;
        const size_t count = session_->GetInputCount();
        names.reserve(count);
        for (size_t i = 0; i < count; ++i)
        {
            Ort::AllocatedStringPtr name = session_->GetInputNameAllocated(i, allocator);
            names.emplace_back(name.get());
        }
        return names;
    }

    std::vector<std::string> ONNXPredictorBase::getOutputNames() const
    {
        Ort::AllocatorWithDefaultOptions allocator;
        std::vector<std::string> names;
        const size_t count = session_->GetOutputCount();
        names.reserve(count);
        for (size_t i = 0; i < count; ++i)
        {
            Ort::AllocatedStringPtr name = session_->GetOutputNameAllocated(i, allocator);
            names.emplace_back(name.get());
        }
        return names;
    }

    std::vector<int64_t> ONNXPredictorBase::getInputShape(size_t input_index) const
    {
        if (input_index >= session_->GetInputCount())
        {
            throw Exception::InvalidValue(
                __FILE__,
                __LINE__,
                OPENMS_PRETTY_FUNCTION,
                "ONNX input index is out of range.",
                std::to_string(input_index));
        }
        Ort::TypeInfo type_info = session_->GetInputTypeInfo(input_index);
        return type_info.GetTensorTypeAndShapeInfo().GetShape();
    }

    size_t ONNXPredictorBase::getInputCount() const
    {
        return session_->GetInputCount();
    }

    size_t ONNXPredictorBase::getOutputCount() const
    {
        return session_->GetOutputCount();
    }

} // namespace OpenMS
