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
#include <algorithm>
#include <exception>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace OpenMS
{
    ONNXPredictorBase::ONNXPredictorBase(const std::string& model_path, int threads) :
        threads_(std::max(1, threads))
    {
        try {
            session_options_ = std::make_unique<Ort::SessionOptions>();

            // One intra-op thread: the predictors run their model runs side by side (forEachCall()).
            // Split over 16 intra-op threads, a run of a batch of PeptDeep's MS2 model kept them busy
            // waiting for each other (all cores used, 2.2x the wall time of 16 runs at once on iPRG2015).
            session_options_->SetIntraOpNumThreads(1);
            session_options_->SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
            // Without ONNX Runtime's memory pattern. The inputs of the PeptDeep models change shape
            // from call to call (batch size, peptide length); for each shape, the pattern plans one
            // block for all intermediate tensors of a Run, which every Run of that shape allocates and
            // frees again. Allocated tensor by tensor instead, the results are the same, inference is
            // as fast and the peak memory is lower. Where ONNX Runtime has no CPU arena (none on
            // 32-bit platforms), these blocks come from the system allocator, and one that does not
            // reuse such large blocks keeps growing the memory of the process.
            session_options_->DisableMemPattern();

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

    int ONNXPredictorBase::threads() const
    {
        return threads_;
    }

    void ONNXPredictorBase::forEachCall(size_t n, const std::function<void(size_t)>& call) const
    {
        std::exception_ptr failure; // exceptions must not leave the parallel region
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, std::min(threads_, static_cast<int>(std::min<size_t>(n, 1024))))) if (n > 1 && threads_ > 1)
        for (std::ptrdiff_t k = 0; k < static_cast<std::ptrdiff_t>(n); ++k)
        {
            try
            {
                call(static_cast<size_t>(k));
            }
            catch (...)
            {
#pragma omp critical (ONNXPredictorBase_forEachCall)
                if (!failure) { failure = std::current_exception(); }
            }
        }
        if (failure) { std::rethrow_exception(failure); }
    }

} // namespace OpenMS
