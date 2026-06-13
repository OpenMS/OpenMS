// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include "ONNXEnvironment.h"
#include "AminoAcidVocabulary.h"
#include <stdexcept>

namespace OpenMS
{
    // 1. THE HIDDEN IMPLEMENTATION STRUCT
    struct PeptDeepRTInference::Impl
    {
        Ort::SessionOptions session_options_;
        std::unique_ptr<Ort::Session> session_;

        // Constructor
        Impl(const std::string& model_path)
        {
            try
            {
                session_options_.SetIntraOpNumThreads(1);
                session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
                session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), session_options_);
            }
            catch (const Ort::Exception& e)
            {
                throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, model_path);
            }
        }

        // Helper Function (Moved from public class to hidden Impl)
        std::vector<int64_t> tokenizePeptides(const std::vector<std::string>& peptides, size_t& max_length)
        {
            max_length = 132;
            std::vector<int64_t> flat_tokens;
            flat_tokens.reserve(peptides.size() * max_length);

            for (const auto& p : peptides) {
                if (p.length() > max_length) {
                    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                        "Peptide sequence exceeds the maximum allowed length of 132 residues.",
                        std::to_string(p.length()));
                }

                for (size_t i = 0; i < max_length; ++i) {
                    if (i < p.length()) {
                        flat_tokens.push_back(ML::getAAIndex(p[i]));
                    } else {
                        flat_tokens.push_back(0); // Padding token
                    }
                }
            }
            return flat_tokens;
        }

        // Core Prediction Logic
        std::vector<float> predictRT(const std::vector<std::string>& peptides)
        {
            if (peptides.empty())
            {
                throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide list cannot be empty.");
            }

            size_t max_seq_len = 132;
            std::vector<int64_t> input_tokens = tokenizePeptides(peptides, max_seq_len);
            std::vector<int64_t> seq_shape = { static_cast<int64_t>(peptides.size()), static_cast<int64_t>(max_seq_len) };

            Ort::TypeInfo mod_type_info = session_->GetInputTypeInfo(1);
            auto mod_tensor_info = mod_type_info.GetTensorTypeAndShapeInfo();
            std::vector<int64_t> mod_shape = mod_tensor_info.GetShape();

            // Guard against unexpected lower-rank shapes
            if (mod_shape.size() < 2)
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "PeptDeep RT model input 'mod_x' must have at least 2 dimensions.",
                    std::to_string(mod_shape.size()));
            }

            mod_shape[0] = peptides.size();
            mod_shape[1] = max_seq_len;

            int64_t total_mod_elements = 1;
            for (int64_t dim : mod_shape)
            {
                // If a dynamic dimension (-1) is encountered, fallback to standard AlphaPeptDeep feature width
                if (dim < 0)
                {
                    total_mod_elements *= 109;
                }
                else
                {
                    total_mod_elements *= dim;
                }
            }

            std::vector<float> mod_x_data(total_mod_elements, 0.0f);
            Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
            std::vector<Ort::Value> input_tensors;

            input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                memory_info, input_tokens.data(), input_tokens.size(), seq_shape.data(), seq_shape.size()
            ));

            input_tensors.push_back(Ort::Value::CreateTensor<float>(
                memory_info, mod_x_data.data(), mod_x_data.size(), mod_shape.data(), mod_shape.size()
            ));

            Ort::AllocatorWithDefaultOptions ort_alloc;
            Ort::AllocatedStringPtr seq_name_ptr = session_->GetInputNameAllocated(0, ort_alloc);
            Ort::AllocatedStringPtr mod_name_ptr = session_->GetInputNameAllocated(1, ort_alloc);
            Ort::AllocatedStringPtr output_name_ptr = session_->GetOutputNameAllocated(0, ort_alloc);

            const char* input_names[] = { seq_name_ptr.get(), mod_name_ptr.get() };
            const char* output_names[] = { output_name_ptr.get() };

            auto output_tensors = session_->Run(
                Ort::RunOptions{nullptr}, input_names, input_tensors.data(), 2, output_names, 1
            );

            size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();
            if (output_count < peptides.size()) {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "ONNX model output shape mismatch: returned fewer elements than requested.",
                    std::to_string(output_count) + " < " + std::to_string(peptides.size()));
            }

            float* floatarr = output_tensors.front().GetTensorMutableData<float>();
            return std::vector<float>(floatarr, floatarr + peptides.size());
        }
    };

    // 2. THE PUBLIC WRAPPER (Forwarding calls to Impl)

    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path)
        : pimpl_(std::make_unique<Impl>(model_path))
    {}

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        return pimpl_->predictRT(peptides);
    }

} // namespace OpenMS