// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PeptDeepUtils.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>
#include <string>

namespace OpenMS
{
    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path)
        : ONNXPredictorBase(model_path) // Handled by base class
    {}

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        if (peptides.empty()) {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide list cannot be empty.");
        }

        // 1. Tokenize peptides
        std::vector<int64_t> input_tokens;
        input_tokens.reserve(peptides.size() * ML::PEPTDEEP_MAX_SEQUENCE_LENGTH);

        for (const auto& p : peptides) {
            ML::validatePeptide(p);
            for (size_t i = 0; i < ML::PEPTDEEP_MAX_SEQUENCE_LENGTH; ++i) {
                if (i < p.length()) {
                    input_tokens.push_back(ML::getAAIndex(p[i]));
                } else {
                    input_tokens.push_back(0); // Pad sequence
                }
            }
        }

        std::vector<int64_t> seq_shape = { static_cast<int64_t>(peptides.size()), static_cast<int64_t>(ML::PEPTDEEP_MAX_SEQUENCE_LENGTH) };

        // 2. Fetch expected mod_shape dynamically from ONNX
        Ort::TypeInfo mod_type_info = session_->GetInputTypeInfo(1);
        auto mod_tensor_info = mod_type_info.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> mod_shape = mod_tensor_info.GetShape();

        if (mod_shape.size() < 2) {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "PeptDeep RT model input 'mod_x' must have at least 2 dimensions.",
                std::to_string(mod_shape.size()));
        }

        mod_shape[0] = peptides.size();
        mod_shape[1] = ML::PEPTDEEP_MAX_SEQUENCE_LENGTH;

        // 3. Shared utility generates the 109-element tensor for the batch
        std::vector<float> mod_x_data = ML::generateUnmodifiedModXTensor(peptides.size(), ML::PEPTDEEP_MAX_SEQUENCE_LENGTH);

        std::vector<Ort::Value> input_tensors;

        input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
            *memory_info_, input_tokens.data(), input_tokens.size(), seq_shape.data(), seq_shape.size()
        ));

        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            *memory_info_, mod_x_data.data(), mod_x_data.size(), mod_shape.data(), mod_shape.size()
        ));

        Ort::AllocatorWithDefaultOptions ort_alloc;
        Ort::AllocatedStringPtr seq_name_ptr = session_->GetInputNameAllocated(0, ort_alloc);
        Ort::AllocatedStringPtr mod_name_ptr = session_->GetInputNameAllocated(1, ort_alloc);
        Ort::AllocatedStringPtr output_name_ptr = session_->GetOutputNameAllocated(0, ort_alloc);

        const char* input_names[] = { seq_name_ptr.get(), mod_name_ptr.get() };
        const char* output_names[] = { output_name_ptr.get() };

        // 4. Run Inference
        auto output_tensors = session_->Run(
            Ort::RunOptions{nullptr}, input_names, input_tensors.data(), 2, output_names, 1
        );

        size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();
        if (output_count < peptides.size()) {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "ONNX model output shape mismatch.",
                std::to_string(output_count));
        }

        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        return std::vector<float>(floatarr, floatarr + peptides.size());
    }

} // namespace OpenMS
