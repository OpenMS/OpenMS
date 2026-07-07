// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PeptDeepInput.h>
#include <onnxruntime_cxx_api.h>
#include <map>
#include <stdexcept>
#include <string>

namespace OpenMS
{
    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path)
        : model_(model_path)
    {}

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        if (peptides.empty())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
        }

        ML::PeptDeepInputConfig input_config;
        const std::vector<int64_t> input_shape = model_.getInputShape(0);
        if (input_shape.size() >= 2 && input_shape[1] > 0)
        {
            input_config.fixed_sequence_length = static_cast<size_t>(input_shape[1]);
        }

        std::vector<std::vector<size_t>> groups;
        if (input_config.fixed_sequence_length > 0)
        {
            groups.push_back({});
            groups.back().reserve(peptides.size());
            for (size_t i = 0; i < peptides.size(); ++i)
            {
                groups.back().push_back(i);
            }
        }
        else
        {
            std::map<size_t, std::vector<size_t>> indices_by_encoded_length;
            for (size_t i = 0; i < peptides.size(); ++i)
            {
                indices_by_encoded_length[peptides[i].size() + 2].push_back(i);
            }
            for (auto& item : indices_by_encoded_length)
            {
                groups.push_back(std::move(item.second));
            }
        }

        std::vector<float> predictions(peptides.size(), 0.0f);

        for (const std::vector<size_t>& group_indices : groups)
        {
            std::vector<std::string> group_peptides;
            group_peptides.reserve(group_indices.size());
            for (size_t idx : group_indices)
            {
                group_peptides.push_back(peptides[idx]);
            }

            ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildUnmodifiedPeptideBatch(group_peptides, input_config);
            std::vector<int64_t> seq_shape = {
                static_cast<int64_t>(batch.batch_size),
                static_cast<int64_t>(batch.sequence_length)
            };

            Ort::TypeInfo mod_type_info = model_.session().GetInputTypeInfo(1);
            auto mod_tensor_info = mod_type_info.GetTensorTypeAndShapeInfo();
            std::vector<int64_t> mod_shape = mod_tensor_info.GetShape();

            if (mod_shape.size() < 2) {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "PeptDeep RT model input 'mod_x' must have at least 2 dimensions.",
                    std::to_string(mod_shape.size()));
            }

            mod_shape[0] = static_cast<int64_t>(batch.batch_size);
            mod_shape[1] = static_cast<int64_t>(batch.sequence_length);

            std::vector<Ort::Value> input_tensors;

            input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                model_.memoryInfo(), batch.aa_indices.data(), batch.aa_indices.size(), seq_shape.data(), seq_shape.size()
            ));

            input_tensors.push_back(Ort::Value::CreateTensor<float>(
                model_.memoryInfo(), batch.mod_x.data(), batch.mod_x.size(), mod_shape.data(), mod_shape.size()
            ));

            Ort::AllocatorWithDefaultOptions ort_alloc;
            Ort::AllocatedStringPtr seq_name_ptr = model_.session().GetInputNameAllocated(0, ort_alloc);
            Ort::AllocatedStringPtr mod_name_ptr = model_.session().GetInputNameAllocated(1, ort_alloc);
            Ort::AllocatedStringPtr output_name_ptr = model_.session().GetOutputNameAllocated(0, ort_alloc);

            const char* input_names[] = { seq_name_ptr.get(), mod_name_ptr.get() };
            const char* output_names[] = { output_name_ptr.get() };

            auto output_tensors = model_.session().Run(
                Ort::RunOptions{nullptr}, input_names, input_tensors.data(), 2, output_names, 1
            );

            size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();
            if (output_count < batch.batch_size) {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "ONNX model output shape mismatch.",
                    std::to_string(output_count));
            }

            float* floatarr = output_tensors.front().GetTensorMutableData<float>();
            for (size_t i = 0; i < group_indices.size(); ++i)
            {
                predictions[group_indices[i]] = floatarr[i];
            }
        }

        return predictions;
    }

} // namespace OpenMS
