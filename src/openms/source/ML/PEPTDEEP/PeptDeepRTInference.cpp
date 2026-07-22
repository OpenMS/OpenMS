// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <onnxruntime_cxx_api.h>
#include <map>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace OpenMS
{
    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path, int intra_op_threads, size_t batch_size)
        : model_(model_path, intra_op_threads), batch_size_(batch_size)
    {}

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        if (peptides.empty())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
        }

        // 1. Grouping by Sequence Length
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
                // FIX: Parse the sequence first so modification strings aren't counted as length
                indices_by_encoded_length[OpenMS::AASequence::fromString(peptides[i]).size() + 2].push_back(i);
            }
            for (auto& item : indices_by_encoded_length)
            {
                groups.push_back(std::move(item.second));
            }
        }

        std::vector<float> predictions(peptides.size(), 0.0f);

        std::vector<std::string> input_names = model_.getInputNames();
        std::vector<const char*> input_names_chars;
        input_names_chars.reserve(input_names.size());
        for (const std::string& name : input_names)
        {
            input_names_chars.push_back(name.c_str());
        }

        Ort::AllocatorWithDefaultOptions ort_alloc;
        Ort::AllocatedStringPtr out_name = model_.session().GetOutputNameAllocated(0, ort_alloc);
        const char* output_names[] = { out_name.get() };

        for (const std::vector<size_t>& group_indices : groups)
        {
            for (size_t chunk_start = 0; chunk_start < group_indices.size(); chunk_start += batch_size_)
            {
                size_t current_chunk_size = std::min(batch_size_, group_indices.size() - chunk_start);

                std::vector<std::string> chunk_peptides;
                chunk_peptides.reserve(current_chunk_size);
                for (size_t j = 0; j < current_chunk_size; ++j)
                {
                    chunk_peptides.push_back(peptides[group_indices[chunk_start + j]]);
                }

                ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildPeptideBatch(chunk_peptides, input_config);
                const int64_t batch_size_cast = static_cast<int64_t>(batch.batch_size);
                const int64_t sequence_length = static_cast<int64_t>(batch.sequence_length);

                std::vector<Ort::Value> input_tensors;
                input_tensors.reserve(input_names.size());

                for (size_t i = 0; i < input_names.size(); ++i)
                {
                    const std::string& name = input_names[i];
                    Ort::TypeInfo type_info = model_.session().GetInputTypeInfo(i);
                    auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
                    std::vector<int64_t> expected_shape = tensor_info.GetShape();
                    if (!expected_shape.empty())
                    {
                        expected_shape[0] = batch_size_cast;
                    }

                    if (name == "aa_indices" || name == "input_sequences")
                    {
                        expected_shape = {batch_size_cast, sequence_length};
                        input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                            model_.memoryInfo(), batch.aa_indices.data(), batch.aa_indices.size(), expected_shape.data(), expected_shape.size()));
                    }
                    else if (name == "mod_x")
                    {
                        expected_shape = {batch_size_cast, sequence_length, ML::PEPTDEEP_MOD_ELEMENTS};
                        input_tensors.push_back(Ort::Value::CreateTensor<float>(
                            model_.memoryInfo(), batch.mod_x.data(), batch.mod_x.size(), expected_shape.data(), expected_shape.size()));
                    }
                    else
                    {
                        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
                    }
                }

                auto output_tensors = model_.session().Run(
                    Ort::RunOptions{nullptr}, input_names_chars.data(), input_tensors.data(), input_tensors.size(), output_names, 1
                );

                float* floatarr = output_tensors.front().GetTensorMutableData<float>();

                for (size_t j = 0; j < current_chunk_size; ++j)
                {
                    predictions[group_indices[chunk_start + j]] = floatarr[j];
                }
            }
        }

        return predictions;
    }

} // namespace OpenMS