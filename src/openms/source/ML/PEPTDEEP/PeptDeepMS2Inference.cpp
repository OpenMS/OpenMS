// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <map>

namespace {
    constexpr float MAX_PREDICTED_INTENSITY = 10.0f; // Cap artifactual spikes from the ONNX model
    constexpr float MIN_INTENSITY_THRESHOLD = 1e-4f; // Floor near-zero noise after normalization
}

namespace OpenMS {

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path, int intra_op_threads, size_t batch_size)
    : model_(model_path, intra_op_threads), batch_size_(batch_size)
{
    if (batch_size_ == 0) {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Batch size cannot be zero.");
    }
}

PeptDeepMS2Inference::~PeptDeepMS2Inference() = default;

std::vector<std::vector<float>> PeptDeepMS2Inference::predictMS2(
    const std::vector<std::string>& peptides,
    const std::vector<float>& charges,
    const std::vector<float>& nces,
    const std::vector<int64_t>& instrument_indices)
{
    if (peptides.empty()) {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
    }
    if (charges.size() != peptides.size() || nces.size() != peptides.size() || instrument_indices.size() != peptides.size()) {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input vectors must have the same size.");
    }

    // 1. Grouping by Sequence Length
    ML::PeptDeepInputConfig input_config;
    const std::vector<int64_t> input_shape = model_.getInputShape(0);
    if (input_shape.size() >= 2 && input_shape[1] > 0) {
        input_config.fixed_sequence_length = static_cast<size_t>(input_shape[1]);
    }

    std::vector<std::vector<size_t>> groups;
    if (input_config.fixed_sequence_length > 0) {
        groups.push_back({});
        groups.back().reserve(peptides.size());
        for (size_t i = 0; i < peptides.size(); ++i) {
            groups.back().push_back(i);
        }
    } else {
        std::map<size_t, std::vector<size_t>> indices_by_encoded_length;
        for (size_t i = 0; i < peptides.size(); ++i) {
            indices_by_encoded_length[OpenMS::AASequence::fromString(peptides[i]).size() + 2].push_back(i);
        }
        for (auto& item : indices_by_encoded_length) {
            groups.push_back(std::move(item.second));
        }
    }

    // 2. Setup master predictions array and ONNX variables
    std::vector<std::vector<float>> predictions(peptides.size());

    std::vector<std::string> input_names = model_.getInputNames();
    std::vector<const char*> input_names_chars;
    input_names_chars.reserve(input_names.size());
    for (const std::string& name : input_names) {
        input_names_chars.push_back(name.c_str());
    }

    Ort::AllocatorWithDefaultOptions ort_alloc;
    Ort::AllocatedStringPtr out_name = model_.session().GetOutputNameAllocated(0, ort_alloc);
    const char* output_names[] = { out_name.get() };

    // 3. Sub-chunking and Execution
    for (const std::vector<size_t>& group_indices : groups)
    {
        for (size_t chunk_start = 0; chunk_start < group_indices.size(); chunk_start += batch_size_)
        {
            size_t current_chunk_size = std::min(batch_size_, group_indices.size() - chunk_start);

            std::vector<std::string> chunk_peptides;
            std::vector<float> chunk_charges, chunk_nces;
            std::vector<int64_t> chunk_instruments;

            chunk_peptides.reserve(current_chunk_size);
            chunk_charges.reserve(current_chunk_size);
            chunk_nces.reserve(current_chunk_size);
            chunk_instruments.reserve(current_chunk_size);

            for (size_t j = 0; j < current_chunk_size; ++j) {
                size_t original_idx = group_indices[chunk_start + j];
                chunk_peptides.push_back(peptides[original_idx]);
                chunk_charges.push_back(charges[original_idx]);
                chunk_nces.push_back(nces[original_idx]);
                chunk_instruments.push_back(instrument_indices[original_idx]);
            }

            ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildProductMetaBatch(
                chunk_peptides, chunk_charges, chunk_nces, chunk_instruments, input_config);

            const int64_t batch_size_cast = static_cast<int64_t>(batch.batch_size);
            const int64_t sequence_length = static_cast<int64_t>(batch.sequence_length);

            std::vector<Ort::Value> input_tensors;
            input_tensors.reserve(input_names.size());

            for (size_t i = 0; i < input_names.size(); ++i) {
                const std::string& name = input_names[i];
                Ort::TypeInfo type_info = model_.session().GetInputTypeInfo(i);
                auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
                std::vector<int64_t> expected_shape = tensor_info.GetShape();
                if (!expected_shape.empty()) {
                    expected_shape[0] = batch_size_cast;
                }

                if (name == "aa_indices" || name == "input_sequences") {
                    expected_shape = {batch_size_cast, sequence_length};
                    input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                        model_.memoryInfo(), batch.aa_indices.data(), batch.aa_indices.size(), expected_shape.data(), expected_shape.size()));
                } else if (name == "mod_x") {
                    expected_shape = {batch_size_cast, sequence_length, ML::PEPTDEEP_MOD_ELEMENTS};
                    input_tensors.push_back(Ort::Value::CreateTensor<float>(
                        model_.memoryInfo(), batch.mod_x.data(), batch.mod_x.size(), expected_shape.data(), expected_shape.size()));
                } else if (name == "charges" || name == "charge") {
                    expected_shape = {batch_size_cast, 1};
                    input_tensors.push_back(Ort::Value::CreateTensor<float>(
                        model_.memoryInfo(), batch.charges.data(), batch.charges.size(), expected_shape.data(), expected_shape.size()));
                } else if (name == "nces" || name == "nce") {
                    expected_shape = {batch_size_cast, 1};
                    input_tensors.push_back(Ort::Value::CreateTensor<float>(
                        model_.memoryInfo(), batch.nces.data(), batch.nces.size(), expected_shape.data(), expected_shape.size()));
                } else if (name == "instrument_indices" || name == "instrument") {
                    expected_shape = {batch_size_cast};
                    input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                        model_.memoryInfo(), batch.instrument_indices.data(), batch.instrument_indices.size(), expected_shape.data(), expected_shape.size()));
                } else {
                    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
                }
            }

            auto output_tensors = model_.session().Run(
                Ort::RunOptions{nullptr}, input_names_chars.data(), input_tensors.data(), input_tensors.size(), output_names, 1);

            // 4. Output Processing and MS2 Base Peak Normalization
            float* floatarr = output_tensors.front().GetTensorMutableData<float>();
            auto out_shape = output_tensors.front().GetTensorTypeAndShapeInfo().GetShape();

            if (out_shape.size() != 3) {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected MS2 output tensor to have a rank of 3.", std::to_string(out_shape.size()));
            }

            int64_t actual_rows = out_shape[1];
            int64_t num_ions = out_shape.back();
            int64_t elements_per_batch = actual_rows * num_ions;

            for (size_t j = 0; j < current_chunk_size; ++j) {
                const std::string& current_peptide = chunk_peptides[j];

                size_t parsed_size = OpenMS::AASequence::fromString(current_peptide).size();
                if (parsed_size < 1) {
                    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Parsed peptide sequence cannot be empty.");
                }

                int64_t valid_fragments = static_cast<int64_t>(parsed_size) - 1;

                int64_t start_row = (actual_rows > valid_fragments) ? (actual_rows - valid_fragments) / 2 : 0;

                std::vector<float> intensities;
                intensities.reserve(valid_fragments * num_ions);

                float apex = 0.0f;
                float* batch_floatarr = floatarr + (j * elements_per_batch);

                for (int64_t i = 0; i < valid_fragments; ++i) {
                    int64_t row_idx = start_row + i;
                    if (row_idx >= actual_rows) break;

                    for (int64_t k_ion = 0; k_ion < num_ions; ++k_ion) {
                        float val = batch_floatarr[(row_idx * num_ions) + k_ion];

                        if (val > MAX_PREDICTED_INTENSITY) {
                            val = 0.0f;
                        }
                        intensities.push_back(val < 0.0f ? 0.0f : val);

                        if (val > apex) apex = val;
                    }
                }

                if (apex <= 0.0f) apex = 1.0f;
                for (float& val : intensities) {
                    val /= apex;
                    if (val < MIN_INTENSITY_THRESHOLD) val = 0.0f;
                }

                predictions[group_indices[chunk_start + j]] = std::move(intensities);
            }
        }
    }

    return predictions;
}

} // namespace OpenMS