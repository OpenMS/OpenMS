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
#include <onnxruntime_cxx_api.h>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace OpenMS {

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path)
    : model_(model_path)
{}

PeptDeepMS2Inference::~PeptDeepMS2Inference() = default;

std::vector<std::vector<float>> PeptDeepMS2Inference::predictMS2(
    const std::vector<std::string>& peptides,
    const std::vector<float>& charges,
    const std::vector<float>& nces,
    const std::vector<int64_t>& instrument_indices)
{
    ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildUnmodifiedInstrumentBatch(
        peptides,
        charges,
        nces,
        instrument_indices);
    const int64_t batch_size = static_cast<int64_t>(batch.batch_size);
    const int64_t padded_length = static_cast<int64_t>(batch.sequence_length);

    // 4. Dynamically map model inputs and their EXACT expected shapes
    Ort::AllocatorWithDefaultOptions ort_alloc;
    size_t input_count = model_.getInputCount();

    std::vector<std::string> input_names = model_.getInputNames();
    std::vector<const char*> input_names_chars;
    input_names_chars.reserve(input_count);
    std::vector<Ort::Value> input_tensors;
    input_tensors.reserve(input_count);

    for (size_t i = 0; i < input_count; i++) {
        const std::string& name = input_names[i];
        input_names_chars.push_back(name.c_str());

        // Ask ONNX what shape it expects for this specific input
        Ort::TypeInfo type_info = model_.session().GetInputTypeInfo(i);
        auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> expected_shape = tensor_info.GetShape();

        // Patch the dynamic batch dimension to match our actual batch size
        if (!expected_shape.empty()) {
            expected_shape[0] = batch_size;
        }

        if (name == "aa_indices") {
            expected_shape = {batch_size, padded_length}; // Explicitly enforce sequence length
            input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(model_.memoryInfo(), batch.aa_indices.data(), batch.aa_indices.size(), expected_shape.data(), expected_shape.size()));
        } else if (name == "mod_x") {
            expected_shape = {batch_size, padded_length, ML::PEPTDEEP_MOD_ELEMENTS};
            input_tensors.push_back(Ort::Value::CreateTensor<float>(model_.memoryInfo(), batch.mod_x.data(), batch.mod_x.size(), expected_shape.data(), expected_shape.size()));
        } else if (name == "charges" || name == "charge") {
            input_tensors.push_back(Ort::Value::CreateTensor<float>(model_.memoryInfo(), batch.charges.data(), batch.charges.size(), expected_shape.data(), expected_shape.size()));
        } else if (name == "nces" || name == "nce") {
            input_tensors.push_back(Ort::Value::CreateTensor<float>(model_.memoryInfo(), batch.nces.data(), batch.nces.size(), expected_shape.data(), expected_shape.size()));
        } else if (name == "instrument_indices" || name == "instrument") {
            input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(model_.memoryInfo(), batch.instrument_indices.data(), batch.instrument_indices.size(), expected_shape.data(), expected_shape.size()));
        } else {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
        }
    }

    Ort::AllocatedStringPtr out_name = model_.session().GetOutputNameAllocated(0, ort_alloc);
    const char* output_names[] = { out_name.get() };

    auto output_tensors = model_.session().Run(
        Ort::RunOptions{nullptr}, input_names_chars.data(), input_tensors.data(), input_tensors.size(), output_names, 1
    );

    float* floatarr = output_tensors.front().GetTensorMutableData<float>();
    auto out_shape = output_tensors.front().GetTensorTypeAndShapeInfo().GetShape();

    int64_t actual_rows = out_shape.size() >= 3 ? out_shape[1] : (out_shape.size() >= 2 ? out_shape[out_shape.size() - 2] : 1);
    int64_t num_ions = out_shape.back();
    int64_t elements_per_batch = actual_rows * num_ions;

    std::vector<std::vector<float>> batch_intensities;
    batch_intensities.reserve(batch_size);

    // 5. Extract and normalize spectra for EACH peptide individually
    for (size_t k = 0; k < static_cast<size_t>(batch_size); ++k) {
        int64_t valid_fragments = peptides[k].size() - 1;
        int64_t start_row = (actual_rows > valid_fragments) ? (actual_rows - valid_fragments) / 2 : 0;

        std::vector<float> intensities;
        intensities.reserve(valid_fragments * num_ions);

        float apex = 0.0f;
        float* batch_floatarr = floatarr + (k * elements_per_batch);

        for (int64_t i = 0; i < valid_fragments; ++i) {
            int64_t row_idx = start_row + i;
            if (row_idx >= actual_rows) break;

            for (int64_t j = 0; j < num_ions; ++j) {
                float val = batch_floatarr[(row_idx * num_ions) + j];

                // Capping massive padding logits to prevent Base Peak Normalization failure
                if (val > 10.0f) {
                    val = 0.0f;
                }
                intensities.push_back(val < 0.0f ? 0.0f : val);

                if (val > apex) apex = val;
            }
        }

        if (apex <= 0.0f) apex = 1.0f;
        for (float& val : intensities) {
            val /= apex;
            if (val < 1e-4f) val = 0.0f;
        }
        if (!intensities.empty()) intensities[0] = 0.0f;

        batch_intensities.push_back(std::move(intensities));
    }

    return batch_intensities;
}

} // namespace OpenMS
