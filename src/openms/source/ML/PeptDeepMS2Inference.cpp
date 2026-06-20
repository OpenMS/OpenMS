// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include "ONNXEnvironment.h"
#include "AminoAcidVocabulary.h"
#include <stdexcept>
#include <algorithm>
#include <string>

namespace OpenMS {

struct PeptDeepMS2Inference::Impl
{
    Ort::SessionOptions session_options_;
    std::unique_ptr<Ort::Session> session_;
    Ort::MemoryInfo memory_info_ = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    Impl(const std::string& model_path)
    {
        try {
            session_options_.SetIntraOpNumThreads(1);
            session_options_.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);

#ifdef _WIN32
            std::wstring w_model_path(model_path.begin(), model_path.end());
            session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), w_model_path.c_str(), session_options_);
#else
            session_ = std::make_unique<Ort::Session>(getONNXEnvironment(), model_path.c_str(), session_options_);
#endif
        } catch (const Ort::Exception& e) {
            throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, model_path);
        }
    }

    std::vector<float> predictMS2(const std::string& peptide_sequence, float charge, float nce, int64_t instrument_index)
    {
        int64_t seq_length = peptide_sequence.size();
        int64_t batch_size = 1;

        if (seq_length == 0) {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide sequence cannot be empty.");
        }

        const std::string valid_aas = "ACDEFGHIKLMNPQRSTVWY";
        if (peptide_sequence.find_first_not_of(valid_aas) != std::string::npos) {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unsupported residue encountered.", peptide_sequence);
        }

        // --- 1. PADDED SEQUENCE ---
        int64_t padded_length = seq_length + 2;
        std::vector<int64_t> aa_indices;
        aa_indices.reserve(padded_length);

        aa_indices.push_back(0); // N-term padding
        for (char aa : peptide_sequence) {
            aa_indices.push_back(ML::getAAIndex(aa));
        }
        aa_indices.push_back(0); // C-term padding

        std::vector<int64_t> aa_shape = {batch_size, padded_length};

        // --- 2. SCALED FEATURES (Timo's Fix) ---
        // charge_factor = 0.1, NCE_factor = 0.01
        std::vector<float> charge_data = {charge * 0.1f};
        std::vector<float> nce_data = {nce * 0.01f};
        std::vector<int64_t> inst_data = {instrument_index};

        std::vector<int64_t> charge_shape = {batch_size, 1};
        std::vector<int64_t> nce_shape = {batch_size, 1};
        std::vector<int64_t> inst_shape = {batch_size, 1};
        std::vector<int64_t> mod_shape;
        std::vector<float> mod_x_data;

        Ort::AllocatorWithDefaultOptions ort_alloc;
        size_t input_count = session_->GetInputCount();
        std::vector<std::string> input_names;
        std::vector<const char*> input_names_chars;

        for (size_t i = 0; i < input_count; i++) {
            Ort::AllocatedStringPtr name_ptr = session_->GetInputNameAllocated(i, ort_alloc);
            input_names.push_back(name_ptr.get());
        }

        for (size_t i = 0; i < input_count; i++) {
            std::string name = input_names[i];
            auto shape = session_->GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();

            if (name == "mod_x") {
                mod_shape = shape;
                if (mod_shape.empty()) mod_shape = {batch_size, padded_length, 109};
                mod_shape[0] = batch_size;
                if (mod_shape.size() > 1) mod_shape[1] = padded_length;

                int64_t total_mod = 1;
                for (size_t d = 0; d < mod_shape.size(); ++d) {
                    if (mod_shape[d] < 0) mod_shape[d] = 109; // fallback
                    total_mod *= mod_shape[d];
                }
                mod_x_data.resize(total_mod, 0.0f);
            } else if (name == "charges" || name == "charge") {
                charge_shape = shape;
                if (charge_shape.empty()) charge_shape = {batch_size, 1};
                charge_shape[0] = batch_size;
                for (size_t d = 0; d < charge_shape.size(); ++d) if (charge_shape[d] < 0) charge_shape[d] = 1;
            } else if (name == "nces" || name == "nce") {
                nce_shape = shape;
                if (nce_shape.empty()) nce_shape = {batch_size, 1};
                nce_shape[0] = batch_size;
                for (size_t d = 0; d < nce_shape.size(); ++d) if (nce_shape[d] < 0) nce_shape[d] = 1;
            } else if (name == "instrument_indices" || name == "instrument") {
                inst_shape = shape;
                if (inst_shape.empty()) inst_shape = {batch_size, 1};
                inst_shape[0] = batch_size;
                for (size_t d = 0; d < inst_shape.size(); ++d) if (inst_shape[d] < 0) inst_shape[d] = 1;
            }
        }

        std::vector<Ort::Value> input_tensors;
        input_tensors.reserve(input_count);

        for (size_t i = 0; i < input_count; i++) {
            input_names_chars.push_back(input_names[i].c_str());
            std::string name = input_names[i];

            if (name == "aa_indices") {
                input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                    memory_info_, aa_indices.data(), aa_indices.size(), aa_shape.data(), aa_shape.size()));
            } else if (name == "mod_x") {
                input_tensors.push_back(Ort::Value::CreateTensor<float>(
                    memory_info_, mod_x_data.data(), mod_x_data.size(), mod_shape.data(), mod_shape.size()));
            } else if (name == "charges" || name == "charge") {
                input_tensors.push_back(Ort::Value::CreateTensor<float>(
                    memory_info_, charge_data.data(), charge_data.size(), charge_shape.data(), charge_shape.size()));
            } else if (name == "nces" || name == "nce") {
                input_tensors.push_back(Ort::Value::CreateTensor<float>(
                    memory_info_, nce_data.data(), nce_data.size(), nce_shape.data(), nce_shape.size()));
            } else if (name == "instrument_indices" || name == "instrument") {
                input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
                    memory_info_, inst_data.data(), inst_data.size(), inst_shape.data(), inst_shape.size()));
            } else {
                throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
            }
        }

        Ort::AllocatedStringPtr out_name = session_->GetOutputNameAllocated(0, ort_alloc);
        const char* output_names[] = { out_name.get() };

        auto output_tensors = session_->Run(
            Ort::RunOptions{nullptr}, input_names_chars.data(), input_tensors.data(), input_tensors.size(), output_names, 1
        );

        // --- 3. TIMO'S POST-PROCESSING ---
        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();

        std::vector<float> intensities;
        intensities.reserve(output_count);

        // 1. Find the apex (max value)
        float apex = 0.0f;
        for (size_t i = 0; i < output_count; ++i) {
            float val = floatarr[i];
            intensities.push_back(val);
            if (val > apex) {
                apex = val;
            }
        }

        // 2. Guard against zero/negative apex
        if (apex <= 0.0f) {
            apex = 1.0f;
        }

        // 3. Normalize and floor noise
        for (float& val : intensities) {
            val /= apex;
            if (val < 1e-4f) {
                val = 0.0f;
            }
        }

        return intensities;
    }
};

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path)
    : pimpl_(std::make_unique<Impl>(model_path))
{}

PeptDeepMS2Inference::~PeptDeepMS2Inference() = default;

std::vector<float> PeptDeepMS2Inference::predictMS2(const std::string& peptide_sequence, float charge, float nce, int64_t instrument_index)
{
    return pimpl_->predictMS2(peptide_sequence, charge, nce, instrument_index);
}

} // namespace OpenMS
