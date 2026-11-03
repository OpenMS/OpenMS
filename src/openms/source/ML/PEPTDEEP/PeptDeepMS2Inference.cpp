// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include "PeptDeepModelRuns.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>
#include <algorithm>
#include <string>

namespace {
    constexpr float MAX_PREDICTED_INTENSITY = 10.0f; // Cap artifactual spikes from the ONNX model
    constexpr float MIN_INTENSITY_THRESHOLD = 1e-4f; // Floor near-zero noise after normalization
}

namespace OpenMS {

PeptDeepMS2Inference::PeptDeepMS2Inference(const std::string& model_path, int threads, size_t batch_size)
    : model_(model_path, threads), batch_size_(batch_size)
{
  if (batch_size_ == 0)
  {
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
    ML::validateInstrumentIndices(instrument_indices);

    // The peptides in groups of one length, and the model runs of batches of at most batch_size_ of them, side by side
    ML::Internal::PreparedPeptides prepared = ML::Internal::preparePeptides(model_, peptides);
    const ML::Internal::ModelRuns runs(model_, {ML::Internal::InputKind::AA_INDICES, ML::Internal::InputKind::MOD_X,
                                               ML::Internal::InputKind::CHARGES, ML::Internal::InputKind::NCES,
                                               ML::Internal::InputKind::INSTRUMENTS});
    std::vector<std::vector<float>> predictions(prepared.peptides.size());
    model_.forEachBatch(prepared.groups, batch_size_, [&](const size_t* indices, size_t n)
    {
        std::vector<OpenMS::AASequence> chunk_peptides;
        std::vector<float> chunk_charges, chunk_nces;
        std::vector<int64_t> chunk_instruments;
        chunk_peptides.reserve(n);
        chunk_charges.reserve(n);
        chunk_nces.reserve(n);
        chunk_instruments.reserve(n);
        for (size_t j = 0; j < n; ++j) {
            const size_t original_idx = indices[j];
            chunk_peptides.push_back(prepared.peptides[original_idx]);
            chunk_charges.push_back(charges[original_idx]);
            chunk_nces.push_back(nces[original_idx]);
            chunk_instruments.push_back(instrument_indices[original_idx]);
        }

        ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildProductMetaBatch(
            chunk_peptides, chunk_charges, chunk_nces, chunk_instruments, prepared.config);
        const int64_t batch_size_cast = static_cast<int64_t>(batch.batch_size);
        auto output_tensors = runs.run(batch);

        // Output Processing and MS2 Base Peak Normalization
        float* floatarr = output_tensors.front().GetTensorMutableData<float>();
        auto out_shape = output_tensors.front().GetTensorTypeAndShapeInfo().GetShape();

        if (out_shape.size() != 3)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Expected MS2 output tensor to have a rank of 3.", std::to_string(out_shape.size()));
        }

        if (out_shape[0] != batch_size_cast)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ONNX MS2 output batch dimension mismatch.", std::to_string(out_shape[0]));
        }

        int64_t actual_rows = out_shape[1];
        int64_t num_ions = out_shape.back();
        int64_t elements_per_batch = actual_rows * num_ions;

        for (size_t j = 0; j < n; ++j) {
            const OpenMS::AASequence& current_seq = chunk_peptides[j];

            size_t parsed_size = current_seq.size();
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

            predictions[indices[j]] = std::move(intensities);
        }
    });

    return predictions;
}

size_t PeptDeepMS2Inference::ionIndex(size_t length, char type, size_t ordinal, int charge, size_t channels)
{
    if (length < 2 || ordinal < 1 || ordinal >= length || (charge != 1 && charge != 2) || channels < Y_Z2 + 1) {
        return NO_ION;
    }
    // cleavage position pos carries b(pos + 1) and y(length - 1 - pos): the model emits the fragments N- to C-terminal
    if (type == 'b') { return (ordinal - 1) * channels + (charge == 1 ? B_Z1 : B_Z2); }
    if (type == 'y') { return (length - 1 - ordinal) * channels + (charge == 1 ? Y_Z1 : Y_Z2); }
    return NO_ION;
}

} // namespace OpenMS