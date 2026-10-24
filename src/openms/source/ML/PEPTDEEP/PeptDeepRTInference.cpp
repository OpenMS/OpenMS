// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include "PeptDeepModelRuns.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <onnxruntime_cxx_api.h>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace OpenMS
{
    PeptDeepRTInference::PeptDeepRTInference(const std::string& model_path, int threads, size_t batch_size)
        : model_(model_path, threads), batch_size_(batch_size)
    {
      if (batch_size_ == 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Batch size cannot be zero.");
      }
    }

    PeptDeepRTInference::~PeptDeepRTInference() = default;

    std::vector<float> PeptDeepRTInference::predictRT(const std::vector<std::string>& peptides)
    {
        if (peptides.empty())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
        }

        // The peptides in groups of one length, and the model runs of batches of at most batch_size_ of them, side by side
        ML::Internal::PreparedPeptides prepared = ML::Internal::preparePeptides(model_, peptides);
        const ML::Internal::ModelRuns runs(model_, {ML::Internal::InputKind::AA_INDICES, ML::Internal::InputKind::MOD_X});
        std::vector<float> predictions(prepared.peptides.size(), 0.0f);
        model_.forEachBatch(prepared.groups, batch_size_, [&](const size_t* indices, size_t n)
        {
            std::vector<OpenMS::AASequence> chunk_peptides;
            chunk_peptides.reserve(n);
            for (size_t j = 0; j < n; ++j)
            {
                chunk_peptides.push_back(prepared.peptides[indices[j]]);
            }
            ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildPeptideBatch(chunk_peptides, prepared.config);
            auto output_tensors = runs.run(batch);
            float* floatarr = output_tensors.front().GetTensorMutableData<float>();
            for (size_t j = 0; j < n; ++j)
            {
                predictions[indices[j]] = floatarr[j];
            }
        });

        return predictions;
    }

} // namespace OpenMS