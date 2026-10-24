// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// PRIVATE OpenMS-internal header of the PeptDeep predictors (source/ML/PEPTDEEP only): like ONNXEnvironment.h, not
// installed, so that the ONNX Runtime types do not leak into the public API.

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/ONNX/ONNXPredictorBase.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <onnxruntime_cxx_api.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace OpenMS::ML::Internal
{
  /// The peptides of a prediction, parsed, and their groups: the peptides of one model run each have one length
  struct PreparedPeptides
  {
    std::vector<AASequence> peptides;
    PeptDeepInputConfig config;
    /// Indices into @p peptides: of equal encoded length (sequence and terminal tokens) for a model with a dynamic
    /// sequence length, in order of length; all in one group for a model with a fixed one
    std::vector<std::vector<size_t>> groups;
  };

  /// Validates and parses @p peptides for @p model (whose first input's shape tells a fixed sequence length) and
  /// groups them by length
  inline PreparedPeptides preparePeptides(const ONNXPredictorBase& model, const std::vector<std::string>& peptides)
  {
    PreparedPeptides prepared;
    prepared.peptides.reserve(peptides.size());
    for (const std::string& peptide : peptides)
    {
      validatePeptide(peptide);
      prepared.peptides.push_back(AASequence::fromString(peptide));
    }
    const std::vector<int64_t> input_shape = model.getInputShape(0);
    if (input_shape.size() >= 2 && input_shape[1] > 0)
    {
      prepared.config.fixed_sequence_length = static_cast<size_t>(input_shape[1]);
    }
    if (prepared.config.fixed_sequence_length > 0)
    {
      prepared.groups.emplace_back(prepared.peptides.size());
      for (size_t i = 0; i < prepared.peptides.size(); ++i) { prepared.groups.back()[i] = i; }
    }
    else
    {
      std::map<size_t, std::vector<size_t>> by_encoded_length;
      for (size_t i = 0; i < prepared.peptides.size(); ++i) { by_encoded_length[prepared.peptides[i].size() + 2].push_back(i); }
      for (auto& [length, indices] : by_encoded_length) { prepared.groups.push_back(std::move(indices)); }
    }
    return prepared;
  }

  /// The inputs of a PeptDeep model, by name
  enum class InputKind { AA_INDICES, MOD_X, CHARGES, NCES, INSTRUMENTS };

  /**
    @brief The runs of a model on batches of peptides: its inputs resolved once (by name; an input the predictor does
    not support is an error), then tensors of a PeptDeepInputBatch per run

    run() may be called concurrently (ONNXPredictorBase::forEachBatch()).
  */
  class ModelRuns
  {
  public:
    /// @throws Exception::IllegalArgument if @p model has an input that is not of @p supported
    ModelRuns(ONNXPredictorBase& model, const std::vector<InputKind>& supported) :
      model_(model),
      input_names_(model.getInputNames()),
      output_name_(model.session().GetOutputNameAllocated(0, allocator_))
    {
      static const std::vector<std::pair<std::string, InputKind>> names = {
        {"aa_indices", InputKind::AA_INDICES}, {"input_sequences", InputKind::AA_INDICES}, {"mod_x", InputKind::MOD_X},
        {"charges", InputKind::CHARGES}, {"charge", InputKind::CHARGES}, {"nces", InputKind::NCES}, {"nce", InputKind::NCES},
        {"instrument_indices", InputKind::INSTRUMENTS}, {"instrument", InputKind::INSTRUMENTS}};
      for (const std::string& name : input_names_)
      {
        const auto found = std::find_if(names.begin(), names.end(), [&name](const auto& n) { return n.first == name; });
        if (found == names.end() || std::find(supported.begin(), supported.end(), found->second) == supported.end())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
        }
        kinds_.push_back(found->second);
        input_name_ptrs_.push_back(name.c_str());
      }
    }

    /// Runs the model on @p batch (its buffers are the input tensors' memory); returns the output tensors
    std::vector<Ort::Value> run(PeptDeepInputBatch& batch) const
    {
      const int64_t n = static_cast<int64_t>(batch.batch_size);
      const int64_t length = static_cast<int64_t>(batch.sequence_length);
      const Ort::MemoryInfo& memory = model_.memoryInfo();
      std::vector<Ort::Value> inputs;
      inputs.reserve(kinds_.size());
      for (const InputKind kind : kinds_)
      {
        switch (kind)
        {
          case InputKind::AA_INDICES:
          {
            const int64_t shape[] = {n, length};
            inputs.push_back(Ort::Value::CreateTensor<int64_t>(memory, batch.aa_indices.data(), batch.aa_indices.size(), shape, 2));
            break;
          }
          case InputKind::MOD_X:
          {
            const int64_t shape[] = {n, length, static_cast<int64_t>(PEPTDEEP_MOD_ELEMENTS)};
            inputs.push_back(Ort::Value::CreateTensor<float>(memory, batch.mod_x.data(), batch.mod_x.size(), shape, 3));
            break;
          }
          case InputKind::CHARGES:
          {
            const int64_t shape[] = {n, 1};
            inputs.push_back(Ort::Value::CreateTensor<float>(memory, batch.charges.data(), batch.charges.size(), shape, 2));
            break;
          }
          case InputKind::NCES:
          {
            const int64_t shape[] = {n, 1};
            inputs.push_back(Ort::Value::CreateTensor<float>(memory, batch.nces.data(), batch.nces.size(), shape, 2));
            break;
          }
          case InputKind::INSTRUMENTS:
          {
            const int64_t shape[] = {n};
            inputs.push_back(Ort::Value::CreateTensor<int64_t>(memory, batch.instrument_indices.data(), batch.instrument_indices.size(), shape, 1));
            break;
          }
        }
      }
      const char* output_names[] = {output_name_.get()};
      return model_.session().Run(Ort::RunOptions{nullptr}, input_name_ptrs_.data(), inputs.data(), inputs.size(), output_names, 1);
    }

  private:
    ONNXPredictorBase& model_;
    std::vector<std::string> input_names_;
    std::vector<const char*> input_name_ptrs_;
    std::vector<InputKind> kinds_;
    Ort::AllocatorWithDefaultOptions allocator_;
    Ort::AllocatedStringPtr output_name_;
  };
} // namespace OpenMS::ML::Internal
