// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepInput.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>

#include <onnxruntime_cxx_api.h>

#include <map>
#include <string>
#include <utility>

namespace OpenMS
{
  PeptDeepCCSInference::PeptDeepCCSInference(const std::string& model_path)
    : model_(model_path)
  {}

  PeptDeepCCSInference::~PeptDeepCCSInference() = default;

  std::vector<float> PeptDeepCCSInference::predictCCS(
    const std::vector<std::string>& peptides,
    const std::vector<float>& charges)
  {
    if (peptides.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide batch cannot be empty.");
    }
    if (charges.size() != peptides.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide and charge input vectors must have the same size.");
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
    const std::vector<std::string> input_names = model_.getInputNames();
    std::vector<const char*> input_names_chars;
    input_names_chars.reserve(input_names.size());
    for (const std::string& name : input_names)
    {
      input_names_chars.push_back(name.c_str());
    }

    Ort::AllocatorWithDefaultOptions ort_alloc;
    Ort::AllocatedStringPtr output_name_ptr = model_.session().GetOutputNameAllocated(0, ort_alloc);
    const char* output_names[] = {output_name_ptr.get()};

    for (const std::vector<size_t>& group_indices : groups)
    {
      std::vector<std::string> group_peptides;
      std::vector<float> group_charges;
      group_peptides.reserve(group_indices.size());
      group_charges.reserve(group_indices.size());
      for (size_t idx : group_indices)
      {
        group_peptides.push_back(peptides[idx]);
        group_charges.push_back(charges[idx]);
      }

      ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildUnmodifiedChargedBatch(group_peptides, group_charges, input_config);
      const int64_t batch_size = static_cast<int64_t>(batch.batch_size);
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
          expected_shape[0] = batch_size;
        }

        if (name == "aa_indices" || name == "input_sequences")
        {
          expected_shape = {batch_size, sequence_length};
          input_tensors.push_back(Ort::Value::CreateTensor<int64_t>(
            model_.memoryInfo(),
            batch.aa_indices.data(),
            batch.aa_indices.size(),
            expected_shape.data(),
            expected_shape.size()));
        }
        else if (name == "mod_x")
        {
          expected_shape = {batch_size, sequence_length, ML::PEPTDEEP_MOD_ELEMENTS};
          input_tensors.push_back(Ort::Value::CreateTensor<float>(
            model_.memoryInfo(),
            batch.mod_x.data(),
            batch.mod_x.size(),
            expected_shape.data(),
            expected_shape.size()));
        }
        else if (name == "charges" || name == "charge")
        {
          expected_shape = {batch_size, 1};
          input_tensors.push_back(Ort::Value::CreateTensor<float>(
            model_.memoryInfo(),
            batch.charges.data(),
            batch.charges.size(),
            expected_shape.data(),
            expected_shape.size()));
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ONNX input: " + name);
        }
      }

      auto output_tensors = model_.session().Run(
        Ort::RunOptions{nullptr},
        input_names_chars.data(),
        input_tensors.data(),
        input_tensors.size(),
        output_names,
        1);

      const size_t output_count = output_tensors.front().GetTensorTypeAndShapeInfo().GetElementCount();
      if (output_count < batch.batch_size)
      {
        throw Exception::InvalidValue(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "ONNX model output shape mismatch.",
          std::to_string(output_count));
      }

      float* output = output_tensors.front().GetTensorMutableData<float>();
      for (size_t i = 0; i < group_indices.size(); ++i)
      {
        predictions[group_indices[i]] = output[i];
      }
    }

    return predictions;
  }
} // namespace OpenMS
