// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACNN.h"
#include <torch/nn/module.h>
#include <torch/nn/modules/conv.h>
#include <torch/nn/modules/linear.h>
#include <torch/nn/modules/activation.h>
#include <torch/nn/modules/pooling.h>
#include <torch/nn/modules/container/sequential.h>
#include <torch/nn/functional/loss.h>
#include <torch/optim/adam.h>
#include <torch/serialize/archive.h>
#include <ATen/Parallel.h>
#include <ATen/Context.h>
#include <torch/utils.h>
#include <torch/version.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr std::size_t WINDOW = 6 * 32 * 8;
    constexpr const char* VERSION = "OpenDIA-CNN-v1\nlog1p-per-window-max-library-sum-v2\nfnv1a64-seed-nul-group-mod10-v1\n6,32,8\n";
    // ponytail: serialize calls that change process-global LibTorch settings; use separate processes for concurrent models.
    std::mutex execution_mutex;
    using Indices = std::vector<std::size_t>;

    struct Model : torch::nn::Module
    {
      torch::nn::Sequential conv;
      torch::nn::Sequential head;
      Model()
      {
        conv = register_module("conv", torch::nn::Sequential(
          torch::nn::Conv2d(torch::nn::Conv2dOptions(6, 8, 3).padding(1)), torch::nn::ReLU(), torch::nn::AvgPool2d(2),
          torch::nn::Conv2d(torch::nn::Conv2dOptions(8, 12, 3).padding(1)), torch::nn::ReLU(), torch::nn::AvgPool2d(2)));
        head = register_module("head", torch::nn::Sequential(torch::nn::Linear(199, 24), torch::nn::ReLU(), torch::nn::Linear(24, 1)));
      }
      torch::Tensor forward(const torch::Tensor& raw, const torch::Tensor& aux)
      {
        return head->forward(torch::cat({conv->forward(raw).flatten(1), aux}, 1)).squeeze(1);
      }
    };

    void require(bool condition, const char* message)
    {
      if (!condition) { throw std::invalid_argument(message); }
    }
    bool validID(const std::string& id)
    {
      return !id.empty() && id.find_first_of("\t\r\n") == std::string::npos && id.find('\0') == std::string::npos &&
        id.find_first_not_of(' ') != std::string::npos;
    }
    void validate(const CNNData& data, const CNNGeometry& geometry)
    {
      for (double v : {geometry.mz_tolerance, geometry.rt_half_width, geometry.im_half_width})
      {
        require(std::isfinite(v) && v > 0, "geometry must be positive and finite");
      }
      require(!data.samples.empty() && data.samples.size() <= std::numeric_limits<std::size_t>::max() / WINDOW,
              "samples must be nonempty with representable shape");
      require(data.intensities && data.intensity_count == data.samples.size() * WINDOW, "raw shape must be [N,6,32,8]");
      std::set<std::string> seen;
      for (std::size_t i = 0; i < data.samples.size(); ++i)
      {
        const auto& s = data.samples[i];
        require(validID(s.sample_id) && validID(s.group_id) && seen.insert(s.sample_id).second, "invalid or duplicate sample ID/group ID");
        require(s.label >= -1 && s.label <= 1 && !(s.is_decoy && s.label == 1), "invalid label or positive decoy");
        for (std::size_t c = 0; c < 6; ++c)
        {
          const double lib = s.library_intensities[c], mz = s.fragment_mz[c];
          require(std::isfinite(lib) && lib >= 0 && std::isfinite(mz) && mz >= 0, "invalid library or fragment m/z");
          require(c == 0 || s.library_intensities[c - 1] >= lib, "library intensities must be sorted descending");
          require(mz != 0 || lib == 0, "padded channels require zero library intensity");
          for (std::size_t j = 0; j < 256; ++j)
          {
            const float value = data.intensities[i * WINDOW + c * 256 + j];
            require(std::isfinite(value) && value >= 0 && (mz != 0 || value == 0), "invalid raw or nonzero padded intensity");
          }
        }
      }
    }
    std::uint64_t hashGroup(const std::string& group, int seed)
    {
      std::string key = std::to_string(seed);
      key.push_back('\0');
      key += group;
      std::uint64_t hash = UINT64_C(14695981039346656037);
      for (unsigned char c : key) { hash = (hash ^ c) * UINT64_C(1099511628211); }
      return hash;
    }
    std::string split(const std::string& group, int seed)
    {
      const auto bucket = hashGroup(group, seed) % 10;
      return bucket < 6 ? "train" : bucket < 8 ? "validation" : "heldout";
    }
    std::string partition(const std::string& group, const CNNOptions& options)
    {
      if (options.outer_folds == 0) { return split(group, options.seed); }
      if (traceCNNFold(group, options.seed, options.outer_folds) == options.heldout_fold) { return "heldout"; }
      return split(group, options.seed) == "train" ? "train" : "validation";
    }
    struct Batch
    {
      torch::Tensor raw;
      torch::Tensor aux;
      torch::Tensor labels;
    };
    Batch batch(const CNNData& data, const Indices& ids, std::size_t begin, std::size_t end)
    {
      const auto n = static_cast<int64_t>(end - begin);
      Batch b{torch::empty({n, 6, 32, 8}, torch::kFloat32), torch::empty({n, 7}, torch::kFloat32), torch::empty({n}, torch::kFloat32)};
      float* x = b.raw.data_ptr<float>();
      float* a = b.aux.data_ptr<float>();
      float* y = b.labels.data_ptr<float>();
      for (std::size_t k = begin; k < end; ++k)
      {
        const auto index = ids[k], row = k - begin;
        const auto& s = data.samples[index];
        float maximum = 1e-12f;
        for (std::size_t j = 0; j < WINDOW; ++j)
        {
          const float value = std::log1p(data.intensities[index * WINDOW + j]);
          x[row * WINDOW + j] = value;
          maximum = std::max(maximum, value);
        }
        for (std::size_t j = 0; j < WINDOW; ++j) { x[row * WINDOW + j] /= maximum; }
        // Scale in double first to avoid overflow for finite double library inputs.
        const double lib_max = std::max(s.library_intensities[0], 1e-12);
        float sum = 0;
        for (std::size_t c = 0; c < 6; ++c)
        {
          a[row * 7 + c] = static_cast<float>(s.library_intensities[c] / lib_max);
          sum += a[row * 7 + c];
        }
        for (std::size_t c = 0; c < 6; ++c) { a[row * 7 + c] /= std::max(sum, 1e-12f); }
        a[row * 7 + 6] = s.has_im ? 1.f : 0.f;
        y[row] = static_cast<float>(s.label);
      }
      return b;
    }
    torch::Tensor loss(Model& model, const Batch& b)
    {
      return torch::nn::functional::binary_cross_entropy_with_logits(model.forward(b.raw, b.aux), b.labels);
    }
    void predict(Model& model, const CNNData& data, CNNResult& result)
    {
      torch::NoGradGuard guard;
      model.eval();
      Indices ids(data.samples.size());
      std::iota(ids.begin(), ids.end(), 0);
      result.logits.clear();
      for (std::size_t start = 0; start < ids.size(); start += 128)
      {
        auto b = batch(data, ids, start, std::min(start + 128, ids.size()));
        auto out = model.forward(b.raw, b.aux).contiguous();
        const auto* p = out.data_ptr<float>();
        result.logits.insert(result.logits.end(), p, p + out.numel());
      }
      std::vector<float> negatives;
      std::vector<float> positives;
      for (std::size_t i = 0; i < ids.size(); ++i)
      {
        require(std::isfinite(result.logits[i]), "model produced nonfinite logits");
        if (result.splits[i] != "heldout") { continue; }
        if (data.samples[i].label == 0) { negatives.push_back(result.logits[i]); }
        if (data.samples[i].label == 1) { positives.push_back(result.logits[i]); }
      }
      if (!negatives.empty() && !positives.empty())
      {
        std::sort(negatives.begin(), negatives.end());
        double wins = 0;
        for (float p : positives)
        {
          wins += 0.5 * ((std::lower_bound(negatives.begin(), negatives.end(), p) - negatives.begin()) +
                        (std::upper_bound(negatives.begin(), negatives.end(), p) - negatives.begin()));
        }
        result.heldout_auc = wins / (static_cast<double>(positives.size()) * negatives.size());
      }
    }
    torch::Tensor textTensor(const std::string& text)
    {
      return torch::from_blob(const_cast<char*>(text.data()), {static_cast<int64_t>(text.size())}, torch::kUInt8).clone();
    }
    std::string tensorText(const torch::Tensor& tensor)
    {
      require(tensor.device().is_cpu() && tensor.scalar_type() == torch::kUInt8 && tensor.dim() == 1 && tensor.is_contiguous(), "invalid checkpoint metadata");
      return std::string(reinterpret_cast<const char*>(tensor.data_ptr<uint8_t>()), tensor.numel());
    }
  }

  int traceCNNFold(const std::string& group_id, int seed, int folds)
  {
    require(validID(group_id) && seed >= 0 && folds >= 2 && folds <= 10, "invalid outer family fold contract");
    // Domain separation avoids reusing the inner train/validation hash buckets.
    return static_cast<int>(hashGroup("outer\t" + group_id, seed) % static_cast<std::uint64_t>(folds));
  }

  CNNResult trainTraceCNN(const CNNData& data, const CNNGeometry& geometry, const CNNOptions& options, const std::string& checkpoint)
  {
    std::lock_guard<std::mutex> lock(execution_mutex);
    validate(data, geometry);
    require(options.seed >= 0 && options.threads > 0 && options.epochs >= 1 && options.epochs <= 100 && options.patience >= 1 && options.patience <= 20 && options.batch_size > 0 &&
            ((options.outer_folds == 0 && options.heldout_fold == -1) ||
             (options.outer_folds >= 2 && options.outer_folds <= 10 && options.heldout_fold >= 0 && options.heldout_fold < options.outer_folds)),
            "invalid threads, epochs, patience or batch size");
    at::set_num_threads(options.threads);
    torch::manual_seed(options.seed);
    at::globalContext().setDeterministicAlgorithms(true, false);
    CNNResult result;
    Indices train, validation;
    std::map<std::string, std::string> provenance;
    bool positive = false, negative = false;
    for (std::size_t i = 0; i < data.samples.size(); ++i)
    {
      const auto& s = data.samples[i];
      const auto role = partition(s.group_id, options);
      provenance[s.group_id] = role;
      result.splits.push_back(role);
      if (s.label < 0) { continue; }
      if (role == "train")
      {
        train.push_back(i);
        positive |= s.label == 1;
        negative |= s.label == 0;
      }
      if (role == "validation") { validation.push_back(i); }
    }
    require(positive && negative && !validation.empty(), "training requires both classes and labeled validation families");
    Model model;
    torch::optim::Adam optimizer(model.parameters(), torch::optim::AdamOptions(0.003));
    double best = std::numeric_limits<double>::infinity();
    int stale = 0;
    std::string best_state;
    for (int epoch = 0; epoch < options.epochs; ++epoch)
    {
      model.train();
      auto permutation = torch::randperm(static_cast<int64_t>(train.size()), torch::kInt64);
      Indices order;
      for (int64_t j = 0; j < permutation.numel(); ++j) { order.push_back(train[permutation.data_ptr<int64_t>()[j]]); }
      for (std::size_t start = 0; start < order.size(); start += options.batch_size)
      {
        auto b = batch(data, order, start, std::min(start + options.batch_size, order.size()));
        optimizer.zero_grad();
        auto value = loss(model, b);
        require(std::isfinite(value.item<double>()), "nonfinite training loss");
        value.backward();
        optimizer.step();
      }
      model.eval();
      double value = 0;
      {
        torch::NoGradGuard guard;
        for (std::size_t start = 0; start < validation.size(); start += options.batch_size)
        {
          const auto end = std::min(start + options.batch_size, validation.size());
          value += loss(model, batch(data, validation, start, end)).item<double>() * (end - start);
        }
      }
      value /= validation.size();
      require(std::isfinite(value), "nonfinite validation loss");
      result.epochs_run = epoch + 1;
      if (value < best - 1e-6)
      {
        best = value;
        stale = 0;
        torch::serialize::OutputArchive snapshot;
        model.save(snapshot);
        std::ostringstream stream(std::ios::binary);
        snapshot.save_to(stream);
        best_state = stream.str();
      }
      else { ++stale; }
      if (stale >= options.patience) { break; }
    }
    // Archives own copied bytes, unlike a shallow Tensor/state_dict snapshot.
    std::istringstream restored(best_state, std::ios::binary);
    torch::serialize::InputArchive snapshot;
    snapshot.load_from(restored, torch::kCPU);
    model.load(snapshot);
    result.validation_loss = best;
    predict(model, data, result);
    torch::serialize::OutputArchive archive, model_archive;
    model.save(model_archive);
    archive.write("model", model_archive);
    std::string metadata = VERSION;
    for (const auto& [group, partition] : provenance) { metadata += partition + "\t" + group + "\n"; }
    archive.write("metadata", textTensor(metadata), true);
    archive.write("geometry", torch::tensor({geometry.mz_tolerance, geometry.rt_half_width, geometry.im_half_width}, torch::kFloat64), true);
    archive.write("training", torch::tensor({static_cast<double>(result.epochs_run), result.validation_loss}, torch::kFloat64), true);
    archive.write("training_options", torch::tensor({static_cast<int64_t>(options.seed), static_cast<int64_t>(options.threads),
                  static_cast<int64_t>(options.epochs), static_cast<int64_t>(options.patience), static_cast<int64_t>(options.batch_size)}, torch::kInt64), true);
    archive.write("torch_version", textTensor(TORCH_VERSION), true);
    archive.write("outer_partition", torch::tensor({static_cast<int64_t>(options.outer_folds),
                  static_cast<int64_t>(options.heldout_fold)}, torch::kInt64), true);
    archive.save_to(checkpoint);
    return result;
  }

  CNNResult scoreTraceCNN(const CNNData& data, const CNNGeometry& geometry, const std::string& checkpoint, int threads, bool allow_training_overlap)
  {
    std::lock_guard<std::mutex> lock(execution_mutex);
    validate(data, geometry);
    require(threads > 0, "threads must be positive");
    at::set_num_threads(threads);
    torch::serialize::InputArchive archive, model_archive;
    archive.load_from(checkpoint, torch::kCPU);
    torch::Tensor metadata, saved_geometry, training, training_options, torch_version;
    archive.read("metadata", metadata, true);
    archive.read("geometry", saved_geometry, true);
    archive.read("training", training, true);
    archive.read("training_options", training_options, true);
    archive.read("torch_version", torch_version, true);
    require(tensorText(torch_version) == TORCH_VERSION, "checkpoint LibTorch version mismatch");
    require(training_options.scalar_type() == torch::kInt64 && training_options.dim() == 1 && training_options.numel() == 5,
            "invalid checkpoint training options");
    const auto seed = training_options[0].item<int64_t>();
    const auto saved_threads = training_options[1].item<int64_t>();
    const auto epoch_cap = training_options[2].item<int64_t>();
    const auto patience = training_options[3].item<int64_t>();
    const auto batch_size = training_options[4].item<int64_t>();
    require(seed >= std::numeric_limits<int>::min() && seed <= std::numeric_limits<int>::max() &&
            saved_threads > 0 && saved_threads <= std::numeric_limits<int>::max() && epoch_cap >= 1 && epoch_cap <= 100 &&
            patience >= 1 && patience <= 20 && batch_size > 0 && batch_size <= std::numeric_limits<int>::max(),
            "invalid checkpoint training option values");
    require(saved_geometry.scalar_type() == torch::kFloat64 && saved_geometry.dim() == 1 && saved_geometry.numel() == 3, "invalid checkpoint geometry");
    require(saved_geometry[0].item<double>() == geometry.mz_tolerance && saved_geometry[1].item<double>() == geometry.rt_half_width &&
            saved_geometry[2].item<double>() == geometry.im_half_width, "checkpoint sampling geometry mismatch");
    CNNOptions saved_options;
    saved_options.seed = static_cast<int>(seed);
    torch::Tensor outer_partition;
    if (archive.try_read("outer_partition", outer_partition, true))
    {
      require(outer_partition.scalar_type() == torch::kInt64 && outer_partition.dim() == 1 && outer_partition.numel() == 2,
              "invalid checkpoint outer partition");
      const auto folds = outer_partition[0].item<int64_t>();
      const auto heldout_fold = outer_partition[1].item<int64_t>();
      require((folds == 0 && heldout_fold == -1) || (folds >= 2 && folds <= 10 && heldout_fold >= 0 && heldout_fold < folds),
              "invalid checkpoint outer fold values");
      saved_options.outer_folds = static_cast<int>(folds);
      saved_options.heldout_fold = static_cast<int>(heldout_fold);
    }
    const auto text = tensorText(metadata);
    const std::string version = VERSION;
    require(text.compare(0, version.size(), version) == 0, "checkpoint version/preprocessing/split mismatch");
    std::map<std::string, std::string> provenance;
    std::istringstream lines(text.substr(version.size()));
    std::string line;
    while (std::getline(lines, line))
    {
      const auto tab = line.find('\t');
      require(tab != std::string::npos, "invalid checkpoint provenance");
      const auto role = line.substr(0, tab), group = line.substr(tab + 1);
      require((role == "train" || role == "validation" || role == "heldout") && validID(group) &&
              role == partition(group, saved_options) && provenance.emplace(group, role).second,
              "invalid checkpoint family provenance");
    }
    require(!provenance.empty(), "checkpoint provenance is empty");
    CNNResult result;
    for (const auto& sample : data.samples)
    {
      const auto found = provenance.find(sample.group_id);
      const std::string partition = found == provenance.end() ? "external" : found->second;
      require(allow_training_overlap || partition != "train", "scoring overlaps training families");
      result.splits.push_back(partition);
    }
    require(training.scalar_type() == torch::kFloat64 && training.dim() == 1 && training.numel() == 2, "invalid checkpoint training metadata");
    const double epochs = training[0].item<double>();
    result.validation_loss = training[1].item<double>();
    require(std::isfinite(epochs) && epochs >= 1 && epochs <= epoch_cap && std::floor(epochs) == epochs &&
            std::isfinite(result.validation_loss) && result.validation_loss >= 0, "invalid checkpoint training values");
    result.epochs_run = static_cast<int>(epochs);
    Model model;
    archive.read("model", model_archive);
    model.load(model_archive);
    predict(model, data, result);
    return result;
  }
}
