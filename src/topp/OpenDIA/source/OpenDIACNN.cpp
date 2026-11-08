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
#include <torch/nn/modules/dropout.h>
#include <torch/nn/functional/loss.h>
#include <torch/optim/adam.h>
#include <torch/optim/adamw.h>
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

    // ---------------------------------------------------------------- TOPAZ-style architecture
    // Follows redeem-topaz (model/topaz.rs, building_blocks/{conv_encoder,coelution,trace_input}.rs).
    constexpr int64_t RT_BINS = 32, IM_BINS = 8, FRAGMENTS = 6, ISOTOPES = 3, CORE_BEGIN = 12, CORE_END = 20;
    constexpr std::size_t MS1_WINDOW = ISOTOPES * RT_BINS;
    constexpr const char* TOPAZ_VERSION = "OpenDIA-CNN-topaz-v1\nxic-sum-im_ms1-m012_xim-core-rt12-19_dual-maxnorm-log1p\n"
                                          "fnv1a64-seed-nul-group-mod10-v1\n6,32,8;3,32\n";
    constexpr double EPS = 1e-6, APEX_BETA = 10.0;
    constexpr int64_t MAX_LAG = 2;

    /// Two views per channel, max-normalised shape and log1p intensity: [n,C,L] -> [n,2C,L].
    torch::Tensor dualView(const torch::Tensor& x)
    {
      return torch::cat({x / (x.amax(2, true) + EPS), torch::log1p(x.clamp_min(0))}, 1);
    }
    torch::Tensor maxNormalise(const torch::Tensor& x) { return x / (x.amax(2, true) + EPS); }
    /// Z-score along the trace axis (unbiased spread, as candle's var).
    torch::Tensor zscore(const torch::Tensor& x)
    {
      const auto centred = x - x.mean(2, true);
      return centred / (torch::sqrt(centred.pow(2).sum(2, true) / static_cast<double>(x.size(2) - 1)) + EPS);
    }
    torch::Tensor unitLength(const torch::Tensor& z) { return z / torch::sqrt(z.pow(2).sum(2, true) + EPS); }
    /// Soft apex position in [0,1] per channel: [n,C,L] -> [n,C].
    torch::Tensor softApex(const torch::Tensor& normalised)
    {
      const auto w = torch::softmax(normalised * APEX_BETA, 2);
      return (w * torch::linspace(0.0, 1.0, normalised.size(2)).view({1, 1, -1})).sum(2);
    }
    torch::Tensor unbiasedStd(const torch::Tensor& x)
    {
      return torch::sqrt((x - x.mean(1, true)).pow(2).sum(1) / static_cast<double>(x.size(1) - 1) + 1e-12);
    }

    /// One 1D convolutional branch over a trace axis (TOPAZ ConvBranch): [n,C,L] -> [n,2w]. TOPAZ uses w = 32.
    struct BranchImpl : torch::nn::Module
    {
      torch::nn::Conv1d conv0{nullptr}, conv1{nullptr}, conv2{nullptr};
      torch::nn::Linear proj{nullptr};
      BranchImpl(int64_t channels, int64_t width)
      {
        conv0 = register_module("conv0", torch::nn::Conv1d(torch::nn::Conv1dOptions(2 * channels, width, 1)));
        conv1 = register_module("conv1", torch::nn::Conv1d(torch::nn::Conv1dOptions(width, 2 * width, 3).padding(1)));
        conv2 = register_module("conv2", torch::nn::Conv1d(torch::nn::Conv1dOptions(2 * width, 2 * width, 3).padding(1)));
        proj = register_module("proj", torch::nn::Linear(4 * width, 2 * width));
      }
      torch::Tensor forward(const torch::Tensor& x)
      {
        auto h = torch::relu(conv0(dualView(x)));
        h = torch::relu(conv1(h));
        h = torch::relu(conv2(h));
        return proj(torch::cat({h.amax(2), h.mean(2)}, 1));
      }
    };
    TORCH_MODULE(Branch);

    /// Explicit co-elution statistics plus a learned pairwise-similarity embedding (TOPAZ CoelutionHead): [n,C,L] -> [n,18].
    struct CoelutionImpl : torch::nn::Module
    {
      torch::Tensor upper;
      torch::nn::Sequential similarity{nullptr};
      explicit CoelutionImpl(int64_t channels)
      {
        const auto pairs = torch::triu_indices(channels, channels, 1);
        upper = register_buffer("upper", pairs[0] * channels + pairs[1]);
        const int64_t count = channels * (channels - 1) / 2;
        similarity = register_module("similarity", torch::nn::Sequential(torch::nn::Linear(2 * count, 32), torch::nn::ReLU(), torch::nn::Linear(32, 8)));
      }
      /// Upper-triangle cosines between the channels of a and of b.
      torch::Tensor pairCosines(const torch::Tensor& a, const torch::Tensor& b) const
      {
        return torch::bmm(unitLength(a), unitLength(b).transpose(1, 2)).flatten(1).index_select(1, upper);
      }
      torch::Tensor forward(const torch::Tensor& x)
      {
        const int64_t length = x.size(2);
        const auto normalised = maxNormalise(x);
        const auto z = zscore(normalised);
        const auto same = pairCosines(z, z);
        // Lag-tolerant cosine: best over shifts of up to MAX_LAG bins.
        std::vector<torch::Tensor> lagged;
        for (int64_t lag = -MAX_LAG; lag <= MAX_LAG; ++lag)
        {
          const int64_t k = std::abs(lag);
          if (k == 0) { lagged.push_back(same); continue; }
          lagged.push_back(lag > 0 ? pairCosines(z.narrow(2, k, length - k), z.narrow(2, 0, length - k))
                                   : pairCosines(z.narrow(2, 0, length - k), z.narrow(2, k, length - k)));
        }
        const auto lag_best = torch::stack(lagged, 2).amax(2);
        const auto apex = softApex(normalised);
        const auto apex_mean = apex.mean(1, true);
        const auto weights = torch::softmax(normalised * APEX_BETA, 2);
        const auto entropy = -(weights * torch::log(weights + EPS)).sum(2);
        const auto peakiness = torch::sigmoid((normalised.amax(2) / (normalised.mean(2) + EPS) - 2.0) * 2.0).mean(1);
        const auto summary = torch::stack({same.mean(1), same.amin(1), lag_best.mean(1), lag_best.amin(1), unbiasedStd(apex),
                                           apex.amax(1) - apex.amin(1), (apex - apex_mean).abs().mean(1), entropy.mean(1), entropy.amin(1), peakiness}, 1);
        return torch::cat({summary, similarity->forward(torch::cat({same, lag_best}, 1))}, 1);
      }
    };
    TORCH_MODULE(Coelution);

    /// MS1-vs-MS2 agreement (TOPAZ ms12_features): mean/max cosine, apex mean difference, apex spread: [n,4].
    torch::Tensor crossModal(const torch::Tensor& ms1, const torch::Tensor& ms2)
    {
      const auto n1 = maxNormalise(ms1), n2 = maxNormalise(ms2);
      const auto cosines = torch::bmm(unitLength(zscore(n1)), unitLength(zscore(n2)).transpose(1, 2));
      const auto apex1 = softApex(n1), apex2 = softApex(n2);
      return torch::stack({cosines.mean({1, 2}), cosines.flatten(1).amax(1), (apex1.mean(1) - apex2.mean(1)).abs(),
                           unbiasedStd(torch::cat({apex1, apex2}, 1))}, 1);
    }
    torch::Tensor rowCosine(const torch::Tensor& a, const torch::Tensor& b)
    {
      return (a * b).sum(1) / (torch::sqrt(a.pow(2).sum(1) + 1e-12) * torch::sqrt(b.pow(2).sum(1) + 1e-12));
    }

    /// Raw-intensity inputs: fragment XIC [n,6,32] (IM summed), MS1 M+0..M+2 [n,3,32], fragment mobilogram over the
    /// core RT bins [n,6,8], aux [n,8] = library intensities / max, has_ms1, has_im.
    struct TopazBatch
    {
      torch::Tensor xic, ms1, xim, aux, labels;
    };

    struct TopazModel : torch::nn::Module
    {
      int64_t embedding;
      Branch xic{nullptr}, ms1{nullptr}, xim{nullptr};
      Coelution xic_coelution{nullptr}, ms1_coelution{nullptr}, xim_coelution{nullptr};
      torch::nn::Sequential head{nullptr};
      /// width 32 is the TOPAZ layout (conv 32-64-64, embeddings 64, MLP 128-64); smaller widths scale every learned layer.
      explicit TopazModel(int64_t width) : embedding(2 * width)
      {
        xic = register_module("xic", Branch(FRAGMENTS, width));
        ms1 = register_module("ms1", Branch(ISOTOPES, width));
        xim = register_module("xim", Branch(FRAGMENTS, width));
        xic_coelution = register_module("xic_coelution", Coelution(FRAGMENTS));
        ms1_coelution = register_module("ms1_coelution", Coelution(ISOTOPES));
        xim_coelution = register_module("xim_coelution", Coelution(FRAGMENTS));
        const int64_t features = embedding + 18 + (embedding + 18 + 4) + (embedding + 18) + 2 + 8;
        head = register_module("head", torch::nn::Sequential(torch::nn::Linear(features, 4 * width), torch::nn::ReLU(), torch::nn::Dropout(0.1),
                                                             torch::nn::Linear(4 * width, 2 * width), torch::nn::ReLU(), torch::nn::Dropout(0.1),
                                                             torch::nn::Linear(2 * width, 1)));
      }
      torch::Tensor forward(const TopazBatch& b)
      {
        const auto library = b.aux.narrow(1, 0, FRAGMENTS);
        const auto has_ms1 = b.aux.narrow(1, FRAGMENTS, 1), has_im = b.aux.narrow(1, FRAGMENTS + 1, 1);
        // Absent modalities contribute exact zeros, not learned biases on zero-filled input. A batch without any
        // sample of a modality skips its branch (same zeros, no compute).
        const auto n = b.xic.size(0);
        const auto ms1_part = has_ms1.max().item<float>() == 0.f ? torch::zeros({n, embedding + 18 + 4})
          : torch::cat({ms1(b.ms1), ms1_coelution(b.ms1), crossModal(b.ms1, b.xic)}, 1) * has_ms1;
        const auto xim_part = has_im.max().item<float>() == 0.f ? torch::zeros({n, embedding + 18})
          : torch::cat({xim(b.xim), xim_coelution(b.xim)}, 1) * has_im;
        const auto core = b.xic.narrow(2, CORE_BEGIN, CORE_END - CORE_BEGIN).sum(2);
        const auto library_fit = torch::stack({rowCosine(torch::sqrt(core), torch::sqrt(library)), rowCosine(core, library)}, 1);
        return head->forward(torch::cat({xic(b.xic), xic_coelution(b.xic), ms1_part, xim_part, library_fit, b.aux}, 1)).squeeze(1);
      }
    };

    TopazBatch topazBatch(const CNNData& data, const std::vector<std::size_t>& ids, std::size_t begin, std::size_t end)
    {
      const auto n = static_cast<int64_t>(end - begin);
      TopazBatch b{torch::zeros({n, FRAGMENTS, RT_BINS}), torch::zeros({n, ISOTOPES, RT_BINS}), torch::zeros({n, FRAGMENTS, IM_BINS}),
                   torch::empty({n, FRAGMENTS + 2}), torch::empty({n})};
      float* xic = b.xic.data_ptr<float>();
      float* ms1 = b.ms1.data_ptr<float>();
      float* xim = b.xim.data_ptr<float>();
      float* aux = b.aux.data_ptr<float>();
      float* y = b.labels.data_ptr<float>();
      for (std::size_t k = begin; k < end; ++k)
      {
        const auto index = ids[k], row = k - begin;
        const auto& s = data.samples[index];
        const float* grid = data.intensities + index * static_cast<std::size_t>(FRAGMENTS * RT_BINS * IM_BINS);
        for (int64_t c = 0; c < FRAGMENTS; ++c)
        {
          for (int64_t rt = 0; rt < RT_BINS; ++rt)
          {
            for (int64_t im = 0; im < IM_BINS; ++im)
            {
              const float value = grid[(c * RT_BINS + rt) * IM_BINS + im];
              xic[(row * FRAGMENTS + c) * RT_BINS + rt] += value;
              if (rt >= CORE_BEGIN && rt < CORE_END) { xim[(row * FRAGMENTS + c) * IM_BINS + im] += value; }
            }
          }
        }
        if (data.ms1 != nullptr) { std::copy_n(data.ms1 + index * MS1_WINDOW, MS1_WINDOW, ms1 + row * MS1_WINDOW); }
        const double lib_max = std::max(s.library_intensities[0], 1e-12);
        for (int64_t c = 0; c < FRAGMENTS; ++c) { aux[row * (FRAGMENTS + 2) + c] = static_cast<float>(s.library_intensities[c] / lib_max); }
        aux[row * (FRAGMENTS + 2) + FRAGMENTS] = s.has_ms1 ? 1.f : 0.f;
        aux[row * (FRAGMENTS + 2) + FRAGMENTS + 1] = s.has_im ? 1.f : 0.f;
        y[row] = static_cast<float>(s.label);
      }
      return b;
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

    void validateMS1(const CNNData& data)
    {
      require(data.ms1 == nullptr ? data.ms1_count == 0 : data.ms1_count == data.samples.size() * MS1_WINDOW, "MS1 traces must be [N,3,32]");
      for (std::size_t i = 0; i < data.samples.size(); ++i)
      {
        require(!data.samples[i].has_ms1 || data.ms1 != nullptr, "has_ms1 requires MS1 traces");
      }
      for (std::size_t j = 0; j < data.ms1_count; ++j)
      {
        require(std::isfinite(data.ms1[j]) && data.ms1[j] >= 0, "invalid MS1 trace intensity");
      }
    }
    void predictTopaz(TopazModel& model, const CNNData& data, CNNResult& result, bool heldout_only = false)
    {
      torch::NoGradGuard guard;
      model.eval();
      Indices ids;
      for (std::size_t i = 0; i < data.samples.size(); ++i)
      {
        if (!heldout_only || result.splits[i] == "heldout") { ids.push_back(i); }
      }
      result.logits.assign(data.samples.size(), std::numeric_limits<float>::quiet_NaN());
      for (std::size_t start = 0; start < ids.size(); start += 512)
      {
        const auto end = std::min(start + 512, ids.size());
        auto out = model.forward(topazBatch(data, ids, start, end)).contiguous();
        const auto* p = out.data_ptr<float>();
        for (std::size_t k = start; k < end; ++k) { result.logits[ids[k]] = p[k - start]; }
      }
      std::vector<float> negatives, positives;
      for (std::size_t i : ids)
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

    /// Candidate bags for multiple-instance training: windows sharing a nonempty bag_id, label = target.
    struct Bag
    {
      Indices members;
      float label = 0;
    };
    /// Masked-max bag logits: each bag scores as its best candidate (TOPAZ forward_bags).
    torch::Tensor bagLogits(TopazModel& model, const CNNData& data, const std::vector<Bag>& bags, const Indices& order, std::size_t begin,
                            std::size_t end, torch::Tensor& labels)
    {
      Indices windows;
      std::vector<int64_t> owner;
      labels = torch::empty({static_cast<int64_t>(end - begin)});
      for (std::size_t k = begin; k < end; ++k)
      {
        const auto& bag = bags[order[k]];
        windows.insert(windows.end(), bag.members.begin(), bag.members.end());
        owner.insert(owner.end(), bag.members.size(), static_cast<int64_t>(k - begin));
        labels[static_cast<int64_t>(k - begin)] = bag.label;
      }
      const auto logits = model.forward(topazBatch(data, windows, 0, windows.size()));
      const auto index = torch::tensor(owner, torch::kInt64);
      return torch::full({static_cast<int64_t>(end - begin)}, -1e9f).scatter_reduce(0, index, logits, "amax", true);
    }

    CNNResult trainTopaz(const CNNData& data, const CNNGeometry& geometry, const CNNOptions& options, const std::string& checkpoint)
    {
      validateMS1(data);
      CNNResult result;
      Indices train, validation;
      std::map<std::string, std::string> provenance;
      bool positive = false, negative = false;
      std::map<std::string, Bag> train_bags, validation_bags;
      std::map<std::string, std::pair<std::string, bool>> bag_identity;
      for (std::size_t i = 0; i < data.samples.size(); ++i)
      {
        const auto& s = data.samples[i];
        const auto role = partition(s.group_id, options);
        provenance[s.group_id] = role;
        result.splits.push_back(role);
        if (options.mil)
        {
          if (s.bag_id.empty()) { continue; }
          const auto [it, inserted] = bag_identity.emplace(s.bag_id, std::make_pair(s.group_id, s.is_decoy));
          require(inserted || it->second == std::make_pair(s.group_id, s.is_decoy), "a MIL bag must have one family and one decoy status");
          if (role == "heldout") { continue; }
          auto& bag = (role == "train" ? train_bags : validation_bags)[s.bag_id];
          bag.members.push_back(i);
          bag.label = s.is_decoy ? 0.f : 1.f;
          continue;
        }
        if (s.label < 0) { continue; }
        if (role == "train")
        {
          train.push_back(i);
          positive |= s.label == 1;
          negative |= s.label == 0;
        }
        if (role == "validation") { validation.push_back(i); }
      }
      // Ordered (std::map) bag lists keep MIL sampling independent of input order.
      std::vector<Bag> fit_bags, check_bags;
      for (auto& [id, bag] : train_bags) { fit_bags.push_back(std::move(bag)); }
      for (auto& [id, bag] : validation_bags) { check_bags.push_back(std::move(bag)); }
      for (const auto& bag : fit_bags) { positive |= bag.label == 1; negative |= bag.label == 0; }
      require(positive && negative && (options.mil ? !check_bags.empty() : !validation.empty()),
              "training requires both classes and labeled validation families");
      // Fixed validation subsample (MIL): the first bags of one seeded permutation.
      Indices check_order(check_bags.size());
      std::iota(check_order.begin(), check_order.end(), 0);
      TopazModel model(options.topaz_width);
      // Weight decay (AdamW) is the compact-variant regulariser; zero keeps plain Adam behaviour.
      torch::optim::AdamW optimizer(model.parameters(), torch::optim::AdamWOptions(1e-3).weight_decay(options.weight_decay));
      if (options.mil && check_order.size() > static_cast<std::size_t>(options.mil_bags_per_epoch))
      {
        auto permutation = torch::randperm(static_cast<int64_t>(check_order.size()), torch::kInt64);
        Indices chosen;
        for (int64_t j = 0; j < options.mil_bags_per_epoch; ++j) { chosen.push_back(static_cast<std::size_t>(permutation.data_ptr<int64_t>()[j])); }
        check_order = chosen;
      }
      double best = std::numeric_limits<double>::infinity();
      int stale = 0;
      std::string best_state;
      for (int epoch = 0; epoch < options.epochs; ++epoch)
      {
        model.train();
        const std::size_t pool = options.mil ? fit_bags.size() : train.size();
        auto permutation = torch::randperm(static_cast<int64_t>(pool), torch::kInt64);
        const std::size_t draw = options.mil ? std::min(pool, static_cast<std::size_t>(options.mil_bags_per_epoch)) : pool;
        Indices order;
        for (std::size_t j = 0; j < draw; ++j)
        {
          const auto picked = static_cast<std::size_t>(permutation.data_ptr<int64_t>()[j]);
          order.push_back(options.mil ? picked : train[picked]);
        }
        for (std::size_t start = 0; start < order.size(); start += options.batch_size)
        {
          const auto end = std::min(start + options.batch_size, order.size());
          optimizer.zero_grad();
          torch::Tensor value;
          if (options.mil)
          {
            torch::Tensor labels;
            value = torch::nn::functional::binary_cross_entropy_with_logits(bagLogits(model, data, fit_bags, order, start, end, labels), labels);
          }
          else
          {
            const auto b = topazBatch(data, order, start, end);
            value = torch::nn::functional::binary_cross_entropy_with_logits(model.forward(b), b.labels);
          }
          require(std::isfinite(value.item<double>()), "nonfinite training loss");
          value.backward();
          optimizer.step();
        }
        model.eval();
        double value = 0;
        const std::size_t checked = options.mil ? check_order.size() : validation.size();
        {
          torch::NoGradGuard guard;
          for (std::size_t start = 0; start < checked; start += options.batch_size)
          {
            const auto end = std::min(start + options.batch_size, checked);
            if (options.mil)
            {
              torch::Tensor labels;
              value += torch::nn::functional::binary_cross_entropy_with_logits(bagLogits(model, data, check_bags, check_order, start, end, labels), labels)
                         .item<double>() * (end - start);
            }
            else
            {
              const auto b = topazBatch(data, validation, start, end);
              value += torch::nn::functional::binary_cross_entropy_with_logits(model.forward(b), b.labels).item<double>() * (end - start);
            }
          }
        }
        value /= checked;
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
      std::istringstream restored(best_state, std::ios::binary);
      torch::serialize::InputArchive snapshot;
      snapshot.load_from(restored, torch::kCPU);
      model.load(snapshot);
      result.validation_loss = best;
      // Outer-fold models score only their held-out families (the only rows crossfit uses); others stay NaN.
      predictTopaz(model, data, result, options.outer_folds > 0);
      torch::serialize::OutputArchive archive, model_archive;
      model.save(model_archive);
      archive.write("model", model_archive);
      std::string metadata = TOPAZ_VERSION;
      for (const auto& [group, role] : provenance) { metadata += role + "\t" + group + "\n"; }
      archive.write("metadata", textTensor(metadata), true);
      archive.write("architecture", textTensor(options.mil ? "topaz-mil" : "topaz"), true);
      archive.write("topaz_width", torch::tensor({static_cast<int64_t>(options.topaz_width)}, torch::kInt64), true);
      archive.write("geometry", torch::tensor({geometry.mz_tolerance, geometry.rt_half_width, geometry.im_half_width}, torch::kFloat64), true);
      archive.write("training", torch::tensor({static_cast<double>(result.epochs_run), result.validation_loss}, torch::kFloat64), true);
      archive.write("training_options", torch::tensor({static_cast<int64_t>(options.seed), static_cast<int64_t>(options.threads),
                    static_cast<int64_t>(options.epochs), static_cast<int64_t>(options.patience), static_cast<int64_t>(options.batch_size)}, torch::kInt64), true);
      archive.write("mil_bags_per_epoch", torch::tensor({static_cast<int64_t>(options.mil_bags_per_epoch)}, torch::kInt64), true);
      archive.write("torch_version", textTensor(TORCH_VERSION), true);
      archive.write("outer_partition", torch::tensor({static_cast<int64_t>(options.outer_folds),
                    static_cast<int64_t>(options.heldout_fold)}, torch::kInt64), true);
      archive.save_to(checkpoint);
      return result;
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
    require(!options.mil || (options.architecture == CNNArchitecture::TOPAZ && options.mil_bags_per_epoch > 0),
            "multiple-instance training requires the TOPAZ architecture and a positive bag count");
    require(options.topaz_width >= 2 && options.topaz_width <= 64 && std::isfinite(options.weight_decay) && options.weight_decay >= 0 &&
            (options.architecture == CNNArchitecture::TOPAZ || (options.topaz_width == 32 && options.weight_decay == 0)),
            "invalid TOPAZ width or weight decay");
    at::set_num_threads(options.threads);
    torch::manual_seed(options.seed);
    at::globalContext().setDeterministicAlgorithms(true, false);
    if (options.architecture == CNNArchitecture::TOPAZ) { return trainTopaz(data, geometry, options, checkpoint); }
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
    torch::Tensor architecture;
    const bool topaz = archive.try_read("architecture", architecture, true);
    if (topaz)
    {
      const auto name = tensorText(architecture);
      require(name == "topaz" || name == "topaz-mil", "unknown checkpoint architecture");
      validateMS1(data);
    }
    const std::string version = topaz ? TOPAZ_VERSION : VERSION;
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
    archive.read("model", model_archive);
    if (topaz)
    {
      torch::Tensor width;
      int64_t saved_width = 32;
      if (archive.try_read("topaz_width", width, true))
      {
        require(width.scalar_type() == torch::kInt64 && width.numel() == 1, "invalid checkpoint TOPAZ width");
        saved_width = width.item<int64_t>();
        require(saved_width >= 2 && saved_width <= 64, "invalid checkpoint TOPAZ width");
      }
      TopazModel model(saved_width);
      model.load(model_archive);
      predictTopaz(model, data, result);
      return result;
    }
    Model model;
    model.load(model_archive);
    predict(model, data, result);
    return result;
  }
}
