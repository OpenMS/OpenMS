// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACNN.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

using namespace OpenMS::OpenDIA;
namespace
{
  constexpr std::size_t WINDOW = 1536;
  void check(bool condition, const char* message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }
  template <typename F> void rejects(F fn)
  {
    bool rejected = false;
    try { fn(); } catch (const std::exception&) { rejected = true; }
    check(rejected, "invalid input was accepted");
  }
  unsigned bucket(const std::string& group)
  {
    std::string key = "42";
    key.push_back('\0');
    key += group;
    std::uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char c : key) { hash = (hash ^ c) * UINT64_C(1099511628211); }
    return hash % 10;
  }
  void binary(const std::filesystem::path& path, const std::vector<float>& values)
  {
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(float)));
    check(static_cast<bool>(out), "failed to write fixture");
  }
}
int main(int argc, char** argv)
{
  try
  {
    const bool parity = argc == 3 && std::string(argv[1]) == "--parity";
    check(argc == 1 || parity, "usage: OpenDIACNN_test [--parity directory]");
    // Two labels in each family ensures both classes in every partition, while
    // identical library/mask and equal total raw intensity prevent aux shortcuts.
    CNNData data;
    std::vector<float> raw;
    std::map<unsigned, int> counts;
    for (int family = 0; counts[0] < 4 || counts[1] < 2 || counts[2] < 2; ++family)
    {
      const auto group = "family-" + std::to_string(family);
      const unsigned b = bucket(group), partition = b < 6 ? 0 : b < 8 ? 1 : 2;
      if (counts[partition] >= (partition == 0 ? 4 : 2)) { continue; }
      ++counts[partition];
      for (int label : {0, 1, -1})
      {
        CNNSample s;
        s.sample_id = group + "-" + std::to_string(label);
        s.group_id = group;
        s.label = label;
        s.has_im = true;
        s.library_intensities = {6, 5, 4, 3, 2, 1};
        s.fragment_mz = {400, 500, 600, 700, 800, 900};
        data.samples.push_back(s);
        const auto offset = raw.size();
        raw.resize(offset + WINDOW, 0);
        for (std::size_t c = 0; c < 6; ++c)
        {
          for (std::size_t rt = (label == 1 ? 4 : 20); rt < (label == 1 ? 12 : 28); ++rt)
          {
            for (std::size_t im = 0; im < 8; ++im) { raw[offset + c * 256 + rt * 8 + im] = 10; }
          }
        }
      }
    }
    data.intensities = raw.data();
    data.intensity_count = raw.size();
    CNNGeometry geometry;
    CNNOptions options;
    options.threads = 1;
    const auto directory = parity ? std::filesystem::path(argv[2]) :
      std::filesystem::temp_directory_path() / ("opendia-cnn-" + std::to_string(std::filesystem::file_time_type::clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    const auto checkpoint = (directory / "native.pt").string();
    auto result = trainTraceCNN(data, geometry, options, checkpoint);
    check(result.logits.size() == data.samples.size(), "wrong logit count");
    check(result.epochs_run > 0 && result.epochs_run <= 40, "epoch cap violated");
    check(result.validation_loss < 0.3 && result.heldout_auc > 0.99, "raw-only signal not learned");
    std::map<std::string, std::string> families;
    for (std::size_t i = 0; i < data.samples.size(); ++i)
    {
      auto [it, inserted] = families.emplace(data.samples[i].group_id, result.splits[i]);
      check(inserted || it->second == result.splits[i], "family leaked between partitions");
      const auto b = bucket(data.samples[i].group_id);
      check(result.splits[i] == (b < 6 ? "train" : b < 8 ? "validation" : "heldout"), "unstable split hash");
    }
    double validation_loss = 0;
    std::size_t validation_count = 0;
    for (std::size_t i = 0; i < data.samples.size(); ++i)
    {
      if (result.splits[i] == "validation" && data.samples[i].label >= 0)
      {
        const double logit = result.logits[i];
        validation_loss += std::max(logit, 0.0) - logit * data.samples[i].label + std::log1p(std::exp(-std::abs(logit)));
        ++validation_count;
      }
    }
    check(std::abs(validation_loss / validation_count - result.validation_loss) < 1e-6, "best validation weights not restored");
    const auto repeat_checkpoint = (directory / "repeat.pt").string();
    auto changed_unknown_raw = raw;
    for (std::size_t i = 0; i < data.samples.size(); ++i)
    {
      if (data.samples[i].label == -1)
      {
        std::fill(changed_unknown_raw.begin() + i * WINDOW, changed_unknown_raw.begin() + (i + 1) * WINDOW, 7.f);
      }
    }
    auto changed_unknown = data;
    changed_unknown.intensities = changed_unknown_raw.data();
    const auto repeated = trainTraceCNN(changed_unknown, geometry, options, repeat_checkpoint);
    check(repeated.epochs_run == result.epochs_run && repeated.validation_loss == result.validation_loss,
          "training is nondeterministic or unknown targets entered loss");
    for (std::size_t i = 0; i < data.samples.size(); ++i)
    {
      if (data.samples[i].label >= 0) { check(repeated.logits[i] == result.logits[i], "unknown rows influenced model"); }
    }
    std::filesystem::remove(repeat_checkpoint);
    rejects([&] { scoreTraceCNN(data, geometry, checkpoint, 1); });
    const auto restored = scoreTraceCNN(data, geometry, checkpoint, 1, true);
    check(restored.splits == result.splits && restored.logits == result.logits, "checkpoint roundtrip changed output");
    auto changed = geometry;
    changed.mz_tolerance *= 2;
    rejects([&] { scoreTraceCNN(data, changed, checkpoint, 1, true); });
    changed = geometry;
    changed.rt_half_width = std::numeric_limits<double>::infinity();
    rejects([&] { scoreTraceCNN(data, changed, checkpoint, 1, true); });
    auto invalid = data;
    --invalid.intensity_count;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[0].sample_id = invalid.samples[1].sample_id;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[0].group_id = "bad\tgroup";
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[1].is_decoy = true;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[0].label = 2;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[0].library_intensities[1] = 20;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    invalid = data;
    invalid.samples[0].fragment_mz[0] = 0;
    rejects([&] { scoreTraceCNN(invalid, geometry, checkpoint, 1, true); });
    raw[0] = -1;
    rejects([&] { scoreTraceCNN(data, geometry, checkpoint, 1, true); });
    raw[0] = std::numeric_limits<float>::quiet_NaN();
    rejects([&] { scoreTraceCNN(data, geometry, checkpoint, 1, true); });
    raw[0] = 0;
    auto bad_options = options;
    bad_options.epochs = 101;
    rejects([&] { trainTraceCNN(data, geometry, bad_options, checkpoint); });
    // Labels and provenance do not enter inference. Unknown rows retain original
    // family provenance even when their sample IDs have never been seen before.
    auto unknown = data;
    for (auto& s : unknown.samples) { s.label = -1; s.sample_id += "-new"; s.is_decoy = true; }
    const auto unknown_result = scoreTraceCNN(unknown, geometry, checkpoint, 1, true);
    check(unknown_result.logits == result.logits && unknown_result.splits == result.splits, "labels/provenance entered network");
    for (auto& s : unknown.samples) { s.group_id = "unseen-family"; }
    const auto external = scoreTraceCNN(unknown, geometry, checkpoint, 1);
    check(external.logits == result.logits && std::all_of(external.splits.begin(), external.splits.end(), [](const auto& s) { return s == "external"; }), "external scoring failed");
    // Outer evaluation families must affect neither fitted weights nor the
    // validation loss used to select an epoch, even when their labels/raw change.
    CNNData outer_data;
    std::vector<float> outer_raw;
    for (int family = 0; family < 60; ++family)
    {
      for (int label : {0, 1})
      {
        auto sample = data.samples.at(static_cast<std::size_t>(label));
        sample.group_id = "outer-family-" + std::to_string(family);
        sample.sample_id = sample.group_id + "-" + std::to_string(label);
        sample.label = label;
        outer_data.samples.push_back(sample);
        outer_raw.insert(outer_raw.end(), raw.begin() + label * WINDOW, raw.begin() + (label + 1) * WINDOW);
      }
    }
    outer_data.intensities = outer_raw.data();
    outer_data.intensity_count = outer_raw.size();
    auto outer_options = options;
    outer_options.outer_folds = 3;
    outer_options.epochs = 3;
    for (int fold = 0; fold < 3; ++fold)
    {
      outer_options.heldout_fold = fold;
      const auto path = (directory / ("outer-" + std::to_string(fold) + ".pt")).string();
      const auto fitted = trainTraceCNN(outer_data, geometry, outer_options, path);
      auto perturbed = outer_data;
      auto perturbed_raw = outer_raw;
      std::size_t heldout_count = 0;
      for (std::size_t i = 0; i < outer_data.samples.size(); ++i)
      {
        const bool heldout = traceCNNFold(outer_data.samples[i].group_id, 42, 3) == fold;
        check((fitted.splits[i] == "heldout") == heldout, "outer family entered train or validation partition");
        if (heldout)
        {
          ++heldout_count;
          perturbed.samples[i].label = 1 - perturbed.samples[i].label;
          std::fill(perturbed_raw.begin() + i * WINDOW, perturbed_raw.begin() + (i + 1) * WINDOW, 3.f);
        }
      }
      check(heldout_count > 0, "outer fixture has no evaluation families");
      perturbed.intensities = perturbed_raw.data();
      const auto repeated_outer = trainTraceCNN(perturbed, geometry, outer_options, repeat_checkpoint);
      check(fitted.validation_loss == repeated_outer.validation_loss && fitted.epochs_run == repeated_outer.epochs_run,
            "outer-heldout data influenced early stopping");
      for (std::size_t i = 0; i < outer_data.samples.size(); ++i)
      {
        if (fitted.splits[i] != "heldout") { check(fitted.logits[i] == repeated_outer.logits[i], "outer-heldout data influenced fitting"); }
      }
      const auto roundtrip = scoreTraceCNN(outer_data, geometry, path, 1, true);
      check(roundtrip.splits == fitted.splits && roundtrip.logits == fitted.logits, "outer checkpoint lost partition or weights");
      std::filesystem::remove(path);
    }
    std::filesystem::remove(repeat_checkpoint);
    rejects([&] { traceCNNFold("family", -1, 3); });
    rejects([&] { traceCNNFold("family", 42, 1); });
    rejects([&] { traceCNNFold("bad\tfamily", 42, 3); });
    outer_options.heldout_fold = 3;
    rejects([&] { trainTraceCNN(outer_data, geometry, outer_options, repeat_checkpoint); });
    if (parity)
    {
      CNNData fixture;
      fixture.samples.assign(data.samples.begin(), data.samples.begin() + 3);
      fixture.intensities = raw.data();
      fixture.intensity_count = 3 * WINDOW;
      fixture.samples[0].has_im = false;
      fixture.samples[2].library_intensities.fill(0);
      const auto scores = scoreTraceCNN(fixture, geometry, checkpoint, 1, true);
      binary(directory / "raw.f32", std::vector<float>(raw.begin(), raw.begin() + 3 * WINDOW));
      std::vector<float> aux;
      for (const auto& s : fixture.samples)
      {
        const float maximum = static_cast<float>(std::max(s.library_intensities[0], 1e-12));
        float sum = 0;
        for (double v : s.library_intensities) { sum += static_cast<float>(v) / maximum; }
        for (double v : s.library_intensities) { aux.push_back((static_cast<float>(v) / maximum) / std::max(sum, 1e-12f)); }
        aux.push_back(s.has_im ? 1.f : 0.f);
      }
      binary(directory / "aux.f32", aux);
      binary(directory / "logits.f32", scores.logits);
    }
    else { std::filesystem::remove_all(directory); }
    std::cout << "Native CNN tests passed, validation loss=" << result.validation_loss << ", heldout AUC=" << result.heldout_auc << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
