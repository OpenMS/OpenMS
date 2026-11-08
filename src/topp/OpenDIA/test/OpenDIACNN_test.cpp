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

    // TOPAZ architecture. Positives co-elute in all six fragments; negatives have the same per-channel
    // intensity but at a different RT per fragment, so only co-elution separates them.
    {
      CNNData topaz;
      std::vector<float> topaz_raw;
      for (int family = 0; family < 80; ++family)
      {
        for (int label : {0, 1})
        {
          CNNSample s;
          s.group_id = "topaz-family-" + std::to_string(family);
          s.sample_id = s.group_id + "-" + std::to_string(label);
          s.label = label;
          s.is_decoy = label == 0;
          s.bag_id = s.sample_id;
          s.library_intensities = {6, 5, 4, 3, 2, 1};
          s.fragment_mz = {400, 500, 600, 700, 800, 900};
          topaz.samples.push_back(s);
          const auto offset = topaz_raw.size();
          topaz_raw.resize(offset + WINDOW, 0);
          for (std::size_t c = 0; c < 6; ++c)
          {
            const std::size_t apex = label == 1 ? 16 : 4 + 4 * c + family % 3;
            for (std::size_t rt = apex - 2; rt <= apex + 2; ++rt)
            {
              topaz_raw[offset + c * 256 + rt * 8 + 4] = 100.f * (3 - std::abs(static_cast<int>(rt) - static_cast<int>(apex)));
            }
          }
        }
      }
      topaz.intensities = topaz_raw.data();
      topaz.intensity_count = topaz_raw.size();
      CNNOptions topaz_options;
      topaz_options.threads = 1;
      topaz_options.epochs = 15;
      topaz_options.architecture = CNNArchitecture::TOPAZ;
      const auto topaz_path = (directory / "topaz.pt").string();
      const auto fitted = trainTraceCNN(topaz, geometry, topaz_options, topaz_path);
      check(fitted.heldout_auc > 0.99, "TOPAZ did not learn co-elution");
      const auto reloaded = scoreTraceCNN(topaz, geometry, topaz_path, 1, true);
      check(reloaded.logits == fitted.logits && reloaded.splits == fitted.splits, "TOPAZ checkpoint roundtrip changed output");
      // Absent modalities: with has_ms1 = has_im = false, MS1 traces and the IM distribution must not change any logit.
      std::vector<float> ms1(topaz.samples.size() * 96, 0.f);
      auto with_ms1 = topaz;
      with_ms1.ms1 = ms1.data();
      with_ms1.ms1_count = ms1.size();
      std::vector<float> noisy_ms1(ms1.size());
      for (std::size_t j = 0; j < noisy_ms1.size(); ++j) { noisy_ms1[j] = static_cast<float>((j * 7919) % 97); }
      auto with_noise = with_ms1;
      with_noise.ms1 = noisy_ms1.data();
      check(scoreTraceCNN(with_noise, geometry, topaz_path, 1, true).logits == fitted.logits, "masked MS1 traces changed logits");
      auto spread_raw = topaz_raw;
      for (std::size_t j = 0; j < spread_raw.size(); j += 8)
      {
        // Move the signal to another IM bin: the XIC (IM sum) is unchanged, only the mobilogram differs.
        std::swap(spread_raw[j + 4], spread_raw[j + 1]);
      }
      auto spread = topaz;
      spread.intensities = spread_raw.data();
      check(scoreTraceCNN(spread, geometry, topaz_path, 1, true).logits == fitted.logits, "masked mobilogram changed logits");
      for (auto& s : spread.samples) { s.has_im = true; }
      auto unmoved = topaz;
      for (auto& s : unmoved.samples) { s.has_im = true; }
      check(scoreTraceCNN(spread, geometry, topaz_path, 1, true).logits != scoreTraceCNN(unmoved, geometry, topaz_path, 1, true).logits,
            "present mobilogram did not reach the network");
      auto bad_ms1 = with_ms1;
      bad_ms1.ms1_count -= 1;
      rejects([&] { scoreTraceCNN(bad_ms1, geometry, topaz_path, 1, true); });
      auto claims_ms1 = topaz;
      claims_ms1.samples[0].has_ms1 = true;
      rejects([&] { scoreTraceCNN(claims_ms1, geometry, topaz_path, 1, true); });
      // A grid2d checkpoint still scores the same windows (it ignores MS1, bags and has_ms1).
      check(scoreTraceCNN(topaz, geometry, checkpoint, 1, true).logits.size() == topaz.samples.size(), "grid2d checkpoint failed on new data");
      // MIL: bags of two windows (the positive and a scrambled copy) per target; decoy bags hold scrambled windows only.
      auto bags = topaz;
      for (std::size_t i = 0; i < bags.samples.size(); ++i)
      {
        bags.samples[i].label = -1;
        bags.samples[i].bag_id = bags.samples[i].group_id + (bags.samples[i].is_decoy ? "-decoy" : "-target");
      }
      // Pair each target bag with a scrambled member: reuse the decoy window as a second target candidate.
      CNNData mil_data = bags;
      std::vector<float> mil_raw = topaz_raw;
      for (std::size_t i = 0; i < bags.samples.size(); ++i)
      {
        if (!bags.samples[i].is_decoy) { continue; }
        auto extra = bags.samples[i];
        extra.is_decoy = false;
        extra.sample_id += "-target-candidate";
        extra.bag_id = extra.group_id + "-target";
        mil_data.samples.push_back(extra);
        mil_raw.insert(mil_raw.end(), topaz_raw.begin() + i * WINDOW, topaz_raw.begin() + (i + 1) * WINDOW);
      }
      mil_data.intensities = mil_raw.data();
      mil_data.intensity_count = mil_raw.size();
      auto mil_options = topaz_options;
      mil_options.mil = true;
      mil_options.mil_bags_per_epoch = 64;
      const auto mil_path = (directory / "topaz-mil.pt").string();
      const auto mil_fit = trainTraceCNN(mil_data, geometry, mil_options, mil_path);
      // Score windows: co-eluting target windows must outrank scrambled windows in held-out families.
      std::vector<float> good, scrambled;
      for (std::size_t i = 0; i < mil_data.samples.size(); ++i)
      {
        if (mil_fit.splits[i] != "heldout") { continue; }
        (mil_data.samples[i].sample_id.ends_with("-1") ? good : scrambled).push_back(mil_fit.logits[i]);
      }
      check(!good.empty() && !scrambled.empty() && *std::min_element(good.begin(), good.end()) > *std::max_element(scrambled.begin(), scrambled.end()),
            "MIL did not rank co-eluting candidates first");
      auto grid_mil = mil_options;
      grid_mil.architecture = CNNArchitecture::GRID2D;
      rejects([&] { trainTraceCNN(mil_data, geometry, grid_mil, mil_path + ".grid"); });
      auto mixed = mil_data;
      mixed.samples.back().is_decoy = true;
      rejects([&] { trainTraceCNN(mixed, geometry, mil_options, mil_path + ".mixed"); });
      // Compact width with weight decay still learns co-elution and roundtrips with its saved width.
      auto compact_options = topaz_options;
      compact_options.topaz_width = 8;
      compact_options.weight_decay = 0.01;
      const auto compact_path = (directory / "topaz-compact.pt").string();
      const auto compact = trainTraceCNN(topaz, geometry, compact_options, compact_path);
      check(compact.heldout_auc > 0.99, "compact TOPAZ did not learn co-elution");
      check(scoreTraceCNN(topaz, geometry, compact_path, 1, true).logits == compact.logits, "compact TOPAZ checkpoint roundtrip changed output");
      auto bad_width = compact_options;
      bad_width.topaz_width = 1;
      rejects([&] { trainTraceCNN(topaz, geometry, bad_width, compact_path + ".bad"); });
      auto grid_decay = options;
      grid_decay.weight_decay = 0.01;
      rejects([&] { trainTraceCNN(data, geometry, grid_decay, compact_path + ".grid"); });
      // TOPAZ-faithful regime: learns co-elution, roundtrips, and (block max normalisation) is invariant to per-window scale.
      auto faithful_options = compact_options;
      faithful_options.topaz_faithful = true;
      const auto faithful_path = (directory / "topaz-faithful.pt").string();
      const auto faithful = trainTraceCNN(topaz, geometry, faithful_options, faithful_path);
      check(faithful.heldout_auc > 0.99, "faithful TOPAZ did not learn co-elution");
      const auto faithful_scored = scoreTraceCNN(topaz, geometry, faithful_path, 1, true);
      check(faithful_scored.logits == faithful.logits, "faithful TOPAZ checkpoint roundtrip changed output");
      auto scaled_raw = topaz_raw;
      for (std::size_t i = 0; i < topaz.samples.size(); ++i)
      {
        for (std::size_t j = 0; j < WINDOW; ++j) { scaled_raw[i * WINDOW + j] *= static_cast<float>(1 << (i % 8)); }
      }
      auto scaled = topaz;
      scaled.intensities = scaled_raw.data();
      const auto scaled_scored = scoreTraceCNN(scaled, geometry, faithful_path, 1, true);
      double worst = 0;
      for (std::size_t i = 0; i < scaled_scored.logits.size(); ++i) { worst = std::max(worst, std::abs(double(scaled_scored.logits[i]) - faithful.logits[i])); }
      check(worst < 1e-3, "faithful TOPAZ is not invariant to per-window intensity scale");
      check(scoreTraceCNN(scaled, geometry, compact_path, 1, true).logits != compact.logits, "raw-scale TOPAZ unexpectedly scale invariant");
      auto grid_faithful = options;
      grid_faithful.topaz_faithful = true;
      rejects([&] { trainTraceCNN(data, geometry, grid_faithful, faithful_path + ".grid"); });
      auto faithful_mil = mil_options;
      faithful_mil.topaz_faithful = true;
      const auto faithful_mil_fit = trainTraceCNN(mil_data, geometry, faithful_mil, faithful_path + ".mil");
      std::vector<float> fgood, fscrambled;
      for (std::size_t i = 0; i < mil_data.samples.size(); ++i)
      {
        if (faithful_mil_fit.splits[i] != "heldout") { continue; }
        (mil_data.samples[i].sample_id.ends_with("-1") ? fgood : fscrambled).push_back(faithful_mil_fit.logits[i]);
      }
      check(*std::min_element(fgood.begin(), fgood.end()) > *std::max_element(fscrambled.begin(), fscrambled.end()),
            "faithful MIL did not rank co-eluting candidates first");
    }
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
