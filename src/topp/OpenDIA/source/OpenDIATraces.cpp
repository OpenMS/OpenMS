// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIATraces.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    void require(bool valid, const std::string& message)
    {
      if (!valid)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
      }
    }

    bool safeText(const std::string& text)
    {
      return !text.empty() && text.find_first_of("\t\r\n") == std::string::npos;
    }
  }

  std::vector<TraceWindow> readTraceWindows(const std::string& path, std::size_t max_samples)
  {
    std::ifstream input(path);
    require(input.good(), "Could not open candidate windows: " + path);
    std::string line;
    require(static_cast<bool>(std::getline(input, line)), "Candidate window table is empty.");
    if (!line.empty() && line.back() == '\r')
    {
      line.pop_back();
    }
    std::vector<std::string> columns;
    StringUtils::split(line, '\t', columns);
    const std::array<std::string, 5> names = {"sample_id", "precursor_id", "group_id", "label", "rt_apex"};
    std::array<std::size_t, 5> indices{};
    for (std::size_t i = 0; i < names.size(); ++i)
    {
      require(std::count(columns.begin(), columns.end(), names[i]) == 1, "Expected one column: " + names[i]);
      indices[i] = std::find(columns.begin(), columns.end(), names[i]) - columns.begin();
    }
    std::vector<TraceWindow> windows;
    while (std::getline(input, line))
    {
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      require(!line.empty(), "Blank candidate window row.");
      std::vector<std::string> values;
      StringUtils::split(line, '\t', values);
      require(values.size() == columns.size(), "Candidate window column count differs from header.");
      require(windows.size() < max_samples, "Candidate windows exceed max_samples. Select a representative subset explicitly.");
      windows.push_back({values[indices[0]], values[indices[1]], values[indices[2]],
                         StringUtils::toInt32(values[indices[3]]), StringUtils::toDouble(values[indices[4]])});
    }
    require(input.eof(), "Error reading candidate windows.");
    return windows;
  }

  TraceCollector::TraceCollector(std::vector<TraceWindow> windows, const OpenSwath::LightTargetedExperiment& library,
                                 double mz_tolerance, double rt_half_width, double im_half_width,
                                 std::size_t max_samples, std::size_t max_tensor_bytes) :
    mz_tolerance_(mz_tolerance), rt_half_width_(rt_half_width), im_half_width_(im_half_width)
  {
    require(std::isfinite(mz_tolerance_) && mz_tolerance_ > 0.0 &&
            std::isfinite(rt_half_width_) && rt_half_width_ > 0.0 &&
            std::isfinite(im_half_width_) && im_half_width_ > 0.0, "Trace tolerances must be finite and positive.");
    // ponytail: sampled tensors are retained under an explicit byte budget.
    // Disk-backed tiled accumulation is only needed if a run cannot fit in RAM.
    require(!windows.empty() && windows.size() <= max_samples &&
            windows.size() <= max_tensor_bytes / (SAMPLE_SIZE * sizeof(float)),
            "Trace sample count is empty or exceeds max_samples / the raw tensor byte budget.");
    std::map<std::string, std::vector<const OpenSwath::LightTransition*>> transitions;
    std::map<std::string, double> mobilities;
    std::map<std::string, int> charges;
    for (const auto& compound : library.getCompounds())
    {
      mobilities[compound.id] = compound.drift_time;
      charges[compound.id] = compound.charge;
    }
    for (const auto& transition : library.getTransitions())
    {
      transitions[transition.peptide_ref].push_back(&transition);
    }
    std::set<std::string> ids;
    for (auto& window : windows)
    {
      require(safeText(window.sample_id) && safeText(window.precursor_id) && safeText(window.group_id) &&
              ids.insert(window.sample_id).second && std::isfinite(window.apex_rt) && window.apex_rt >= 0.0 &&
              window.label >= -1 && window.label <= 1, "Invalid/duplicate candidate sample, group, label or RT.");
      const auto found = transitions.find(window.precursor_id);
      require(found != transitions.end() && !found->second.empty(), "Unknown precursor: " + window.precursor_id);
      auto& fragments = found->second;
      for (const auto* fragment : fragments)
      {
        require(std::isfinite(fragment->product_mz) && fragment->product_mz > 0.0 &&
                std::isfinite(fragment->precursor_mz) && fragment->precursor_mz > 0.0 &&
                std::isfinite(fragment->library_intensity) && fragment->library_intensity >= 0.0 &&
                fragment->precursor_mz == fragments.front()->precursor_mz &&
                fragment->getDecoy() == fragments.front()->getDecoy(), "Invalid or inconsistent library fragment.");
      }
      std::stable_sort(fragments.begin(), fragments.end(), [](const auto* a, const auto* b)
      {
        if (a->library_intensity != b->library_intensity)
        {
          return a->library_intensity > b->library_intensity;
        }
        return a->product_mz < b->product_mz;
      });
      Sample sample;
      sample.window = std::move(window);
      sample.precursor_mz = fragments.front()->precursor_mz;
      sample.is_decoy = fragments.front()->getDecoy();
      require(!sample.is_decoy || sample.window.label != 1, "A decoy cannot be a positive seed.");
      const auto mobility = mobilities.find(sample.window.precursor_id);
      if (mobility != mobilities.end() && std::isfinite(mobility->second) && mobility->second > 0.0)
      {
        sample.precursor_im = mobility->second;
      }
      const auto charge = charges.find(sample.window.precursor_id);
      sample.precursor_charge = charge != charges.end() ? charge->second : 0;
      for (std::size_t i = 0; i < std::min(FRAGMENTS, fragments.size()); ++i)
      {
        sample.fragment_mz[i] = fragments[i]->product_mz;
        sample.library_intensity[i] = fragments[i]->library_intensity;
      }
      samples_.push_back(std::move(sample));
    }
    intensities_.resize(samples_.size() * SAMPLE_SIZE, 0.0f);
    by_rt_.resize(samples_.size());
    std::iota(by_rt_.begin(), by_rt_.end(), 0);
    std::stable_sort(by_rt_.begin(), by_rt_.end(), [this](std::size_t a, std::size_t b)
    { return samples_[a].window.apex_rt < samples_[b].window.apex_rt; });
  }

  void TraceCollector::enableEvidence(double ms1_ppm, double core_half_width_s)
  {
    require(std::isfinite(ms1_ppm) && ms1_ppm > 0.0 && std::isfinite(core_half_width_s) && core_half_width_s > 0.0,
            "Evidence tolerances must be finite and positive.");
    ms1_ppm_ = ms1_ppm;
    core_half_width_ = core_half_width_s;
    evidence_.assign(samples_.size(), TraceEvidence{});
  }

  std::size_t TraceCollector::rtBin_(double rt, double apex_rt) const
  {
    return std::min(RT_BINS - 1, static_cast<std::size_t>(
      std::max(0.0, (rt - apex_rt + rt_half_width_) / (2.0 * rt_half_width_)) * RT_BINS));
  }

  void TraceCollector::consumeMS1_(MSSpectrum& spectrum)
  {
    const double rt = spectrum.getRT();
    auto first = std::lower_bound(by_rt_.begin(), by_rt_.end(), rt - rt_half_width_,
                                  [this](std::size_t i, double value) { return samples_[i].window.apex_rt < value; });
    if (first == by_rt_.end() || samples_[*first].window.apex_rt > rt + rt_half_width_)
    {
      return;
    }
    if (!spectrum.isSorted())
    {
      spectrum.sortByPosition();
    }
    // MS1 points are summed over ion mobility; the precursor m/z window is narrow enough for that on DIA data.
    const auto sum_window = [&](double target, double& ppm_sum, double& weight)
    {
      const double half = target * ms1_ppm_ * 1e-6;
      double total = 0.0;
      for (auto point = spectrum.MZBegin(target - half); point != spectrum.end() && point->getMZ() <= target + half; ++point)
      {
        const double intensity = point->getIntensity();
        if (!std::isfinite(intensity) || intensity <= 0.0) continue;
        total += intensity;
        ppm_sum += intensity * (point->getMZ() - target) / target * 1e6;
        weight += intensity;
      }
      return total;
    };
    for (auto active = first; active != by_rt_.end() && samples_[*active].window.apex_rt <= rt + rt_half_width_; ++active)
    {
      const Sample& sample = samples_[*active];
      TraceEvidence& evidence = evidence_[*active];
      const std::size_t bin = rtBin_(rt, sample.window.apex_rt);
      evidence.ms1_mask |= std::uint32_t{1} << bin;
      const bool core = std::abs(rt - sample.window.apex_rt) <= core_half_width_;
      double ppm_sum = 0.0, weight = 0.0;
      const double mono = sum_window(sample.precursor_mz, ppm_sum, weight);
      evidence.ms1[bin] += static_cast<float>(mono);
      if (core)
      {
        evidence.ms1_ppm_sum += ppm_sum;
        evidence.ms1_weight += weight;
        if (sample.precursor_charge > 0)
        {
          double unused_sum = 0.0, unused_weight = 0.0;
          evidence.ms1_isotope += static_cast<float>(
            sum_window(sample.precursor_mz + 1.0033548 / sample.precursor_charge, unused_sum, unused_weight));
        }
      }
    }
  }

  void TraceCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() == 1 && !evidence_.empty() && !spectrum.empty() && std::isfinite(spectrum.getRT()))
    {
      consumeMS1_(spectrum);
      return;
    }
    if (spectrum.getMSLevel() != 2 || spectrum.empty() || spectrum.getPrecursors().empty() ||
        !std::isfinite(spectrum.getRT()))
    {
      return;
    }
    const double rt = spectrum.getRT();
    auto first = std::lower_bound(by_rt_.begin(), by_rt_.end(), rt - rt_half_width_,
                                  [this](std::size_t i, double value) { return samples_[i].window.apex_rt < value; });
    if (first == by_rt_.end() || samples_[*first].window.apex_rt > rt + rt_half_width_)
    {
      return;
    }
    if (!spectrum.isSorted())
    {
      spectrum.sortByPosition();
    }

    const DataArrays::FloatDataArray* mobilities = nullptr;
    if (spectrum.containsIMData())
    {
      const auto [index, unit] = spectrum.getIMData();
      if (spectrum.getFloatDataArrays()[index].size() == spectrum.size() &&
          (unit == DriftTimeUnit::VSSC || (ion_mobility_array_as_vssc_ && unit != DriftTimeUnit::CCS)))
      {
        mobilities = &spectrum.getFloatDataArrays()[index];
      }
    }
    const double scalar_im = spectrum.getDriftTime();
    const bool scalar_valid = spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC &&
                               std::isfinite(scalar_im) && scalar_im > 0.0;
    for (auto active = first; active != by_rt_.end() && samples_[*active].window.apex_rt <= rt + rt_half_width_; ++active)
    {
      const std::size_t index = *active;
      Sample& sample = samples_[index];
      bool isolated = false;
      for (const auto& precursor : spectrum.getPrecursors())
      {
        if (precursor.getIsolationWindowLowerOffset() >= 0.0 && precursor.getIsolationWindowUpperOffset() >= 0.0 &&
            precursor.getIsolationWindowLowerOffset() + precursor.getIsolationWindowUpperOffset() > 0.0 &&
            precursor.getMZ() - precursor.getIsolationWindowLowerOffset() <= sample.precursor_mz &&
            precursor.getMZ() + precursor.getIsolationWindowUpperOffset() >= sample.precursor_mz)
        {
          isolated = true;
          break;
        }
      }
      if (!isolated)
      {
        continue;
      }
      const bool compatible_im = sample.precursor_im > 0.0 && (mobilities != nullptr || scalar_valid);
      // Mixed availability cannot be represented by one sample mask. Collapse
      // such a sample to an XIC instead of inventing an IM position for missing data.
      if (!compatible_im && !sample.saw_missing_im)
      {
        for (std::size_t offset = index * SAMPLE_SIZE; offset < (index + 1) * SAMPLE_SIZE; offset += IM_BINS)
        {
          const double sum = std::accumulate(intensities_.begin() + offset, intensities_.begin() + offset + IM_BINS, 0.0);
          require(sum <= std::numeric_limits<float>::max(), "Collapsed trace intensity overflow.");
          std::fill_n(intensities_.begin() + offset, IM_BINS, 0.0f);
          intensities_[offset + IM_BINS / 2] = static_cast<float>(sum);
        }
        sample.saw_missing_im = true;
        sample.has_im = false;
      }
      if (compatible_im && !sample.saw_missing_im)
      {
        sample.has_im = true;
      }
      const auto rt_bin = std::min(RT_BINS - 1, static_cast<std::size_t>(
        std::max(0.0, (rt - sample.window.apex_rt + rt_half_width_) / (2.0 * rt_half_width_)) * RT_BINS));
      TraceEvidence* evidence = evidence_.empty() ? nullptr : &evidence_[index];
      const bool core = evidence != nullptr && std::abs(rt - sample.window.apex_rt) <= core_half_width_;
      if (evidence != nullptr)
      {
        evidence->ms2_mask |= std::uint32_t{1} << rt_bin;
      }
      for (std::size_t fragment = 0; fragment < FRAGMENTS; ++fragment)
      {
        const double mz = sample.fragment_mz[fragment];
        if (mz == 0.0)
        {
          continue;
        }
        for (auto point = spectrum.MZBegin(mz - mz_tolerance_);
             point != spectrum.end() && point->getMZ() <= mz + mz_tolerance_; ++point)
        {
          const double intensity = point->getIntensity();
          if (!std::isfinite(intensity) || intensity <= 0.0)
          {
            continue;
          }
          std::size_t im_bin = IM_BINS / 2;
          if (compatible_im)
          {
            const double im = mobilities != nullptr ? (*mobilities)[point - spectrum.begin()] : scalar_im;
            if (!std::isfinite(im) || im <= 0.0 || std::abs(im - sample.precursor_im) > im_half_width_)
            {
              continue;
            }
            if (!sample.saw_missing_im)
            {
              im_bin = std::min(IM_BINS - 1, static_cast<std::size_t>(
                std::max(0.0, (im - sample.precursor_im + im_half_width_) / (2.0 * im_half_width_)) * IM_BINS));
            }
          }
          float& value = intensities_[((index * FRAGMENTS + fragment) * RT_BINS + rt_bin) * IM_BINS + im_bin];
          require(intensity <= std::numeric_limits<float>::max() - static_cast<double>(value), "Trace intensity overflow.");
          value += static_cast<float>(intensity);
          if (core)
          {
            evidence->ms2_ppm_sum += intensity * (point->getMZ() - mz) / mz * 1e6;
            evidence->ms2_weight += intensity;
          }
        }
      }
    }
  }

  void TraceCollector::consumeChromatogram(MSChromatogram&) {}
  void TraceCollector::setExpectedSize(std::size_t, std::size_t) {}
  void TraceCollector::setExperimentalSettings(const ExperimentalSettings&) {}

  void TraceCollector::store(const std::string& prefix) const
  {
    require(std::endian::native == std::endian::little && sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
            "Raw tensor export requires an IEEE754 little-endian float32 platform.");
    std::ofstream data(prefix + ".traces.npy", std::ios::binary);
    std::ofstream metadata(prefix + ".samples.tsv");
    std::ofstream geometry(prefix + ".geometry.json");
    require(data.good() && metadata.good() && geometry.good(), "Could not open tensor output prefix: " + prefix);
    geometry << std::setprecision(17) << "{\"version\":1,\"fragments\":6,\"rt_bins\":32,\"im_bins\":8,"
             << "\"rt_half_width_seconds\":" << rt_half_width_ << ",\"im_half_width_vssc\":" << im_half_width_
             << ",\"mz_tolerance_da\":" << mz_tolerance_ << "}\n";
    // NumPy v1.0: explicit dtype/shape, 16-byte-aligned ASCII header, C-order data.
    std::string header = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + std::to_string(samples_.size()) +
                         ", 6, 32, 8), }";
    header.append((16 - (10 + header.size() + 1) % 16) % 16, ' ');
    header += '\n';
    const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0,
                          static_cast<char>(header.size() & 255), static_cast<char>((header.size() >> 8) & 255)};
    data.write(magic, sizeof(magic));
    data << header;
    data.write(reinterpret_cast<const char*>(intensities_.data()), static_cast<std::streamsize>(intensities_.size() * sizeof(float)));
    metadata << "sample_id\tprecursor_id\tgroup_id\tis_decoy\tlabel\trt_apex\thas_im\tlibrary_intensities\tfragment_mz\n"
             << std::setprecision(17);
    for (const Sample& sample : samples_)
    {
      metadata << sample.window.sample_id << '\t' << sample.window.precursor_id << '\t' << sample.window.group_id << '\t'
               << sample.is_decoy << '\t' << sample.window.label << '\t' << sample.window.apex_rt << '\t' << sample.has_im << '\t';
      for (std::size_t i = 0; i < FRAGMENTS; ++i)
      {
        metadata << (i ? "," : "") << sample.library_intensity[i];
      }
      metadata << '\t';
      for (std::size_t i = 0; i < FRAGMENTS; ++i)
      {
        metadata << (i ? "," : "") << sample.fragment_mz[i];
      }
      metadata << '\n';
    }
    data.flush();
    metadata.flush();
    geometry.flush();
    require(data.good() && metadata.good() && geometry.good(), "Failed to write raw tensor dataset: " + prefix);
  }
}
