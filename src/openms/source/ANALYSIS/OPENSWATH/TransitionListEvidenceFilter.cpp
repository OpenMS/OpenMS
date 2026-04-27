// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflowScheduler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  namespace
  {
    struct PeakEvidence
    {
      double mz{0.0};
      double intensity{0.0};
      double im{-1.0};
    };

    struct FragmentCandidate
    {
      double product_mz{0.0};
      double library_intensity{0.0};
      Size slot{0};
      std::string transition_name;
    };

    struct PrecursorCandidate
    {
      std::string compound_id;
      OpenSwath::LightCompound compound;
      double precursor_mz{0.0};
      double precursor_im{-1.0};
      std::vector<OpenSwath::LightTransition> transitions;
      std::vector<FragmentCandidate> fragments;
    };

    struct PrecursorIMTransform
    {
      double scale{1.0};
      bool multiply_by_charge{false};

      bool isIdentity() const
      {
        return scale == 1.0 && !multiply_by_charge;
      }

      double factor(const PrecursorCandidate& candidate) const
      {
        double factor = scale;
        if (multiply_by_charge && candidate.compound.charge > 0)
        {
          factor *= candidate.compound.charge;
        }
        return factor;
      }
    };

    struct MzIndexEntry
    {
      double mz{0.0};
      Size candidate_index{0};
    };

    struct ProductIndexEntry
    {
      double mz{0.0};
      Size candidate_index{0};
      Size fragment_slot{0};
    };

    struct EvidenceAccum
    {
      Size ms1_hit_count{0};
      double ms1_max_intensity{0.0};
      double ms1_sum_intensity{0.0};
      double ms1_best_rt{-1.0};
      Size ms2_hit_count{0};
      Size ms2_best_fragment_hits{0};
      double ms2_max_intensity{0.0};
      double ms2_sum_intensity{0.0};
      double ms2_best_rt{-1.0};

      void merge(const EvidenceAccum& rhs)
      {
        ms1_hit_count += rhs.ms1_hit_count;
        ms1_sum_intensity += rhs.ms1_sum_intensity;
        if (rhs.ms1_max_intensity > ms1_max_intensity)
        {
          ms1_max_intensity = rhs.ms1_max_intensity;
          ms1_best_rt = rhs.ms1_best_rt;
        }

        ms2_hit_count += rhs.ms2_hit_count;
        ms2_sum_intensity += rhs.ms2_sum_intensity;
        ms2_best_fragment_hits = std::max(ms2_best_fragment_hits, rhs.ms2_best_fragment_hits);
        if (rhs.ms2_max_intensity > ms2_max_intensity)
        {
          ms2_max_intensity = rhs.ms2_max_intensity;
          ms2_best_rt = rhs.ms2_best_rt;
        }
      }
    };

    struct SpectrumMS2Hit
    {
      std::uint64_t fragment_mask{0};
      double sum_intensity{0.0};
      double max_intensity{0.0};
    };

    Size countBits_(std::uint64_t value)
    {
      Size count = 0;
      while (value != 0)
      {
        value &= (value - 1);
        ++count;
      }
      return count;
    }

    int resolveThreadCount_(int threads)
    {
      if (threads == 0)
      {
        return 1;
      }
#ifdef _OPENMP
      const int max_threads = std::max(1, omp_get_max_threads());
      if (threads < 0)
      {
        return max_threads;
      }
      return std::max(1, std::min(threads, max_threads));
#else
      return 1;
#endif
    }

    double halfWindow_(double center_mz, const ChromExtractParams& params)
    {
      if (params.mz_extraction_window <= 0.0)
      {
        return 0.0;
      }
      if (params.ppm)
      {
        return center_mz * params.mz_extraction_window * 0.5e-6;
      }
      return params.mz_extraction_window * 0.5;
    }

    bool withinMzTolerance_(double observed_mz, double target_mz, const ChromExtractParams& params)
    {
      if (params.mz_extraction_window <= 0.0)
      {
        return false;
      }
      const double delta = std::fabs(observed_mz - target_mz);
      if (params.ppm)
      {
        return delta <= target_mz * params.mz_extraction_window * 0.5e-6;
      }
      return delta <= params.mz_extraction_window * 0.5;
    }

    bool withinIMTolerance_(double observed_im, double target_im, double extraction_window)
    {
      if (extraction_window <= 0.0 || observed_im < 0.0 || target_im < 0.0)
      {
        return true;
      }
      return std::fabs(observed_im - target_im) <= extraction_window * 0.5;
    }

    void keepTopPeaks_(std::vector<PeakEvidence>& peaks, Size top_n)
    {
      if (top_n > 0 && peaks.size() > top_n)
      {
        std::nth_element(peaks.begin(), peaks.begin() + static_cast<std::ptrdiff_t>(top_n), peaks.end(),
                         [](const PeakEvidence& lhs, const PeakEvidence& rhs)
                         {
                           return lhs.intensity > rhs.intensity;
                         });
        peaks.resize(top_n);
      }
      std::sort(peaks.begin(), peaks.end(),
                [](const PeakEvidence& lhs, const PeakEvidence& rhs)
                {
                  return lhs.mz < rhs.mz;
                });
    }

    std::vector<PeakEvidence> extractCentroidPeaks_(const OpenSwath::SpectrumPtr& spectrum, Size top_n)
    {
      std::vector<PeakEvidence> peaks;
      if (!spectrum)
      {
        return peaks;
      }

      OpenSwath::BinaryDataArrayPtr mz_array = spectrum->getMZArray();
      OpenSwath::BinaryDataArrayPtr intensity_array = spectrum->getIntensityArray();
      if (!mz_array || !intensity_array)
      {
        return peaks;
      }

      const Size n = std::min(mz_array->data.size(), intensity_array->data.size());
      peaks.reserve(top_n == 0 ? n : std::min(n, top_n));
      OpenSwath::BinaryDataArrayPtr im_array = spectrum->getDriftTimeArray();
      const bool has_im = im_array && im_array->data.size() >= n;
      for (Size i = 0; i < n; ++i)
      {
        const double intensity = intensity_array->data[i];
        if (intensity > 0.0 && std::isfinite(intensity) && std::isfinite(mz_array->data[i]))
        {
          PeakEvidence peak;
          peak.mz = mz_array->data[i];
          peak.intensity = intensity;
          peak.im = has_im ? im_array->data[i] : -1.0;
          peaks.push_back(peak);
        }
      }
      keepTopPeaks_(peaks, top_n);
      return peaks;
    }

    std::vector<PeakEvidence> extractPickedPeaks_(const OpenSwath::SpectrumPtr& spectrum,
                                                  Size top_n,
                                                  const Param& params,
                                                  bool use_gauss)
    {
      std::vector<PeakEvidence> peaks;
      if (!spectrum)
      {
        return peaks;
      }

      MSSpectrum openms_spectrum;
      OpenSwathDataAccessHelper::convertToOpenMSSpectrum(spectrum, openms_spectrum);
      if (openms_spectrum.empty())
      {
        return peaks;
      }
      openms_spectrum.sortByPosition();

      MSSpectrum smoothed_spectrum = openms_spectrum;
      if (use_gauss)
      {
        GaussFilter gauss;
        gauss.filter(smoothed_spectrum);
      }
      else
      {
        SavitzkyGolayFilter sgolay;
        sgolay.filter(smoothed_spectrum);
      }

      Param picker_params = PeakPickerHiRes().getDefaults();
      picker_params.update(params);
      picker_params.setValue("spacing_difference", 0.0);
      picker_params.setValue("spacing_difference_gap", 0.0);

      PeakPickerHiRes picker;
      picker.setParameters(picker_params);
      MSSpectrum picked_spectrum;
      picker.pick(smoothed_spectrum, picked_spectrum);

      peaks.reserve(top_n == 0 ? picked_spectrum.size() : std::min(static_cast<Size>(picked_spectrum.size()), top_n));
      for (const auto& peak : picked_spectrum)
      {
        if (peak.getIntensity() > 0.0 && std::isfinite(peak.getIntensity()) && std::isfinite(peak.getMZ()))
        {
          peaks.push_back({peak.getMZ(), peak.getIntensity(), -1.0});
        }
      }
      keepTopPeaks_(peaks, top_n);
      return peaks;
    }

    std::vector<PeakEvidence> extractPeaks_(const OpenSwath::SpectrumPtr& spectrum,
                                            Size top_n,
                                            bool peak_picking_enabled,
                                            const Param& picker_params,
                                            bool use_gauss)
    {
      if (peak_picking_enabled)
      {
        return extractPickedPeaks_(spectrum, top_n, picker_params, use_gauss);
      }
      return extractCentroidPeaks_(spectrum, top_n);
    }

    void updateMS1Evidence_(std::unordered_map<Size, EvidenceAccum>& local_evidence,
                            Size candidate_index,
                            double intensity,
                            double rt)
    {
      EvidenceAccum& evidence = local_evidence[candidate_index];
      ++evidence.ms1_hit_count;
      evidence.ms1_sum_intensity += intensity;
      if (intensity > evidence.ms1_max_intensity)
      {
        evidence.ms1_max_intensity = intensity;
        evidence.ms1_best_rt = rt;
      }
    }

    void updateMS2Evidence_(std::unordered_map<Size, EvidenceAccum>& local_evidence,
                            Size candidate_index,
                            Size distinct_fragment_hits,
                            double sum_intensity,
                            double max_intensity,
                            double rt)
    {
      EvidenceAccum& evidence = local_evidence[candidate_index];
      evidence.ms2_hit_count += distinct_fragment_hits;
      evidence.ms2_sum_intensity += sum_intensity;
      evidence.ms2_best_fragment_hits = std::max(evidence.ms2_best_fragment_hits, distinct_fragment_hits);
      if (max_intensity > evidence.ms2_max_intensity)
      {
        evidence.ms2_max_intensity = max_intensity;
        evidence.ms2_best_rt = rt;
      }
    }

    bool hasDecoyPrefix_(const std::string& id)
    {
      return id.find("DECOY") == 0 || id.find("Decoy") == 0 || id.find("decoy") == 0;
    }

    std::vector<PrecursorCandidate> buildCandidates_(const OpenSwath::LightTargetedExperiment& transition_exp,
                                                     Size top_transitions_per_precursor)
    {
      std::unordered_map<std::string, const OpenSwath::LightCompound*> compounds_by_id;
      compounds_by_id.reserve(transition_exp.compounds.size());
      for (const auto& compound : transition_exp.compounds)
      {
        compounds_by_id[compound.id] = &compound;
      }

      std::vector<PrecursorCandidate> candidates;
      std::unordered_map<std::string, Size> candidate_by_id;
      candidate_by_id.reserve(transition_exp.compounds.size());

      for (const auto& transition : transition_exp.transitions)
      {
        if (transition.getDecoy() || hasDecoyPrefix_(transition.getPeptideRef()))
        {
          continue;
        }

        const std::string compound_id = transition.getPeptideRef();
        if (compound_id.empty())
        {
          continue;
        }

        auto candidate_it = candidate_by_id.find(compound_id);
        if (candidate_it == candidate_by_id.end())
        {
          PrecursorCandidate candidate;
          candidate.compound_id = compound_id;
          const auto compound_it = compounds_by_id.find(compound_id);
          if (compound_it != compounds_by_id.end())
          {
            candidate.compound = *(compound_it->second);
          }
          else
          {
            candidate.compound.id = compound_id;
            candidate.compound.rt = -1.0;
            candidate.compound.sequence = compound_id;
          }
          candidate.precursor_mz = transition.getPrecursorMZ();
          candidate.precursor_im = candidate.compound.getDriftTime() >= 0.0 ? candidate.compound.getDriftTime() : transition.getPrecursorIM();
          candidates.push_back(std::move(candidate));
          const Size new_index = candidates.size() - 1;
          candidate_by_id[compound_id] = new_index;
          candidate_it = candidate_by_id.find(compound_id);
        }

        PrecursorCandidate& candidate = candidates[candidate_it->second];
        candidate.transitions.push_back(transition);
        if (candidate.precursor_mz <= 0.0)
        {
          candidate.precursor_mz = transition.getPrecursorMZ();
        }
        if (candidate.precursor_im < 0.0 && transition.getPrecursorIM() >= 0.0)
        {
          candidate.precursor_im = transition.getPrecursorIM();
        }

        if (transition.isDetectingTransition() || transition.isQuantifyingTransition())
        {
          FragmentCandidate fragment;
          fragment.product_mz = transition.getProductMZ();
          fragment.library_intensity = transition.getLibraryIntensity();
          fragment.transition_name = transition.getNativeID();
          candidate.fragments.push_back(std::move(fragment));
        }
      }

      for (auto& candidate : candidates)
      {
        std::sort(candidate.fragments.begin(), candidate.fragments.end(),
                  [](const FragmentCandidate& lhs, const FragmentCandidate& rhs)
                  {
                    if (lhs.library_intensity == rhs.library_intensity)
                    {
                      return lhs.transition_name < rhs.transition_name;
                    }
                    return lhs.library_intensity > rhs.library_intensity;
                  });
        if (top_transitions_per_precursor > 0 && candidate.fragments.size() > top_transitions_per_precursor)
        {
          candidate.fragments.resize(top_transitions_per_precursor);
        }
        for (Size i = 0; i < candidate.fragments.size(); ++i)
        {
          candidate.fragments[i].slot = i;
        }
      }

      return candidates;
    }

    std::vector<MzIndexEntry> buildPrecursorIndex_(const std::vector<PrecursorCandidate>& candidates)
    {
      std::vector<MzIndexEntry> index;
      index.reserve(candidates.size());
      for (Size i = 0; i < candidates.size(); ++i)
      {
        if (candidates[i].precursor_mz > 0.0)
        {
          index.push_back({candidates[i].precursor_mz, i});
        }
      }
      std::sort(index.begin(), index.end(),
                [](const MzIndexEntry& lhs, const MzIndexEntry& rhs)
                {
                  if (lhs.mz == rhs.mz)
                  {
                    return lhs.candidate_index < rhs.candidate_index;
                  }
                  return lhs.mz < rhs.mz;
                });
      return index;
    }

    bool candidateInSwathMap_(const PrecursorCandidate& candidate,
                              const OpenSwath::SwathMap& map,
                              const ChromExtractParams& params,
                              bool pasef)
    {
      if (!(map.lower < candidate.precursor_mz && candidate.precursor_mz < map.upper))
      {
        return false;
      }
      if (std::fabs(map.upper - candidate.precursor_mz) < params.min_upper_edge_dist)
      {
        return false;
      }
      if (pasef)
      {
        return candidate.precursor_im >= 0.0 &&
               map.imLower < candidate.precursor_im &&
               candidate.precursor_im < map.imUpper;
      }
      return true;
    }

    Size countPasefWindowMatches_(const std::vector<PrecursorCandidate>& candidates,
                                  const std::vector<OpenSwath::SwathMap>& swath_maps,
                                  const ChromExtractParams& params,
                                  const PrecursorIMTransform& transform)
    {
      Size matches = 0;
      for (const auto& candidate : candidates)
      {
        if (candidate.precursor_mz <= 0.0 || candidate.precursor_im < 0.0)
        {
          continue;
        }
        const double scaled_im = candidate.precursor_im * transform.factor(candidate);
        for (const auto& map : swath_maps)
        {
          if (map.ms1 || map.imLower < 0.0 || map.imUpper < 0.0)
          {
            continue;
          }
          if (map.lower < candidate.precursor_mz &&
              candidate.precursor_mz < map.upper &&
              std::fabs(map.upper - candidate.precursor_mz) >= params.min_upper_edge_dist &&
              map.imLower < scaled_im &&
              scaled_im < map.imUpper)
          {
            ++matches;
            break;
          }
        }
      }
      return matches;
    }

    PrecursorIMTransform determinePrecursorIMTransform_(const std::vector<PrecursorCandidate>& candidates,
                                                        const std::vector<OpenSwath::SwathMap>& swath_maps,
                                                        const ChromExtractParams& params,
                                                        bool pasef)
    {
      if (!pasef)
      {
        return {};
      }

      const bool has_pasef_windows = std::any_of(swath_maps.begin(), swath_maps.end(),
                                                 [](const OpenSwath::SwathMap& map)
                                                 {
                                                   return !map.ms1 && map.imLower >= 0.0 && map.imUpper >= 0.0;
                                                 });
      if (!has_pasef_windows)
      {
        return {};
      }

      const double scales[] = {1.0, 1000.0, 0.001};
      PrecursorIMTransform best_transform;
      Size best_matches = 0;
      Size unscaled_matches = 0;
      for (double scale : scales)
      {
        for (bool multiply_by_charge : {false, true})
        {
          PrecursorIMTransform transform{scale, multiply_by_charge};
          const Size matches = countPasefWindowMatches_(candidates, swath_maps, params, transform);
          if (transform.isIdentity())
          {
            unscaled_matches = matches;
          }
          if (matches > best_matches)
          {
            best_matches = matches;
            best_transform = transform;
          }
        }
      }

      if (!best_transform.isIdentity() && best_matches > 0 &&
          (unscaled_matches == 0 || best_matches > unscaled_matches * 10))
      {
        return best_transform;
      }
      return {};
    }

    void applyPrecursorIMTransform_(std::vector<PrecursorCandidate>& candidates,
                                    const PrecursorIMTransform& transform)
    {
      if (transform.isIdentity())
      {
        return;
      }

      for (auto& candidate : candidates)
      {
        const double factor = transform.factor(candidate);
        if (candidate.precursor_im >= 0.0)
        {
          candidate.precursor_im *= factor;
          candidate.compound.drift_time = candidate.precursor_im;
        }
        for (auto& transition : candidate.transitions)
        {
          if (transition.precursor_im >= 0.0)
          {
            transition.precursor_im *= factor;
          }
          else if (candidate.precursor_im >= 0.0)
          {
            transition.precursor_im = candidate.precursor_im;
          }
        }
      }
    }

    std::vector<ProductIndexEntry> buildProductIndexForMap_(const std::vector<PrecursorCandidate>& candidates,
                                                            const std::vector<MzIndexEntry>& precursor_index,
                                                            const OpenSwath::SwathMap& map,
                                                            const ChromExtractParams& params,
                                                            bool pasef)
    {
      std::vector<ProductIndexEntry> index;
      const auto lower_it = std::lower_bound(precursor_index.begin(), precursor_index.end(), map.lower,
                                             [](const MzIndexEntry& entry, double value)
                                             {
                                               return entry.mz < value;
                                             });
      const auto upper_it = std::lower_bound(precursor_index.begin(), precursor_index.end(), map.upper,
                                             [](const MzIndexEntry& entry, double value)
                                             {
                                               return entry.mz < value;
                                             });
      for (auto precursor_it = lower_it; precursor_it != upper_it; ++precursor_it)
      {
        const Size i = precursor_it->candidate_index;
        const auto& candidate = candidates[i];
        if (!candidateInSwathMap_(candidate, map, params, pasef))
        {
          continue;
        }
        for (const auto& fragment : candidate.fragments)
        {
          if (fragment.product_mz > 0.0 && fragment.slot < 64)
          {
            index.push_back({fragment.product_mz, i, fragment.slot});
          }
        }
      }
      std::sort(index.begin(), index.end(),
                [](const ProductIndexEntry& lhs, const ProductIndexEntry& rhs)
                {
                  if (lhs.mz == rhs.mz)
                  {
                    if (lhs.candidate_index == rhs.candidate_index)
                    {
                      return lhs.fragment_slot < rhs.fragment_slot;
                    }
                    return lhs.candidate_index < rhs.candidate_index;
                  }
                  return lhs.mz < rhs.mz;
                });
      return index;
    }

    void mergeLocalEvidence_(std::vector<EvidenceAccum>& merged,
                             const std::vector<std::unordered_map<Size, EvidenceAccum>>& local_results)
    {
      for (const auto& local : local_results)
      {
        for (const auto& item : local)
        {
          merged[item.first].merge(item.second);
        }
      }
    }

    void scanMS1MapSerial_(const OpenSwath::SwathMap& map,
                           const std::vector<PrecursorCandidate>& candidates,
                           const std::vector<MzIndexEntry>& precursor_index,
                           const ChromExtractParams& params,
                           Size top_peaks_per_spectrum,
                           bool peak_picking_enabled,
                           const Param& picker_params,
                           bool peak_picking_use_gauss,
                           std::unordered_map<Size, EvidenceAccum>& local_evidence)
    {
      if (!map.sptr || precursor_index.empty() || map.sptr->getNrSpectra() == 0)
      {
        return;
      }

      OpenSwath::SpectrumAccessPtr access = map.sptr->lightClone();
      const Size nr_spectra = map.sptr->getNrSpectra();
      for (Size spectrum_index = 0; spectrum_index < nr_spectra; ++spectrum_index)
      {
        OpenSwath::SpectrumPtr spectrum = access->getSpectrumById(static_cast<int>(spectrum_index));
        const double rt = access->getSpectrumMetaById(static_cast<int>(spectrum_index)).RT;
        const std::vector<PeakEvidence> peaks = extractPeaks_(
          spectrum, top_peaks_per_spectrum, peak_picking_enabled, picker_params, peak_picking_use_gauss);
        for (const auto& peak : peaks)
        {
          const double abs_window = halfWindow_(peak.mz, params);
          const auto lower_it = std::lower_bound(precursor_index.begin(), precursor_index.end(), peak.mz - abs_window,
                                                 [](const MzIndexEntry& entry, double value)
                                                 {
                                                   return entry.mz < value;
                                                 });
          const auto upper_it = std::upper_bound(precursor_index.begin(), precursor_index.end(), peak.mz + abs_window,
                                                 [](double value, const MzIndexEntry& entry)
                                                 {
                                                   return value < entry.mz;
                                                 });
          for (auto it = lower_it; it != upper_it; ++it)
          {
            const auto& candidate = candidates[it->candidate_index];
            if (withinMzTolerance_(peak.mz, it->mz, params) &&
                withinIMTolerance_(peak.im, candidate.precursor_im, params.im_extraction_window))
            {
              updateMS1Evidence_(local_evidence, it->candidate_index, peak.intensity, rt);
            }
          }
        }
      }
    }

    void scanMS2MapSerial_(const OpenSwath::SwathMap& map,
                           const std::vector<ProductIndexEntry>& product_index,
                           const ChromExtractParams& params,
                           Size top_peaks_per_spectrum,
                           Size min_fragment_hits,
                           bool peak_picking_enabled,
                           const Param& picker_params,
                           bool peak_picking_use_gauss,
                           std::unordered_map<Size, EvidenceAccum>& local_evidence)
    {
      if (!map.sptr || product_index.empty() || map.sptr->getNrSpectra() == 0)
      {
        return;
      }

      OpenSwath::SpectrumAccessPtr access = map.sptr->lightClone();
      const Size nr_spectra = map.sptr->getNrSpectra();
      for (Size spectrum_index = 0; spectrum_index < nr_spectra; ++spectrum_index)
      {
        OpenSwath::SpectrumPtr spectrum = access->getSpectrumById(static_cast<int>(spectrum_index));
        const double rt = access->getSpectrumMetaById(static_cast<int>(spectrum_index)).RT;
        const std::vector<PeakEvidence> peaks = extractPeaks_(
          spectrum, top_peaks_per_spectrum, peak_picking_enabled, picker_params, peak_picking_use_gauss);
        std::unordered_map<Size, SpectrumMS2Hit> spectrum_hits;
        for (const auto& peak : peaks)
        {
          const double abs_window = halfWindow_(peak.mz, params);
          const auto lower_it = std::lower_bound(product_index.begin(), product_index.end(), peak.mz - abs_window,
                                                 [](const ProductIndexEntry& entry, double value)
                                                 {
                                                   return entry.mz < value;
                                                 });
          const auto upper_it = std::upper_bound(product_index.begin(), product_index.end(), peak.mz + abs_window,
                                                 [](double value, const ProductIndexEntry& entry)
                                                 {
                                                   return value < entry.mz;
                                                 });
          for (auto it = lower_it; it != upper_it; ++it)
          {
            if (withinMzTolerance_(peak.mz, it->mz, params))
            {
              SpectrumMS2Hit& hit = spectrum_hits[it->candidate_index];
              hit.fragment_mask |= (std::uint64_t{1} << it->fragment_slot);
              hit.sum_intensity += peak.intensity;
              hit.max_intensity = std::max(hit.max_intensity, peak.intensity);
            }
          }
        }
        for (const auto& hit : spectrum_hits)
        {
          const Size distinct_fragments = countBits_(hit.second.fragment_mask);
          if (distinct_fragments >= min_fragment_hits)
          {
            updateMS2Evidence_(local_evidence, hit.first, distinct_fragments, hit.second.sum_intensity, hit.second.max_intensity, rt);
          }
        }
      }
    }

    void scanMS1Map_(const OpenSwath::SwathMap& map,
                     const std::vector<PrecursorCandidate>& candidates,
                     const std::vector<MzIndexEntry>& precursor_index,
                     const ChromExtractParams& params,
                     Size top_peaks_per_spectrum,
                     bool peak_picking_enabled,
                     const Param& picker_params,
                     bool peak_picking_use_gauss,
                     int threads,
                     std::vector<EvidenceAccum>& merged_evidence)
    {
      if (!map.sptr || precursor_index.empty() || map.sptr->getNrSpectra() == 0)
      {
        return;
      }

      const int thread_count = resolveThreadCount_(threads);
      const Size nr_spectra = map.sptr->getNrSpectra();
      std::vector<std::unordered_map<Size, EvidenceAccum>> local_results(static_cast<Size>(thread_count));

#ifdef _OPENMP
      if (thread_count > 1 && !omp_in_parallel())
      {
#pragma omp parallel num_threads(thread_count)
        {
          const int thread_id = omp_get_thread_num();
          OpenSwath::SpectrumAccessPtr access = map.sptr->lightClone();
          auto& local_evidence = local_results[static_cast<Size>(thread_id)];
#pragma omp for schedule(dynamic, 16)
          for (SignedSize spectrum_index = 0; spectrum_index < static_cast<SignedSize>(nr_spectra); ++spectrum_index)
          {
            OpenSwath::SpectrumPtr spectrum = access->getSpectrumById(static_cast<int>(spectrum_index));
            const double rt = access->getSpectrumMetaById(static_cast<int>(spectrum_index)).RT;
            const std::vector<PeakEvidence> peaks = extractPeaks_(spectrum, top_peaks_per_spectrum, peak_picking_enabled, picker_params, peak_picking_use_gauss);
            for (const auto& peak : peaks)
            {
              const double abs_window = halfWindow_(peak.mz, params);
              const auto lower_it = std::lower_bound(precursor_index.begin(), precursor_index.end(), peak.mz - abs_window,
                                                     [](const MzIndexEntry& entry, double value)
                                                     {
                                                       return entry.mz < value;
                                                     });
              const auto upper_it = std::upper_bound(precursor_index.begin(), precursor_index.end(), peak.mz + abs_window,
                                                     [](double value, const MzIndexEntry& entry)
                                                     {
                                                       return value < entry.mz;
                                                     });
              for (auto it = lower_it; it != upper_it; ++it)
              {
                const auto& candidate = candidates[it->candidate_index];
                if (withinMzTolerance_(peak.mz, it->mz, params) &&
                    withinIMTolerance_(peak.im, candidate.precursor_im, params.im_extraction_window))
                {
                  updateMS1Evidence_(local_evidence, it->candidate_index, peak.intensity, rt);
                }
              }
            }
          }
        }
      }
      else
#endif
      {
        scanMS1MapSerial_(map, candidates, precursor_index, params, top_peaks_per_spectrum,
                          peak_picking_enabled, picker_params, peak_picking_use_gauss, local_results[0]);
      }

      mergeLocalEvidence_(merged_evidence, local_results);
    }

    void scanMS2Map_(const OpenSwath::SwathMap& map,
                     const std::vector<ProductIndexEntry>& product_index,
                     const ChromExtractParams& params,
                     Size top_peaks_per_spectrum,
                     Size min_fragment_hits,
                     bool peak_picking_enabled,
                     const Param& picker_params,
                     bool peak_picking_use_gauss,
                     int threads,
                     std::vector<EvidenceAccum>& merged_evidence)
    {
      if (!map.sptr || product_index.empty() || map.sptr->getNrSpectra() == 0)
      {
        return;
      }

      const int thread_count = resolveThreadCount_(threads);
      const Size nr_spectra = map.sptr->getNrSpectra();
      std::vector<std::unordered_map<Size, EvidenceAccum>> local_results(static_cast<Size>(thread_count));

#ifdef _OPENMP
      if (thread_count > 1 && !omp_in_parallel())
      {
#pragma omp parallel num_threads(thread_count)
        {
          const int thread_id = omp_get_thread_num();
          OpenSwath::SpectrumAccessPtr access = map.sptr->lightClone();
          auto& local_evidence = local_results[static_cast<Size>(thread_id)];
#pragma omp for schedule(dynamic, 16)
          for (SignedSize spectrum_index = 0; spectrum_index < static_cast<SignedSize>(nr_spectra); ++spectrum_index)
          {
            OpenSwath::SpectrumPtr spectrum = access->getSpectrumById(static_cast<int>(spectrum_index));
            const double rt = access->getSpectrumMetaById(static_cast<int>(spectrum_index)).RT;
            const std::vector<PeakEvidence> peaks = extractPeaks_(spectrum, top_peaks_per_spectrum, peak_picking_enabled, picker_params, peak_picking_use_gauss);
            std::unordered_map<Size, SpectrumMS2Hit> spectrum_hits;
            for (const auto& peak : peaks)
            {
              const double abs_window = halfWindow_(peak.mz, params);
              const auto lower_it = std::lower_bound(product_index.begin(), product_index.end(), peak.mz - abs_window,
                                                     [](const ProductIndexEntry& entry, double value)
                                                     {
                                                       return entry.mz < value;
                                                     });
              const auto upper_it = std::upper_bound(product_index.begin(), product_index.end(), peak.mz + abs_window,
                                                     [](double value, const ProductIndexEntry& entry)
                                                     {
                                                       return value < entry.mz;
                                                     });
              for (auto it = lower_it; it != upper_it; ++it)
              {
                if (withinMzTolerance_(peak.mz, it->mz, params))
                {
                  SpectrumMS2Hit& hit = spectrum_hits[it->candidate_index];
                  hit.fragment_mask |= (std::uint64_t{1} << it->fragment_slot);
                  hit.sum_intensity += peak.intensity;
                  hit.max_intensity = std::max(hit.max_intensity, peak.intensity);
                }
              }
            }
            for (const auto& hit : spectrum_hits)
            {
              const Size distinct_fragments = countBits_(hit.second.fragment_mask);
              if (distinct_fragments >= min_fragment_hits)
              {
                updateMS2Evidence_(local_evidence, hit.first, distinct_fragments, hit.second.sum_intensity, hit.second.max_intensity, rt);
              }
            }
          }
        }
      }
      else
#endif
      {
        scanMS2MapSerial_(map, product_index, params, top_peaks_per_spectrum, min_fragment_hits,
                          peak_picking_enabled, picker_params, peak_picking_use_gauss, local_results[0]);
      }

      mergeLocalEvidence_(merged_evidence, local_results);
    }

    void scanMS2MapsWithWaveScheduler_(const std::vector<OpenSwath::SwathMap>& swath_maps,
                                       const std::vector<PrecursorCandidate>& candidates,
                                       const std::vector<MzIndexEntry>& precursor_index,
                                       const ChromExtractParams& params,
                                       Size top_peaks_per_spectrum,
                                       Size min_fragment_hits,
                                       bool peak_picking_enabled,
                                       const Param& picker_params,
                                       bool peak_picking_use_gauss,
                                       int threads,
                                       bool pasef,
                                       std::vector<EvidenceAccum>& merged_evidence)
    {
      const int thread_count = resolveThreadCount_(threads);
      if (thread_count <= 1 || swath_maps.size() <= 1)
      {
        for (const auto& map : swath_maps)
        {
          const std::vector<ProductIndexEntry> product_index =
            buildProductIndexForMap_(candidates, precursor_index, map, params, pasef);
          scanMS2Map_(map, product_index, params, top_peaks_per_spectrum, min_fragment_hits,
                      peak_picking_enabled, picker_params, peak_picking_use_gauss, 1, merged_evidence);
        }
        return;
      }

      OpenSwathWorkflowScheduler::Options scheduler_options;
      scheduler_options.osw_buffer_bytes = 0ULL;
      scheduler_options.bytes_per_spectrum = 0ULL;
      // TransitionListEvidenceFilter streams one spectrum at a time, so only
      // a small per-SWATH working set must be kept active concurrently.
      scheduler_options.per_swath_overhead_bytes = 8ULL * 1024ULL * 1024ULL;
      scheduler_options.max_concurrent_swaths = thread_count;
      scheduler_options.scoring_threads = static_cast<Size>(thread_count);

      const OpenSwathWorkflowScheduler::ConcurrencyEstimate concurrency_estimate =
        OpenSwathWorkflowScheduler::estimateConcurrency(swath_maps, scheduler_options);
      const std::vector<OpenSwathWorkflowScheduler::Wave> swath_waves =
        OpenSwathWorkflowScheduler::planWaves(swath_maps, scheduler_options, concurrency_estimate);

      OPENMS_LOG_DEBUG << "TransitionListEvidenceFilter uses SWATH wave scheduler with "
                       << swath_waves.size() << " waves and max_concurrent_swaths="
                       << concurrency_estimate.max_concurrent_swaths << ".\n";

      std::vector<std::unordered_map<Size, EvidenceAccum>> local_results(static_cast<Size>(thread_count));
      for (const auto& wave : swath_waves)
      {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) num_threads(thread_count)
#endif
        for (SignedSize wave_pos = 0; wave_pos < static_cast<SignedSize>(wave.swath_indices.size()); ++wave_pos)
        {
          const OpenSwath::SwathMap& map = swath_maps[wave.swath_indices[static_cast<Size>(wave_pos)]];
          const std::vector<ProductIndexEntry> product_index =
            buildProductIndexForMap_(candidates, precursor_index, map, params, pasef);
#ifdef _OPENMP
          auto& local_evidence = local_results[static_cast<Size>(omp_get_thread_num())];
#else
          auto& local_evidence = local_results[0];
#endif
          scanMS2MapSerial_(map, product_index, params, top_peaks_per_spectrum, min_fragment_hits,
                            peak_picking_enabled, picker_params, peak_picking_use_gauss, local_evidence);
        }
      }

      mergeLocalEvidence_(merged_evidence, local_results);
    }

    OpenSwath::LightTargetedExperiment buildFilteredExperiment_(const OpenSwath::LightTargetedExperiment& transition_exp,
                                                                const std::vector<PrecursorCandidate>& candidates,
                                                                const std::vector<bool>& supported)
    {
      OpenSwath::LightTargetedExperiment filtered_exp;
      std::unordered_set<std::string> matching_proteins;

      for (Size i = 0; i < candidates.size(); ++i)
      {
        if (supported[i])
        {
          filtered_exp.compounds.push_back(candidates[i].compound);
          for (const auto& protein_ref : candidates[i].compound.protein_refs)
          {
            matching_proteins.insert(protein_ref);
          }
        }
      }

      for (Size i = 0; i < candidates.size(); ++i)
      {
        if (supported[i])
        {
          for (const auto& transition : candidates[i].transitions)
          {
            filtered_exp.transitions.push_back(transition);
          }
        }
      }

      for (const auto& protein : transition_exp.proteins)
      {
        if (matching_proteins.find(protein.id) != matching_proteins.end())
        {
          filtered_exp.proteins.push_back(protein);
        }
      }

      return filtered_exp;
    }
  }

  TransitionListEvidenceFilter::TransitionListEvidenceFilter() :
    DefaultParamHandler("TransitionListEvidenceFilter"),
    ProgressLogger()
  {
    defaults_.setValue("enabled", "true", "Enable raw-data evidence prefiltering.");
    defaults_.setValidStrings("enabled", {"true", "false"});

    defaults_.setValue("evidence_sources", "hybrid", "Evidence source used to keep precursor candidates: 'ms1', 'ms2', or 'hybrid' (MS1 OR MS2).");
    defaults_.setValidStrings("evidence_sources", {"ms1", "ms2", "hybrid"});

    defaults_.setValue("ms1_top_peaks_per_spectrum", 1000, "Number of most intense MS1 peaks to keep per spectrum before precursor m/z matching.");
    defaults_.setMinInt("ms1_top_peaks_per_spectrum", 1);

    defaults_.setValue("ms2_top_peaks_per_spectrum", 1000, "Number of most intense MS2 peaks to keep per spectrum before fragment m/z matching.");
    defaults_.setMinInt("ms2_top_peaks_per_spectrum", 1);

    defaults_.setValue("ms2_top_transitions_per_precursor", 6, "Number of highest library-intensity fragment transitions to index per precursor for MS2 evidence.");
    defaults_.setMinInt("ms2_top_transitions_per_precursor", 1);
    defaults_.setMaxInt("ms2_top_transitions_per_precursor", 64);

    defaults_.setValue("ms2_min_fragment_hits", 4, "Minimum number of distinct indexed fragments that must match in one MS2 spectrum to support a precursor.");
    defaults_.setMinInt("ms2_min_fragment_hits", 1);
    defaults_.setMaxInt("ms2_min_fragment_hits", 64);

    defaults_.setValue("min_supported_precursors", 3, "Minimum number of supported precursors required when the prefilter is enabled.");
    defaults_.setMinInt("min_supported_precursors", 1);

    defaults_.setValue("peak_picking:enabled", "false", "Run PeakPickerHiRes on spectra before top-N matching. This is intended for profile data and is disabled by default.");
    defaults_.setValidStrings("peak_picking:enabled", {"true", "false"});

    defaults_.setValue("peak_picking:use_gauss", "true", "Use Gaussian smoothing before PeakPickerHiRes. If false, use Savitzky-Golay smoothing.");
    defaults_.setValidStrings("peak_picking:use_gauss", {"true", "false"});

    defaults_.insert("peak_picking:PeakPickerHiRes:", PeakPickerHiRes().getDefaults());
    subsections_.emplace_back("peak_picking:PeakPickerHiRes");

    defaultsToParam_();
  }

  void TransitionListEvidenceFilter::updateMembers_()
  {
    enabled_ = param_.getValue("enabled").toBool();
    evidence_sources_ = param_.getValue("evidence_sources").toString();
    ms1_top_peaks_per_spectrum_ = static_cast<Size>(param_.getValue("ms1_top_peaks_per_spectrum"));
    ms2_top_peaks_per_spectrum_ = static_cast<Size>(param_.getValue("ms2_top_peaks_per_spectrum"));
    ms2_top_transitions_per_precursor_ = static_cast<Size>(param_.getValue("ms2_top_transitions_per_precursor"));
    ms2_min_fragment_hits_ = static_cast<Size>(param_.getValue("ms2_min_fragment_hits"));
    min_supported_precursors_ = static_cast<Size>(param_.getValue("min_supported_precursors"));
    peak_picking_enabled_ = param_.getValue("peak_picking:enabled").toBool();
    peak_picking_use_gauss_ = param_.getValue("peak_picking:use_gauss").toBool();
  }

  TransitionListEvidenceFilter::Result TransitionListEvidenceFilter::filter(
    const std::vector<OpenSwath::SwathMap>& swath_maps,
    const OpenSwath::LightTargetedExperiment& transition_exp,
    const ChromExtractParams& ms1_params,
    const ChromExtractParams& ms2_params,
    bool pasef,
    int threads) const
  {
    Result result;

    const bool use_ms1 = evidence_sources_ == "ms1" || evidence_sources_ == "hybrid";
    const bool use_ms2 = evidence_sources_ == "ms2" || evidence_sources_ == "hybrid";
    if (!use_ms1 && !use_ms2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "TransitionListEvidenceFilter evidence_sources must be 'ms1', 'ms2', or 'hybrid'.");
    }
    if (evidence_sources_ == "ms1" && ms1_params.mz_extraction_window <= 0.0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "TransitionListEvidenceFilter requires a positive MS1 m/z extraction window for MS1 evidence.");
    }
    if (evidence_sources_ == "ms2" && ms2_params.mz_extraction_window <= 0.0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "TransitionListEvidenceFilter requires a positive MS2/iRT m/z extraction window for MS2 evidence.");
    }

    startProgress(0, 1, "Filtering Transition List");
    bool progress_started = true;
    try
    {
      std::vector<PrecursorCandidate> candidates = buildCandidates_(transition_exp, ms2_top_transitions_per_precursor_);
      const PrecursorIMTransform precursor_im_transform = determinePrecursorIMTransform_(candidates, swath_maps, ms2_params, pasef);
      result.precursor_im_scale = precursor_im_transform.scale;
      result.precursor_im_scaled_by_charge = precursor_im_transform.multiply_by_charge;
      if (!precursor_im_transform.isIdentity())
      {
        OPENMS_LOG_WARN << "TransitionListEvidenceFilter detected a precursor ion mobility transform mismatch between the transition library and PASEF windows. "
                        << "Applying scale factor " << result.precursor_im_scale;
        if (result.precursor_im_scaled_by_charge)
        {
          OPENMS_LOG_WARN << " and multiplying by precursor charge";
        }
        OPENMS_LOG_WARN << " to precursor ion mobility values for matching and filtered output.\n";
        applyPrecursorIMTransform_(candidates, precursor_im_transform);
      }
      result.total_target_precursors = candidates.size();
      result.evidence.resize(candidates.size());
      for (Size i = 0; i < candidates.size(); ++i)
      {
        result.evidence[i].compound_id = candidates[i].compound_id;
        result.evidence[i].sequence = candidates[i].compound.sequence;
        result.evidence[i].precursor_mz = candidates[i].precursor_mz;
        result.evidence[i].precursor_im = candidates[i].precursor_im;
      }

      if (candidates.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "TransitionListEvidenceFilter found no non-decoy target precursors in the transition experiment.");
      }

      const Param picker_params = param_.copy("peak_picking:PeakPickerHiRes:", true);
      std::vector<EvidenceAccum> evidence(candidates.size());

      bool has_ms1_map = false;
      bool has_ms2_map = false;

      const std::vector<MzIndexEntry> precursor_index = buildPrecursorIndex_(candidates);

      if (use_ms1 && ms1_params.mz_extraction_window > 0.0)
      {
        for (const auto& map : swath_maps)
        {
          if (map.ms1 && map.sptr && map.sptr->getNrSpectra() > 0)
          {
            has_ms1_map = true;
            scanMS1Map_(map, candidates, precursor_index, ms1_params,
                        ms1_top_peaks_per_spectrum_, peak_picking_enabled_,
                        picker_params, peak_picking_use_gauss_, threads, evidence);
          }
        }
      }

      std::vector<OpenSwath::SwathMap> active_ms2_maps;
      if (use_ms2 && ms2_params.mz_extraction_window > 0.0)
      {
        for (const auto& map : swath_maps)
        {
          if (!map.ms1 && map.sptr && map.sptr->getNrSpectra() > 0)
          {
            has_ms2_map = true;
            active_ms2_maps.push_back(map);
          }
        }

        scanMS2MapsWithWaveScheduler_(active_ms2_maps, candidates, precursor_index, ms2_params,
                                      ms2_top_peaks_per_spectrum_, ms2_min_fragment_hits_,
                                      peak_picking_enabled_, picker_params,
                                      peak_picking_use_gauss_, threads, pasef, evidence);
      }

      if (evidence_sources_ == "ms1" && !has_ms1_map)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "TransitionListEvidenceFilter was configured for MS1 evidence but no MS1 spectra were available.");
      }
      if (evidence_sources_ == "ms2" && !has_ms2_map)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "TransitionListEvidenceFilter was configured for MS2 evidence but no MS2 spectra were available.");
      }
      if (evidence_sources_ == "hybrid" && !has_ms1_map && !has_ms2_map)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "TransitionListEvidenceFilter was configured for hybrid evidence but no usable MS1 or MS2 spectra were available.");
      }

      std::vector<bool> supported(candidates.size(), false);
      for (Size i = 0; i < candidates.size(); ++i)
      {
        auto& out = result.evidence[i];
        out.supported_ms1 = evidence[i].ms1_hit_count > 0;
        out.supported_ms2 = evidence[i].ms2_best_fragment_hits >= ms2_min_fragment_hits_;
        out.ms1_hit_count = evidence[i].ms1_hit_count;
        out.ms1_max_intensity = evidence[i].ms1_max_intensity;
        out.ms1_sum_intensity = evidence[i].ms1_sum_intensity;
        out.ms1_best_rt = evidence[i].ms1_best_rt;
        out.ms2_hit_count = evidence[i].ms2_hit_count;
        out.ms2_best_fragment_hits = evidence[i].ms2_best_fragment_hits;
        out.ms2_max_intensity = evidence[i].ms2_max_intensity;
        out.ms2_sum_intensity = evidence[i].ms2_sum_intensity;
        out.ms2_best_rt = evidence[i].ms2_best_rt;

        if (out.supported_ms1)
        {
          ++result.ms1_supported;
        }
        if (out.supported_ms2)
        {
          ++result.ms2_supported;
        }
        if (out.supported_ms1 && out.supported_ms2)
        {
          ++result.hybrid_supported;
        }

        if (evidence_sources_ == "ms1")
        {
          supported[i] = out.supported_ms1;
        }
        else if (evidence_sources_ == "ms2")
        {
          supported[i] = out.supported_ms2;
        }
        else
        {
          supported[i] = out.supported_ms1 || out.supported_ms2;
        }
        if (supported[i])
        {
          ++result.supported_precursors;
        }
      }

      if (enabled_ && result.supported_precursors < min_supported_precursors_)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "TransitionListEvidenceFilter retained only " + String(result.supported_precursors) +
                                         " supported target precursors, fewer than min_supported_precursors=" + String(min_supported_precursors_) +
                                         ". Disable Calibration:auto_irt:prefilter or loosen its parameters.");
      }

      result.filtered_targets = buildFilteredExperiment_(transition_exp, candidates, supported);

      std::ostringstream summary;
      summary << "TransitionListEvidenceFilter retained " << result.supported_precursors
              << " of " << result.total_target_precursors << " target precursors"
              << " (MS1: " << result.ms1_supported
              << ", MS2: " << result.ms2_supported
              << ", both: " << result.hybrid_supported
              << ", transitions: " << result.filtered_targets.transitions.size()
              << ", compounds: " << result.filtered_targets.compounds.size()
              << ").";
      result.summary = summary.str();

      OPENMS_LOG_INFO << result.summary << std::endl;
      endProgress();
      progress_started = false;
      return result;
    }
    catch (...)
    {
      if (progress_started)
      {
        endProgress();
      }
      throw;
    }
  }
}
