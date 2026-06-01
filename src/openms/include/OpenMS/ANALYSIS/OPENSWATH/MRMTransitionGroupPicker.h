// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

// Cross-correlation
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>

#include <numeric>

//#define DEBUG_TRANSITIONGROUPPICKER

namespace OpenMS
{

  // Shared thread-local pools used by MRMTransitionGroupPicker to avoid
  // repeated allocations across nested helpers (pickTransitionGroup, computeQuality_)
  struct PickerPoolShared
  {
    std::vector<MSChromatogram> picked_chroms;
    std::vector<MSChromatogram> smoothed_chroms;
    std::vector<const MSChromatogram*> picked_input_chromatograms;
    std::vector<const MSChromatogram*> transition_chromatograms;
    std::vector<MRMFeature> features;
    std::vector<double> left_edges;
    std::vector<double> right_edges;
    std::vector<std::vector<double>> all_ints;
    std::vector<double> mean_shapes;
    std::vector<double> mean_coel;
    std::vector<double> left_borders;
    std::vector<double> right_borders;
    std::vector<int> largest_peak_indices;
    std::vector<double> transition_total_xics;

    void reset()
    {
      picked_chroms.clear();
      smoothed_chroms.clear();
      picked_input_chromatograms.clear();
      transition_chromatograms.clear();
      features.clear();
      left_edges.clear();
      right_edges.clear();
      all_ints.clear();
      mean_shapes.clear();
      mean_coel.clear();
      left_borders.clear();
      right_borders.clear();
      largest_peak_indices.clear();
      transition_total_xics.clear();
    }
  };

  // Per-thread instance of the shared pools. Wrapped in a function so the
  // function-local `thread_local static` is the single point of definition
  // across all translation units / DLLs and we don't have to deal with
  // platform-specific behavior of namespace-scope `inline thread_local`
  // (duplicate TLS init wrapper on macOS clang) or `extern thread_local`
  // exported from a DLL (MSVC LNK2019 on TOPP-tool consumers).
  inline PickerPoolShared& g_picker_pool_shared()
  {
    thread_local PickerPoolShared instance;
    return instance;
  }

  namespace MRMTransitionGroupPickerMeta
  {
    inline UInt areaBackgroundLevel()
    {
      static const UInt id = MetaInfo::registry().registerName("area_background_level");
      return id;
    }

    inline UInt asymmetryFactor()
    {
      static const UInt id = MetaInfo::registry().registerName("asymmetry_factor");
      return id;
    }

    inline UInt baselineDelta2Height()
    {
      static const UInt id = MetaInfo::registry().registerName("baseline_delta_2_height");
      return id;
    }

    inline UInt endPositionAt5()
    {
      static const UInt id = MetaInfo::registry().registerName("end_position_at_5");
      return id;
    }

    inline UInt endPositionAt10()
    {
      static const UInt id = MetaInfo::registry().registerName("end_position_at_10");
      return id;
    }

    inline UInt endPositionAt50()
    {
      static const UInt id = MetaInfo::registry().registerName("end_position_at_50");
      return id;
    }

    inline UInt initialPeakQuality()
    {
      static const UInt id = MetaInfo::registry().registerName("initialPeakQuality");
      return id;
    }

    inline UInt leftWidth()
    {
      static const UInt id = MetaInfo::registry().registerName("leftWidth");
      return id;
    }

    inline UInt nativeID()
    {
      static const UInt id = MetaInfo::registry().registerName("native_id");
      return id;
    }

    inline UInt noiseBackgroundLevel()
    {
      static const UInt id = MetaInfo::registry().registerName("noise_background_level");
      return id;
    }

    inline UInt peakApexInt()
    {
      static const UInt id = MetaInfo::registry().registerName("peak_apex_int");
      return id;
    }

    inline UInt peakApexPosition()
    {
      static const UInt id = MetaInfo::registry().registerName("peak_apex_position");
      return id;
    }

    inline UInt peakApicesSum()
    {
      static const UInt id = MetaInfo::registry().registerName("peak_apices_sum");
      return id;
    }

    inline UInt peptideRef()
    {
      static const UInt id = MetaInfo::registry().registerName("PeptideRef");
      return id;
    }

    inline UInt pointsAcrossBaseline()
    {
      static const UInt id = MetaInfo::registry().registerName("points_across_baseline");
      return id;
    }

    inline UInt pointsAcrossHalfHeight()
    {
      static const UInt id = MetaInfo::registry().registerName("points_across_half_height");
      return id;
    }

    inline UInt potentialOutlier()
    {
      static const UInt id = MetaInfo::registry().registerName("potentialOutlier");
      return id;
    }

    inline UInt rightWidth()
    {
      static const UInt id = MetaInfo::registry().registerName("rightWidth");
      return id;
    }

    inline UInt slopeOfBaseline()
    {
      static const UInt id = MetaInfo::registry().registerName("slope_of_baseline");
      return id;
    }

    inline UInt startPositionAt5()
    {
      static const UInt id = MetaInfo::registry().registerName("start_position_at_5");
      return id;
    }

    inline UInt startPositionAt10()
    {
      static const UInt id = MetaInfo::registry().registerName("start_position_at_10");
      return id;
    }

    inline UInt startPositionAt50()
    {
      static const UInt id = MetaInfo::registry().registerName("start_position_at_50");
      return id;
    }

    inline UInt tailingFactor()
    {
      static const UInt id = MetaInfo::registry().registerName("tailing_factor");
      return id;
    }

    inline UInt totalMI()
    {
      static const UInt id = MetaInfo::registry().registerName("total_mi");
      return id;
    }

    inline UInt totalWidth()
    {
      static const UInt id = MetaInfo::registry().registerName("total_width");
      return id;
    }

    inline UInt totalXIC()
    {
      static const UInt id = MetaInfo::registry().registerName("total_xic");
      return id;
    }

    inline UInt widthAt5()
    {
      static const UInt id = MetaInfo::registry().registerName("width_at_5");
      return id;
    }

    inline UInt widthAt10()
    {
      static const UInt id = MetaInfo::registry().registerName("width_at_10");
      return id;
    }

    inline UInt widthAt50()
    {
      static const UInt id = MetaInfo::registry().registerName("width_at_50");
      return id;
    }
  }


  /**

    @brief The MRMTransitionGroupPicker finds peaks in chromatograms that belong to the same precursors.

    @htmlinclude OpenMS_MRMTransitionGroupPicker.parameters

    It is called through pickTransitionGroup which will accept an
    MRMTransitionGroup filled with n chromatograms and perform the following steps:
     - Step 1: find features (peaks) in individual chromatograms
     - Step 2: merge these features to consensus features that span multiple chromatograms

    Step 1 is performed by smoothing the individual chromatogram and applying the
    PeakPickerHiRes.

    Step 2 is performed by finding the largest peak overall and use this to
    create a feature, propagating this through all chromatograms.

  */
  class OPENMS_DLLAPI MRMTransitionGroupPicker :
    public DefaultParamHandler
  {

public:

    //@{
    /// Constructor
    MRMTransitionGroupPicker();

    /// Destructor
    ~MRMTransitionGroupPicker() override;
    //@}

    /**
      @brief Pick a group of chromatograms belonging to the same peptide

      Will identify peaks in a set of chromatograms that belong to the same
      peptide. The chromatograms are given in the MRMTransitionGroup container
      which also contains the mapping of the chromatograms to their metadata.
      Only chromatograms from detecting transitions are used for peak picking.
      Identifying transitions will be processed alongside but do not contribute
      to the meta-data, e.g. total_xic or peak_apices_sum.

      The resulting features are added to the MRMTransitionGroup. Each feature contains the following meta-data:

      - PeptideRef
      - leftWidth
      - rightWidth
      - total_xic (fragment trace XIC sum)
      - peak_apices_sum

    */
    template <typename SpectrumT, typename TransitionT>
    void pickTransitionGroup(MRMTransitionGroup<SpectrumT, TransitionT>& transition_group)
    {
      OPENMS_PRECONDITION(transition_group.isInternallyConsistent(), "Consistent state required")
      OPENMS_PRECONDITION(transition_group.chromatogramIdsMatch(), "Chromatogram native IDs need to match keys in transition group")

      // Reuse thread-local pools for large temporary containers to avoid
      // per-call allocations in multi-threaded runs.
      // Use shared per-thread pool declared at file scope to allow nested
      // helpers (computeQuality_, pickFragmentChromatograms, ...) to reuse
      // the same buffers and avoid repeated allocations.
      auto &picker_pool = g_picker_pool_shared();
      picker_pool.reset();
      auto &picked_chroms = picker_pool.picked_chroms;
      auto &smoothed_chroms = picker_pool.smoothed_chroms;
      auto &picked_input_chromatograms = picker_pool.picked_input_chromatograms;
      const Size expected_chrom_count = transition_group.getChromatograms().size() +
        (use_precursors_ ? transition_group.getPrecursorChromatograms().size() : 0);
      picked_chroms.reserve(expected_chrom_count);
      smoothed_chroms.reserve(expected_chrom_count);
      picked_input_chromatograms.reserve(expected_chrom_count);

      // Pick fragment ion chromatograms
      for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
      {
        MSChromatogram& chromatogram = transition_group.getChromatograms()[k];
        String native_id = chromatogram.getNativeID();

        // only pick detecting transitions (skip all others)
        if (transition_group.getTransitions().size() > 0 && 
            transition_group.hasTransition(native_id)  && 
            !transition_group.getTransition(native_id).isDetectingTransition() )
        {
          continue;
        }

        MSChromatogram picked_chrom, smoothed_chrom;
        smoothed_chrom.setNativeID(native_id);
        picker_.pickChromatogram(chromatogram, picked_chrom, smoothed_chrom);
        picked_chrom.sortByIntensity();
        picked_input_chromatograms.push_back(&chromatogram);
        picked_chroms.push_back(std::move(picked_chrom));
        smoothed_chroms.push_back(std::move(smoothed_chrom));
      }

      // Pick precursor chromatograms
      if (use_precursors_)
      {
        for (Size k = 0; k < transition_group.getPrecursorChromatograms().size(); k++)
        {
          SpectrumT picked_chrom, smoothed_chrom;
          SpectrumT& chromatogram = transition_group.getPrecursorChromatograms()[k];

          picker_.pickChromatogram(chromatogram, picked_chrom, smoothed_chrom);
          picked_chrom.sortByIntensity();
          picked_input_chromatograms.push_back(&chromatogram);
          picked_chroms.push_back(std::move(picked_chrom));
          smoothed_chroms.push_back(std::move(smoothed_chrom));
        }
      }

      auto& transition_total_xics = picker_pool.transition_total_xics;
      transition_total_xics.reserve(transition_group.getTransitions().size());
      auto& transition_chromatograms = picker_pool.transition_chromatograms;
      transition_chromatograms.reserve(transition_group.getTransitions().size());
      double detecting_total_xic = 0.0;
      for (const TransitionT& transition : transition_group.getTransitions())
      {
        const SpectrumT& chromatogram = selectChromHelper_(transition_group, transition.getNativeID());
        transition_chromatograms.push_back(&chromatogram);
        const double transition_total_xic = std::accumulate(chromatogram.begin(), chromatogram.end(), 0.0,
                                                            [](double total, const auto& peak)
                                                            {
                                                              return total + peak.getIntensity();
                                                            });
        transition_total_xics.push_back(transition_total_xic);
        if (transition.isDetectingTransition())
        {
          detecting_total_xic += transition_total_xic;
        }
      }

      auto& largest_peak_indices = picker_pool.largest_peak_indices;
      if (boundary_selection_method_ == "largest")
      {
        largest_peak_indices.reserve(picked_chroms.size());
        for (const SpectrumT& chromatogram : picked_chroms)
        {
          largest_peak_indices.push_back(chromatogram.empty() ? -1 : static_cast<int>(chromatogram.size() - 1));
        }
      }

      // Find features (peak groups) in this group of transitions.
      // While there are still peaks left, one will be picked and used to create
      // a feature. Whenever we run out of peaks, we will get -1 back as index
      // and terminate.
      int chr_idx, peak_idx, cnt = 0;
      auto &features = picker_pool.features;
      while (true)
      {
        chr_idx = -1; peak_idx = -1;

        if (boundary_selection_method_ == "largest")
        {
          findLargestPeakCached_(picked_chroms, largest_peak_indices, chr_idx, peak_idx);
        }
        else if (boundary_selection_method_ == "widest")
        {
          findWidestPeakIndices(picked_chroms, chr_idx, peak_idx);
        }

        if (chr_idx == -1 && peak_idx == -1)
        {
          OPENMS_LOG_DEBUG << "**** No more peaks left. Nr. chroms: " << picked_chroms.size() << std::endl;
          break;
        }

        // Compute a feature from the individual chromatograms and add non-zero features
        MRMFeature mrm_feature = createMRMFeature(transition_group, picked_chroms, smoothed_chroms,
                                                  transition_total_xics, detecting_total_xic,
                                                  transition_chromatograms,
                                                  picked_input_chromatograms,
                                                  chr_idx, peak_idx);
        double total_xic = 0;
        double intensity = mrm_feature.getIntensity();
        if (intensity > 0)
        {
          total_xic = mrm_feature.getMetaValue(MRMTransitionGroupPickerMeta::totalXIC());
          features.push_back(std::move(mrm_feature));
          cnt++;
        }

        if (stop_after_feature_ > 0 && cnt > stop_after_feature_) {
          // If you set this, you only expect one feature anyway. No logging necessary why it stopped.
          break;
        }
        if (intensity > 0 && intensity / total_xic < stop_after_intensity_ratio_)
        {
          OPENMS_LOG_DEBUG << "**** Minimum intensity ratio reached. Nr. chroms: " << picked_chroms.size() << std::endl;
          break;
        }
      }

      // Check for completely overlapping features
      for (Size i = 0; i < features.size(); i++)
      {
        MRMFeature& mrm_feature = features[i];
        bool skip = false;
        for (Size j = 0; j < i; j++)
        {
          if ((double)mrm_feature.getMetaValue(MRMTransitionGroupPickerMeta::leftWidth()) >=  (double)features[j].getMetaValue(MRMTransitionGroupPickerMeta::leftWidth()) &&
              (double)mrm_feature.getMetaValue(MRMTransitionGroupPickerMeta::rightWidth()) <= (double)features[j].getMetaValue(MRMTransitionGroupPickerMeta::rightWidth()))
          { skip = true; }
        }
        if (mrm_feature.getIntensity() > 0 && !skip)
        {
          transition_group.addFeature(mrm_feature);
        }
      }

    }

    /**
      @brief Create a feature from picked chromatograms and a specified peak.

      Compatibility overload for callers that do not precompute transition XICs
      and chromatogram pointers. The optimized scoring path calls the overload
      below directly to avoid rebuilding these inputs for every peak group.

      @param[in] transition_group The transition group that owns the raw chromatograms
      @param[in,out] picked_chroms Picked chromatograms; selected peaks are marked as consumed
      @param[in] smoothed_chroms Smoothed chromatograms used for integration
      @param[in] chr_idx Index of the picked chromatogram containing the seed peak
      @param[in] peak_idx Index of the seed peak within @p picked_chroms[chr_idx]

      @return The created MRM feature, or an empty feature with zero intensity if
              filtering rejects it
    */
    template <typename SpectrumT, typename TransitionT>
    MRMFeature createMRMFeature(const MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                                std::vector<SpectrumT>& picked_chroms,
                                const std::vector<SpectrumT>& smoothed_chroms,
                                const int chr_idx,
                                const int peak_idx)
    {
      std::vector<double> transition_total_xics;
      transition_total_xics.reserve(transition_group.getTransitions().size());
      std::vector<const SpectrumT*> transition_chromatograms;
      transition_chromatograms.reserve(transition_group.getTransitions().size());
      double detecting_total_xic = 0.0;
      for (const TransitionT& transition : transition_group.getTransitions())
      {
        const SpectrumT& chromatogram = selectChromHelper_(transition_group, transition.getNativeID());
        transition_chromatograms.push_back(&chromatogram);
        const double transition_total_xic = std::accumulate(chromatogram.begin(), chromatogram.end(), 0.0,
                                                            [](double total, const auto& peak)
                                                            {
                                                              return total + peak.getIntensity();
                                                            });
        transition_total_xics.push_back(transition_total_xic);
        if (transition.isDetectingTransition())
        {
          detecting_total_xic += transition_total_xic;
        }
      }

      std::vector<const SpectrumT*> picked_input_chromatograms;
      picked_input_chromatograms.reserve(picked_chroms.size());
      for (const SpectrumT& picked_chrom : picked_chroms)
      {
        picked_input_chromatograms.push_back(&selectChromHelper_(transition_group, picked_chrom.getNativeID()));
      }

      return createMRMFeature(transition_group, picked_chroms, smoothed_chroms,
                              transition_total_xics, detecting_total_xic,
                              transition_chromatograms, picked_input_chromatograms,
                              chr_idx, peak_idx);
    }

    /**
      @brief Create a feature from picked chromatograms and precomputed transition data.

      This overload is used by the optimized peak picking path. The caller passes
      total XIC values and raw chromatogram pointers so repeated peak groups can
      reuse the same data without repeated lookups and allocations.

      @param[in] transition_group The transition group that owns the raw chromatograms
      @param[in,out] picked_chroms Picked chromatograms; selected peaks are marked as consumed
      @param[in] smoothed_chroms Smoothed chromatograms used for integration
      @param[in] transition_total_xics Precomputed total XIC per transition
      @param[in] detecting_total_xic Sum of total XIC values for detecting transitions
      @param[in] transition_chromatograms Raw chromatograms corresponding to transitions
      @param[in] picked_input_chromatograms Raw chromatograms corresponding to @p picked_chroms
      @param[in] chr_idx Index of the picked chromatogram containing the seed peak
      @param[in] peak_idx Index of the seed peak within @p picked_chroms[chr_idx]

      @return The created MRM feature, or an empty feature with zero intensity if
              filtering rejects it
    */
    template <typename SpectrumT, typename TransitionT>
    MRMFeature createMRMFeature(const MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                                std::vector<SpectrumT>& picked_chroms,
                                const std::vector<SpectrumT>& smoothed_chroms,
                                const std::vector<double>& transition_total_xics,
                                double detecting_total_xic,
                                const std::vector<const SpectrumT*>& transition_chromatograms,
                                const std::vector<const SpectrumT*>& picked_input_chromatograms,
                                const int chr_idx,
                                const int peak_idx)
    {
      OPENMS_PRECONDITION(transition_group.isInternallyConsistent(), "Consistent state required")
      OPENMS_PRECONDITION(transition_group.chromatogramIdsMatch(), "Chromatogram native IDs need to match keys in transition group")

      MRMFeature mrmFeature;
      mrmFeature.setIntensity(0.0);
      double best_left = picked_chroms[chr_idx].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][peak_idx];
      double best_right = picked_chroms[chr_idx].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][peak_idx];
      double peak_apex = picked_chroms[chr_idx][peak_idx].getRT();
      OPENMS_LOG_DEBUG << "**** Creating MRMFeature for peak " << peak_idx << " in chrom. " << chr_idx << " with " <<
        picked_chroms[chr_idx][peak_idx] << " and borders " << best_left << " " <<
        best_right << " (" << best_right - best_left << ")" << std::endl;

      if (use_consensus_ && recalculate_peaks_)
      {
        // This may change best_left / best_right
        recalculatePeakBorders_(picked_chroms, best_left, best_right, recalculate_peaks_max_z_);
        if (peak_apex < best_left || peak_apex > best_right)
        {
          // apex fell out of range, lets correct it
          peak_apex = (best_left + best_right) / 2.0;
        }
      }

      // Use a small thread-local pool for these per-feature vectors to avoid
      // frequent allocations when called repeatedly.
      struct CreateFeaturePool
      {
        std::vector<double> left_edges;
        std::vector<double> right_edges;
        void reset() { left_edges.clear(); right_edges.clear(); }
      };
      thread_local CreateFeaturePool create_feat_pool;
      create_feat_pool.reset();
      auto &left_edges = create_feat_pool.left_edges;
      auto &right_edges = create_feat_pool.right_edges;
      double min_left = best_left;
      double max_right = best_right;
      if (use_consensus_)
      {
        // Remove other, overlapping, picked peaks (in this and other
        // chromatograms) and then ensure that at least one peak is set to zero
        // (the currently best peak).
        remove_overlapping_features(picked_chroms, best_left, best_right);
      }
      else
      {
        pickApex(picked_chroms, best_left, best_right, peak_apex,
                 min_left, max_right, left_edges, right_edges);

      } // end !use_consensus_
      picked_chroms[chr_idx][peak_idx].setIntensity(0.0); // ensure that we set at least one peak to zero

      // Check for minimal peak width -> return empty feature (Intensity zero)
      if (use_consensus_)
      {
        if (min_peak_width_ > 0.0 && std::fabs(best_right - best_left) < min_peak_width_) 
        {
          return mrmFeature;
        }

        if (compute_peak_quality_)
        {
          String outlier = "none";
          double qual = computeQuality_(transition_group, picked_chroms, picked_input_chromatograms,
                                        chr_idx, best_left, best_right, outlier);
          if (qual < min_qual_) 
          {
            return mrmFeature;
          }
          mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::potentialOutlier(), outlier);
          mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::initialPeakQuality(), qual);
          mrmFeature.setOverallQuality(qual);
        }
      }

      // Prepare linear resampling of all the chromatograms, here creating the
      // empty master_peak_container with the same RT (m/z) values as the
      // reference chromatogram. We use the overall minimal left boundary and
      // maximal right boundary to prepare the container.
      SpectrumT master_peak_container;
      const SpectrumT& ref_chromatogram = *picked_input_chromatograms[chr_idx];
      prepareMasterContainer_(ref_chromatogram, master_peak_container, min_left, max_right);

      // Iterate over initial transitions / chromatograms (note that we may
      // have a different number of picked chromatograms than total transitions
      // as not all are detecting transitions).
      double total_intensity = 0; double total_peak_apices = 0; double total_xic = 0; double total_mi = 0;
      pickFragmentChromatograms(transition_group, picked_chroms, mrmFeature, smoothed_chroms,
                                best_left, best_right, use_consensus_,
                                total_intensity, total_xic, total_mi, total_peak_apices,
                                master_peak_container, left_edges, right_edges,
                                transition_total_xics, detecting_total_xic,
                                transition_chromatograms,
                                chr_idx, peak_idx);

      // Also pick the precursor chromatogram(s); note total_xic is not
      // extracted here, only for fragment traces
      pickPrecursorChromatograms(transition_group,
                                picked_chroms, mrmFeature, smoothed_chroms,
                                best_left, best_right, use_consensus_,
                                total_intensity, master_peak_container, left_edges, right_edges,
                                chr_idx, peak_idx);

      mrmFeature.setRT(peak_apex);
      mrmFeature.setIntensity(total_intensity);
      mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::peptideRef(), transition_group.getTransitionGroupID());
      mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::leftWidth(), best_left);
      mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::rightWidth(), best_right);
      mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::totalXIC(), total_xic);
      if (compute_total_mi_)
      {
        mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::totalMI(), total_mi);
      }
      mrmFeature.setMetaValue(MRMTransitionGroupPickerMeta::peakApicesSum(), total_peak_apices);

      mrmFeature.ensureUniqueId();
      return mrmFeature;
    }

    /** 
     
      @brief Apex-based peak picking

      Pick the peak with the closest apex to the consensus apex for each
      chromatogram.  Use the closest peak for the current peak. 
      
      Note that we will only set the closest peak per chromatogram to zero, so
      if there are two peaks for some transitions, we will have to get to them
      later.  If there is no peak, then we transfer transition boundaries from
      "master" peak.
    */
    template <typename SpectrumT>
    void pickApex(std::vector<SpectrumT>& picked_chroms,
                  const double best_left, const double best_right, const double peak_apex,
                  double &min_left, double &max_right, 
                  std::vector< double > & left_edges, std::vector< double > & right_edges)
    {
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        double peak_apex_dist_min = std::numeric_limits<double>::max();
        int min_dist = -1;
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          PeakIntegrator::PeakArea pa_tmp = pi_.integratePeak(  // get the peak apex
              picked_chroms[k],
              picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][i], 
              picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][i]); 
          if (pa_tmp.apex_pos > 1e-11 && std::fabs(pa_tmp.apex_pos - peak_apex) < peak_apex_dist_min)
          { // update best candidate
            peak_apex_dist_min = std::fabs(pa_tmp.apex_pos - peak_apex);
            min_dist = (int)i;
          }
        }

        // Select master peak boundaries, or in the case we found at least one peak, the local peak boundaries 
        double l = best_left;
        double r = best_right;
        if (min_dist >= 0)
        {
          l = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][min_dist];
          r = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][min_dist];
          picked_chroms[k][min_dist].setIntensity(0.0); // only remove one peak per transition
        }

        left_edges.push_back(l);
        right_edges.push_back(r);
        // ensure we remember the overall maxima / minima
        if (l < min_left) {min_left = l;}
        if (r > max_right) {max_right = r;}
      }
    }

    template <typename SpectrumT, typename TransitionT>
    void pickFragmentChromatograms(const MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                                    const std::vector<SpectrumT>& picked_chroms,
                                    MRMFeature& mrmFeature,
                                    const std::vector<SpectrumT>& smoothed_chroms,
                                    const double best_left, const double best_right,
                                    const bool use_consensus_,
                                    double & total_intensity,
                                    double & total_xic,
                                    double & total_mi,
                                    double & total_peak_apices,
                                    const SpectrumT & master_peak_container,
                                    const std::vector< double > & left_edges,
                                    const std::vector< double > & right_edges,
                                    const std::vector<double>& transition_total_xics,
                                    double detecting_total_xic,
                                    const std::vector<const SpectrumT*>& transition_chromatograms,
                                    const int chr_idx,
                                    const int peak_idx)
    {
      total_xic += detecting_total_xic;
      for (Size k = 0; k < transition_group.getTransitions().size(); k++)
      {

        double local_left = best_left;
        double local_right = best_right;
        if (!use_consensus_)
        {
          // We cannot have any non-detecting transitions (otherwise we have
          // too few left / right edges) as we skipped those when doing peak
          // picking and smoothing.
          if (!transition_group.getTransitions()[k].isDetectingTransition())
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "When using non-consensus peak picker, all transitions need to be detecting transitions.");
          }
          local_left = left_edges[k];
          local_right = right_edges[k];
        }

        const SpectrumT& chromatogram = *transition_chromatograms[k];

        // Compute total intensity on transition-level
        const double transition_total_xic = transition_total_xics[k];

        // Compute total mutual information on transition-level.
        double transition_total_mi = 0;
        if (compute_total_mi_)
        {
          // reuse thread-local temporaries to avoid repeated allocations
          thread_local std::vector<unsigned int> chrom_vect_id_ranked;
          thread_local std::vector<unsigned int> chrom_vect_det_ranked;
          thread_local std::vector<double> chrom_vect_id;
          thread_local std::vector<double> chrom_vect_det;
          chrom_vect_id.clear(); chrom_vect_det.clear(); chrom_vect_id_ranked.clear(); chrom_vect_det_ranked.clear();
          for (typename SpectrumT::const_iterator it = chromatogram.begin(); it != chromatogram.end(); it++)
          {
            chrom_vect_id.push_back(it->getIntensity());
          }
          unsigned int max_rank_det = OpenSwath::Scoring::computeAndAppendRank(chrom_vect_id, chrom_vect_det_ranked);
          // compute baseline mutual information
          int transition_total_mi_norm = 0;
          for (Size m = 0; m < transition_group.getTransitions().size(); m++)
          {
            if (transition_group.getTransitions()[m].isDetectingTransition())
            {
              const SpectrumT& chromatogram_det = *transition_chromatograms[m];
              chrom_vect_det.clear();
              for (typename SpectrumT::const_iterator it = chromatogram_det.begin(); it != chromatogram_det.end(); it++)
              {
                chrom_vect_det.push_back(it->getIntensity());
              }
              unsigned int max_rank_id = OpenSwath::Scoring::computeAndAppendRank(chrom_vect_det, chrom_vect_id_ranked);
              transition_total_mi += OpenSwath::Scoring::rankedMutualInformation(chrom_vect_id_ranked, chrom_vect_det_ranked, max_rank_id, max_rank_det);
              transition_total_mi_norm++;
            }
          }
          if (transition_total_mi_norm > 0) { transition_total_mi /= transition_total_mi_norm; }

          if (transition_group.getTransitions()[k].isDetectingTransition())
          {
            // sum up all transition-level total MI and divide by the number of detection transitions to have peak group level total MI
            total_mi += transition_total_mi / transition_total_mi_norm;
          }
        }

        SpectrumT used_chromatogram;
        // resample the current chromatogram
        if (peak_integration_ == "original")
        {
          used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, local_left, local_right);
        }
        else if (peak_integration_ == "smoothed")
        {
          if (smoothed_chroms.size() <= k) 
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "Tried to calculate peak area and height without any smoothed chromatograms");
          }
          used_chromatogram = resampleChromatogram_(smoothed_chroms[k], master_peak_container, local_left, local_right);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            String("Peak integration chromatogram ") + peak_integration_ + " is not a valid method for MRMTransitionGroupPicker");
        } 

        Feature f;
        double quality = 0;
        f.setQuality(0, quality);
        f.setOverallQuality(quality);

        PeakIntegrator::PeakArea pa = pi_.integratePeak(used_chromatogram, local_left, local_right);
        double peak_integral = pa.area;
        double peak_apex_int = pa.height;
        f.setMetaValue(MRMTransitionGroupPickerMeta::peakApexPosition(), pa.apex_pos);

        // In non-consensus mode, override hull_points with the full resampled
        // chromatogram so all transitions produce equal-sized intensity vectors
        // for cross-correlation scoring. Area/height/apex are already correct
        // from the per-transition integration above.
        if (!use_consensus_)
        {
          pa.hull_points.clear();
          for (const auto& p : used_chromatogram)
          {
            pa.hull_points.push_back(DPosition<2>(p.getPos(), p.getIntensity()));
          }
        }

        if (background_subtraction_ != "none")
        {
          double background{0};
          double avg_noise_level{0};
          if (background_subtraction_ == "original")
          {
            const double intensity_left = chromatogram.PosBegin(local_left)->getIntensity();
            const double intensity_right = (chromatogram.PosEnd(local_right) - 1)->getIntensity();
            const UInt n_points = std::distance(chromatogram.PosBegin(local_left), chromatogram.PosEnd(local_right));
            avg_noise_level = (intensity_right + intensity_left) / 2;
            background = avg_noise_level * n_points;
          }
          else if (background_subtraction_ == "exact")
          {
            PeakIntegrator::PeakBackground pb = pi_.estimateBackground(used_chromatogram, local_left, local_right, pa.apex_pos);
            background = pb.area;
            avg_noise_level = pb.height;
          }
          peak_integral -= background;
          peak_apex_int -= avg_noise_level;
          if (peak_integral < 0) {peak_integral = 0;}
          if (peak_apex_int < 0) {peak_apex_int = 0;}

          f.setMetaValue(MRMTransitionGroupPickerMeta::areaBackgroundLevel(), background);
          f.setMetaValue(MRMTransitionGroupPickerMeta::noiseBackgroundLevel(), avg_noise_level);
        } // end background

        f.setRT(picked_chroms[chr_idx][peak_idx].getPos());
        f.setIntensity(peak_integral);
        ConvexHull2D hull;
        hull.setHullPoints(pa.hull_points);
        f.getConvexHulls().push_back(hull);

        f.setMZ(chromatogram.getProduct().getMZ());
        mrmFeature.setMZ(chromatogram.getPrecursor().getMZ());

        if (chromatogram.metaValueExists("product_mz")) // legacy code (ensures that old tests still work)
        {
          f.setMetaValue("MZ", chromatogram.getMetaValue("product_mz"));
          f.setMZ(chromatogram.getMetaValue("product_mz"));
        }

        f.setMetaValue(MRMTransitionGroupPickerMeta::nativeID(), chromatogram.getNativeID());
        f.setMetaValue(MRMTransitionGroupPickerMeta::peakApexInt(), peak_apex_int);
        f.setMetaValue(MRMTransitionGroupPickerMeta::totalXIC(), transition_total_xic);
        if (compute_total_mi_)
        {
          f.setMetaValue(MRMTransitionGroupPickerMeta::totalMI(), transition_total_mi);
        }

        if (transition_group.getTransitions()[k].isQuantifyingTransition())
        {
          total_intensity += peak_integral;
          total_peak_apices += peak_apex_int;
        }

        // for backwards compatibility with TOPP tests
        // Calculate peak shape metrics that will be used for later QC
        PeakIntegrator::PeakShapeMetrics psm = pi_.calculatePeakShapeMetrics(used_chromatogram, local_left, local_right, peak_apex_int, pa.apex_pos);
        f.setMetaValue(MRMTransitionGroupPickerMeta::widthAt50(), psm.width_at_50);
        if (compute_peak_shape_metrics_)
        {
          f.setMetaValue(MRMTransitionGroupPickerMeta::widthAt5(), psm.width_at_5);
          f.setMetaValue(MRMTransitionGroupPickerMeta::widthAt10(), psm.width_at_10);
          f.setMetaValue(MRMTransitionGroupPickerMeta::startPositionAt5(), psm.start_position_at_5);
          f.setMetaValue(MRMTransitionGroupPickerMeta::startPositionAt10(), psm.start_position_at_10);
          f.setMetaValue(MRMTransitionGroupPickerMeta::startPositionAt50(), psm.start_position_at_50);
          f.setMetaValue(MRMTransitionGroupPickerMeta::endPositionAt5(), psm.end_position_at_5);
          f.setMetaValue(MRMTransitionGroupPickerMeta::endPositionAt10(), psm.end_position_at_10);
          f.setMetaValue(MRMTransitionGroupPickerMeta::endPositionAt50(), psm.end_position_at_50);
          f.setMetaValue(MRMTransitionGroupPickerMeta::totalWidth(), psm.total_width);
          f.setMetaValue(MRMTransitionGroupPickerMeta::tailingFactor(), psm.tailing_factor);
          f.setMetaValue(MRMTransitionGroupPickerMeta::asymmetryFactor(), psm.asymmetry_factor);
          f.setMetaValue(MRMTransitionGroupPickerMeta::slopeOfBaseline(), psm.slope_of_baseline);
          f.setMetaValue(MRMTransitionGroupPickerMeta::baselineDelta2Height(), psm.baseline_delta_2_height);
          f.setMetaValue(MRMTransitionGroupPickerMeta::pointsAcrossBaseline(), psm.points_across_baseline);
          f.setMetaValue(MRMTransitionGroupPickerMeta::pointsAcrossHalfHeight(), psm.points_across_half_height);
        }

        mrmFeature.addFeature(f, chromatogram.getNativeID()); //map index and feature
      }
    }

    template <typename SpectrumT, typename TransitionT>
    void pickPrecursorChromatograms(const MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                                    const std::vector<SpectrumT>& picked_chroms,
                                    MRMFeature& mrmFeature,
                                    const std::vector<SpectrumT>& smoothed_chroms,
                                    const double best_left, const double best_right,
                                    const bool use_consensus_,
                                    double & total_intensity,
                                    const SpectrumT & master_peak_container,
                                    const std::vector< double > & left_edges,
                                    const std::vector< double > & right_edges,
                                    const int chr_idx,
                                    const int peak_idx)
    {
      for (Size k = 0; k < transition_group.getPrecursorChromatograms().size(); k++)
      {
        const SpectrumT& chromatogram = transition_group.getPrecursorChromatograms()[k];

        // Identify precursor index
        // note: this is only valid if all transitions are detecting transitions
        Size prec_idx = transition_group.getChromatograms().size() + k;

        double local_left = best_left;
        double local_right = best_right;
        if (!use_consensus_ && right_edges.size() > prec_idx && left_edges.size() > prec_idx)
        {
          local_left = left_edges[prec_idx];
          local_right = right_edges[prec_idx];
        }

        SpectrumT used_chromatogram;
        // resample the current chromatogram
        if (peak_integration_ == "original")
        {
          used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, local_left, local_right);
          // const SpectrumT& used_chromatogram = chromatogram; // instead of resampling
        }
        else if (peak_integration_ == "smoothed" && smoothed_chroms.size() <= prec_idx)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Tried to calculate peak area and height without any smoothed chromatograms for precursors");
        }
        else if (peak_integration_ == "smoothed")
        {
          used_chromatogram = resampleChromatogram_(smoothed_chroms[prec_idx], master_peak_container, local_left, local_right);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            String("Peak integration chromatogram ") + peak_integration_ + " is not a valid method for MRMTransitionGroupPicker");
        }

        Feature f;
        double quality = 0;
        f.setQuality(0, quality);
        f.setOverallQuality(quality);

        PeakIntegrator::PeakArea pa = pi_.integratePeak(used_chromatogram, local_left, local_right);
        double peak_integral = pa.area;
        double peak_apex_int = pa.height;

        // In non-consensus mode, override hull_points with the full resampled
        // chromatogram so all precursors produce equal-sized intensity vectors
        // for cross-correlation scoring (see pickFragmentChromatograms).
        if (!use_consensus_)
        {
          pa.hull_points.clear();
          for (const auto& p : used_chromatogram)
          {
            pa.hull_points.push_back(DPosition<2>(p.getPos(), p.getIntensity()));
          }
        }

        if (background_subtraction_ != "none")
        {
          double background{0};
          double avg_noise_level{0};
          if ((peak_integration_ == "smoothed") && smoothed_chroms.size() <= prec_idx)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Tried to calculate background estimation without any smoothed chromatograms");
          }
          else if (background_subtraction_ == "original")
          {
            const double intensity_left = chromatogram.PosBegin(local_left)->getIntensity();
            const double intensity_right = (chromatogram.PosEnd(local_right) - 1)->getIntensity();
            const UInt n_points = std::distance(chromatogram.PosBegin(local_left), chromatogram.PosEnd(local_right));
            avg_noise_level = (intensity_right + intensity_left) / 2;
            background = avg_noise_level * n_points;
          }
          else if (background_subtraction_ == "exact")
          {
            PeakIntegrator::PeakBackground pb = pi_.estimateBackground(used_chromatogram, local_left, local_right, pa.apex_pos);
            background = pb.area;
            avg_noise_level = pb.height;
          }
          peak_integral -= background;
          peak_apex_int -= avg_noise_level;
          if (peak_integral < 0) {peak_integral = 0;}
          if (peak_apex_int < 0) {peak_apex_int = 0;}

          f.setMetaValue(MRMTransitionGroupPickerMeta::areaBackgroundLevel(), background);
          f.setMetaValue(MRMTransitionGroupPickerMeta::noiseBackgroundLevel(), avg_noise_level);
        }

        f.setMZ(chromatogram.getPrecursor().getMZ());
        if (k == 0) {mrmFeature.setMZ(chromatogram.getPrecursor().getMZ());} // only use m/z if first (monoisotopic) isotope

        if (chromatogram.metaValueExists("precursor_mz")) // legacy code (ensures that old tests still work)
        {
          f.setMZ(chromatogram.getMetaValue("precursor_mz"));
          if (k == 0) {mrmFeature.setMZ(chromatogram.getMetaValue("precursor_mz"));} // only use m/z if first (monoisotopic) isotope
        }

        f.setRT(picked_chroms[chr_idx][peak_idx].getPos());
        f.setIntensity(peak_integral);
        ConvexHull2D hull;
        hull.setHullPoints(pa.hull_points);
        f.getConvexHulls().push_back(hull);
        f.setMetaValue(MRMTransitionGroupPickerMeta::nativeID(), chromatogram.getNativeID());
        f.setMetaValue(MRMTransitionGroupPickerMeta::peakApexInt(), peak_apex_int);

        if (use_precursors_ && transition_group.getTransitions().empty())
        {
          total_intensity += peak_integral;
        }

        mrmFeature.addPrecursorFeature(f, chromatogram.getNativeID());
      }
    }

    // maybe private, but we have tests

    /**
      @brief Remove overlapping features.

      Remove features that are within the current seed (between best_left and
      best_right) or overlap with it. An overlapping feature is defined as a
      feature that has either of its borders within the border of the current
      peak

      Directly adjacent features are allowed, e.g. they can share one
      border.
    */
    template <typename SpectrumT>
    void remove_overlapping_features(std::vector<SpectrumT>& picked_chroms, double best_left, double best_right)
    {
      // delete all seeds that lie within the current seed
      Size count_inside = 0;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getPos() >= best_left && picked_chroms[k][i].getPos() <= best_right)
          {
            picked_chroms[k][i].setIntensity(0.0);
            count_inside++;
          }
        }
      }

      Size count_overlap = 0;
      // delete all seeds that overlap within the current seed
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getIntensity() <= 1e-11) {continue; }

          double left = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][i];
          double right = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][i];
          if ((left > best_left && left < best_right)
             || (right > best_left && right < best_right))
          {
            picked_chroms[k][i].setIntensity(0.0);
            count_overlap++;
          }
        }
      }
      OPENMS_LOG_DEBUG << " ** Removed " << count_inside << " peaks enclosed in and " <<
       count_overlap << " peaks overlapping with added feature." << std::endl;
    }

    /// Find largest peak in a vector of chromatograms
    void findLargestPeak(const std::vector<MSChromatogram >& picked_chroms, int& chr_idx, int& peak_idx);

    template <typename SpectrumT>
    void findLargestPeakCached_(const std::vector<SpectrumT>& picked_chroms,
                                std::vector<int>& largest_peak_indices,
                                int& chr_idx,
                                int& peak_idx) const
    {
      OPENMS_PRECONDITION(largest_peak_indices.size() == picked_chroms.size(),
                          "Cached largest peak indices must match picked chromatogram count.");
      double largest = 0.0;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        int candidate_idx = largest_peak_indices[k];
        while (candidate_idx >= 0 && picked_chroms[k][candidate_idx].getIntensity() <= 0.0)
        {
          --candidate_idx;
        }
        largest_peak_indices[k] = candidate_idx;
        if (candidate_idx < 0)
        {
          continue;
        }
        const double intensity = picked_chroms[k][candidate_idx].getIntensity();
        if (intensity > largest)
        {
          largest = intensity;
          chr_idx = static_cast<int>(k);
          peak_idx = candidate_idx;
        }
      }
    }

    /**
      @brief Given a vector of chromatograms, find the indices of the chromatogram
      containing the widest peak and of the position of highest intensity.

      @param[in] picked_chroms The vector of chromatograms
      @param[out] chrom_idx The index of the chromatogram containing the widest peak
      @param[out] point_idx The index of the point with highest intensity
    */
    void findWidestPeakIndices(const std::vector<MSChromatogram>& picked_chroms, Int& chrom_idx, Int& point_idx) const;

protected:

    /// Synchronize members with param class
    void updateMembers_() override;

    /// Assignment operator is protected for algorithm
    MRMTransitionGroupPicker& operator=(const MRMTransitionGroupPicker& rhs);

    /**
      @brief Select matching precursor or fragment ion chromatogram
    */
    template <typename SpectrumT, typename TransitionT>
    const SpectrumT& selectChromHelper_(const MRMTransitionGroup<SpectrumT, TransitionT>& transition_group, const String& native_id)
    {
      if (transition_group.hasChromatogram(native_id))
      {
        return transition_group.getChromatogram(native_id);
      }
      else if (transition_group.hasPrecursorChromatogram(native_id))
      {
        return transition_group.getPrecursorChromatogram(native_id);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Did not find chromatogram for id '" + native_id + "'.");
      }
    }

    /**
      @brief Compute transition group quality (higher score is better)

      This is only based on the co-elution of the chromatograms and internal
      consistency without any library information.

      For the final score (larger is better), consider these scores:
      - missing_peaks (the more peaks are missing, the worse)
      - multiple_peaks
      - mean of the shapes (1 is very good, 0 is bad)
      - mean of the coelutions (0 is good, 1 is ok, above 1 is pretty bad)

      These scores are similar to the ones computed by MRMFeatureFinderScoring
      and a simple sum of these scores is returned.

    */
    template <typename SpectrumT, typename TransitionT>
    double computeQuality_(const MRMTransitionGroup<SpectrumT, TransitionT>&,
                           const std::vector<SpectrumT>& picked_chroms,
                           const std::vector<const SpectrumT*>& picked_input_chromatograms,
                           const int chr_idx,
                           const double best_left,
                           const double best_right,
                           String& outlier)
    {
      // Resample all chromatograms around the current estimated peak and
      // collect the raw intensities. For resampling, use a bit more on either
      // side to correctly identify shoulders etc.
      double resample_boundary = resample_boundary_; // sample 15 seconds more on each side
      SpectrumT master_peak_container;
      const SpectrumT& ref_chromatogram = *picked_input_chromatograms[chr_idx];
      prepareMasterContainer_(ref_chromatogram, master_peak_container, best_left - resample_boundary, best_right + resample_boundary);
      // Reuse the shared per-thread pool's containers for quality computation
      auto &pool = g_picker_pool_shared();
      auto &all_ints = pool.all_ints;
      all_ints.clear();
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        const SpectrumT& chromatogram = *picked_input_chromatograms[k];
        const SpectrumT used_chromatogram = resampleChromatogram_(chromatogram, 
            master_peak_container, best_left - resample_boundary, best_right + resample_boundary);

        std::vector<double> int_here;
        for (const auto& peak : used_chromatogram) int_here.push_back(peak.getIntensity());
        // Remove chromatograms without a single peak
        double tic = std::accumulate(int_here.begin(), int_here.end(), 0.0);
        if (tic > 1e-11) all_ints.push_back(int_here);
      }

      // Compute the cross-correlation for the collected intensities
      auto &mean_shapes = pool.mean_shapes;
      auto &mean_coel = pool.mean_coel;
      mean_shapes.clear();
      mean_coel.clear();
      for (Size k = 0; k < all_ints.size(); k++)
      {
        std::vector<double> shapes;
        std::vector<double> coel;
        for (Size i = 0; i < all_ints.size(); i++)
        {
          if (i == k) {continue;}
          OpenSwath::Scoring::XCorrArrayType res = OpenSwath::Scoring::normalizedCrossCorrelation(
              all_ints[k], all_ints[i], boost::numeric_cast<int>(all_ints[i].size()), 1);

          // the first value is the x-axis (retention time) and should be an int -> it show the lag between the two
          double res_coelution = std::abs(OpenSwath::Scoring::xcorrArrayGetMaxPeak(res)->first);
          double res_shape = std::abs(OpenSwath::Scoring::xcorrArrayGetMaxPeak(res)->second);

          shapes.push_back(res_shape);
          coel.push_back(res_coelution);
        }

        // We have computed the cross-correlation of chromatogram k against
        // all others. Use the mean of these computations as the value for k.
        OpenSwath::mean_and_stddev msc;
        msc = std::for_each(shapes.begin(), shapes.end(), msc);
        double shapes_mean = msc.mean();
        msc = std::for_each(coel.begin(), coel.end(), msc);
        double coel_mean = msc.mean();

        // mean shape scores below 0.5-0.6 should be a real sign of trouble ... !
        // mean coel scores above 3.5 should be a real sign of trouble ... !
        mean_shapes.push_back(shapes_mean);
        mean_coel.push_back(coel_mean);
      }

      // find the chromatogram with the minimal shape score and the maximal
      // coelution score -> if it is the same chromatogram, the chance is
      // pretty good that it is different from the others...
      int min_index_shape = std::distance(mean_shapes.begin(), std::min_element(mean_shapes.begin(), mean_shapes.end()));
      int max_index_coel = std::distance(mean_coel.begin(), std::max_element(mean_coel.begin(), mean_coel.end()));

      // Look at the picked peaks that are within the current left/right borders
      int missing_peaks = 0;
      int multiple_peaks = 0;

      // collect all seeds that lie within the current seed
      auto &left_borders = pool.left_borders;
      auto &right_borders = pool.right_borders;
      left_borders.clear();
      right_borders.clear();
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        double l_tmp;
        double r_tmp;
        double max_int = -1;

        int pfound = 0;
        l_tmp = -1;
        r_tmp = -1;
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getPos() >= best_left && picked_chroms[k][i].getPos() <= best_right)
          {
            pfound++;
            if (picked_chroms[k][i].getIntensity() > max_int)
            {
              max_int = picked_chroms[k][i].getIntensity() > max_int;
              l_tmp = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][i];
              r_tmp = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][i];
            }
          }
        }

        if (l_tmp > 1e-11) left_borders.push_back(l_tmp);
        if (r_tmp > 1e-11) right_borders.push_back(r_tmp);

        if (pfound == 0) missing_peaks++;
        if (pfound > 1) multiple_peaks++;
      }

      // Check how many chromatograms had exactly one peak picked between our
      // current left/right borders -> this would be a sign of consistency.
      OPENMS_LOG_DEBUG << " Overall found missing : " << missing_peaks << " and multiple : " << multiple_peaks << std::endl;

      /// left_borders / right_borders might not have the same length since we might have peaks missing!!

      // Is there one transitions that is very different from the rest (e.g.
      // the same element has a bad shape and a bad coelution score) -> potential outlier
      if (min_index_shape == max_index_coel)
      {
        OPENMS_LOG_DEBUG << " Element " << min_index_shape << " is a candidate for removal ... " << std::endl;
        outlier = String(picked_chroms[min_index_shape].getNativeID());
      }
      else
      {
        outlier = "none";
      }

      // For the final score (larger is better), consider these scores:
      // - missing_peaks (the more peaks are missing, the worse)
      // - multiple_peaks
      // - mean of the shapes (1 is very good, 0 is bad)
      // - mean of the co-elution scores (0 is good, 1 is ok, above 1 is pretty bad)
      double shape_score = std::accumulate(mean_shapes.begin(), mean_shapes.end(), 0.0) / mean_shapes.size();
      double coel_score = std::accumulate(mean_coel.begin(), mean_coel.end(), 0.0) / mean_coel.size();
      coel_score = (coel_score - 1.0) / 2.0;

      double score = shape_score - coel_score - 1.0 * missing_peaks / picked_chroms.size();

      OPENMS_LOG_DEBUG << " Computed score  " << score << " (from " <<  shape_score << 
        " - " << coel_score << " - " << 1.0 * missing_peaks / picked_chroms.size() << ")" << std::endl;

      return score;
    }

    /**
      @brief Recalculate the borders of the peak

      By collecting all left and right borders of contained peaks, a consensus
      peak is computed.  By looking at the means and standard deviations of all
      the peak borders it is estimated whether the proposed peak border
      deviates too much from the consensus one. If the deviation is too high
      (in this case), then we fall back to the "consensus" (a median here).
    */
    template <typename SpectrumT>
    void recalculatePeakBorders_(const std::vector<SpectrumT>& picked_chroms, double& best_left, double& best_right, double max_z)
    {
      // 1. Collect all seeds that lie within the current seed 
      // - Per chromatogram only the most intense one counts, otherwise very
      // - low intense peaks can contribute disproportionally to the voting
      // - procedure.
      // TODO Especially with DDA MS1 "transitions" from FFID, you often have exactly the same
      // borders and the logs get confusing and you get -Nan CVs. You also might save time if you use a set here.
      std::vector<double> left_borders;
      std::vector<double> right_borders;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        double max_int = -1;
        double left = -1;
        double right = -1;
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getPos() >= best_left && picked_chroms[k][i].getPos() <= best_right)
          {
            if (picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][i] > max_int)
            {
              max_int = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE][i];
              left = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][i];
              right = picked_chroms[k].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][i];
            }
          }
        }
        if (max_int > -1 )
        {
          left_borders.push_back(left);
          right_borders.push_back(right);
          OPENMS_LOG_DEBUG << " * peak " << k << " left boundary " << left_borders.back()   <<  " with inty " << max_int << std::endl;
          OPENMS_LOG_DEBUG << " * peak " << k << " right boundary " << right_borders.back() <<  " with inty " << max_int << std::endl;
        }
      }

      // Return for empty peak list
      if (right_borders.empty())
      {
        return;
      }

      // FEATURE IDEA: instead of Z-score use modified Z-score for small data sets 
      // http://d-scholarship.pitt.edu/7948/1/Seo.pdf
      // http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
      // 1. calculate median
      // 2. MAD = calculate difference to median for each value -> take median of that
      // 3. Mi = 0.6745*(xi - median) / MAD

      // 2. Calculate mean and standard deviation
      // If the coefficient of variation is too large for one border, we use a
      // "pseudo-median" instead of the border of the most intense peak.
      double mean, stdev;

      // Right borders
      mean = std::accumulate(right_borders.begin(), right_borders.end(), 0.0) / (double) right_borders.size();
      stdev = std::sqrt(std::inner_product(right_borders.begin(), right_borders.end(), right_borders.begin(), 0.0)
                               / right_borders.size() - mean * mean);
      std::sort(right_borders.begin(), right_borders.end());

      OPENMS_LOG_DEBUG << " - Recalculating right peak boundaries " << mean << " mean / best " 
                << best_right << " std " << stdev << " : "  << std::fabs(best_right - mean) / stdev 
                << " coefficient of variation" << std::endl;

      // Compare right borders of best transition with the mean
      if (std::fabs(best_right - mean) / stdev > max_z)
      {
        best_right = right_borders[right_borders.size() / 2]; // pseudo median
        OPENMS_LOG_DEBUG << " - CV too large: correcting right boundary to  " << best_right << std::endl;
      }

      // Left borders
      mean = std::accumulate(left_borders.begin(), left_borders.end(), 0.0) / (double) left_borders.size();
      stdev = std::sqrt(std::inner_product(left_borders.begin(), left_borders.end(), left_borders.begin(), 0.0)
                        / left_borders.size() - mean * mean);
      std::sort(left_borders.begin(), left_borders.end());

      OPENMS_LOG_DEBUG << " - Recalculating left peak boundaries " << mean << " mean / best " 
                << best_left << " std " << stdev << " : "  << std::fabs(best_left - mean) / stdev 
                << " coefficient of variation" << std::endl;

      // Compare left borders of best transition with the mean
      if (std::fabs(best_left - mean)  / stdev > max_z)
      {
        best_left = left_borders[left_borders.size() / 2]; // pseudo median
        OPENMS_LOG_DEBUG << " - CV too large: correcting left boundary to  " << best_left << std::endl;
      }

    }

    /// @name Resampling methods
    //@{

    /**
      @brief Create an empty master peak container that has the correct mz / RT values set

      The empty master peak container fill be filled with mz / RT values at the
      positions where the reference chromatogram has values. The container will
      only be populated between the boundaries given. The output container
      will contain peaks with mz / RT values but all intensity values will be zero.

      @param[in] ref_chromatogram Reference chromatogram containing mz / RT values (possibly beyond the desired range)
      @param[out] master_peak_container Output container to be populated
      @param[in] left_boundary Left boundary of values the container should be populated with
      @param[in] right_boundary Right boundary of values the container should be populated with

    */
    template <typename PeakContainerT>
    void prepareMasterContainer_(const PeakContainerT& ref_chromatogram,
                                 PeakContainerT& master_peak_container,
                                 double left_boundary,
                                 double right_boundary)
    {
      OPENMS_PRECONDITION(master_peak_container.empty(), "Master peak container must be empty")

      // get the start / end point of this chromatogram => then add one more
      // point beyond the two boundaries to make the resampling accurate also
      // at the edge.
      auto begin = ref_chromatogram.begin();
      while (begin != ref_chromatogram.end() && begin->getPos() < left_boundary) {begin++; }
      if (begin != ref_chromatogram.begin()) { begin--; }

      auto end = begin;
      while (end != ref_chromatogram.end() && end->getPos() < right_boundary) {end++; }
      if (end != ref_chromatogram.end()) { end++; }

      // resize the master container and set the m/z values to the ones of the master container
      master_peak_container.resize(distance(begin, end)); // initialize to zero
      auto it = master_peak_container.begin();
      for (auto chrom_it = begin; chrom_it != end; chrom_it++, it++)
      {
        it->setPos(chrom_it->getPos());
      }
    }

    /**
      @brief Resample a container at the positions indicated by the master peak container

      @param[in] chromatogram Container with the input data
      @param[in] master_peak_container Container with the mz / RT values at which to resample
      @param[in] left_boundary Left boundary of values the container should be resampled
      @param[in] right_boundary Right boundary of values the container should be resampled

      @return A container which contains the data from the input chromatogram resampled at the positions of the master container
    */
    template <typename PeakContainerT>
    PeakContainerT resampleChromatogram_(const PeakContainerT& chromatogram,
                                    const PeakContainerT& master_peak_container,
                                    double left_boundary,
                                    double right_boundary)
    {
      // get the start / end point of this chromatogram => then add one more
      // point beyond the two boundaries to make the resampling accurate also
      // at the edge.
      auto begin = chromatogram.begin();
      while (begin != chromatogram.end() && begin->getPos() < left_boundary) {begin++;}
      if (begin != chromatogram.begin()) { begin--; }

      auto end = begin;
      while (end != chromatogram.end() && end->getPos() < right_boundary) {end++;}
      if (end != chromatogram.end()) { end++; }

      auto resampled_peak_container = master_peak_container; // copy the master container, which contains the RT values
      LinearResamplerAlign lresampler;
      lresampler.raster(begin, end, resampled_peak_container.begin(), resampled_peak_container.end());

      return resampled_peak_container;
    }

    //@}

    // Members
    String peak_integration_;
    String background_subtraction_;
    bool recalculate_peaks_;
    bool use_precursors_;
    bool use_consensus_;
    bool compute_peak_quality_;
    bool compute_peak_shape_metrics_;
    bool compute_total_mi_;
    double min_qual_;

    int stop_after_feature_;
    double stop_after_intensity_ratio_;
    double min_peak_width_;
    double recalculate_peaks_max_z_;
    double resample_boundary_;

    /**
      @brief Which method to use for selecting peaks' boundaries

      Valid values are: "largest", "widest"
    */
    String boundary_selection_method_;

    PeakPickerChromatogram picker_;
    PeakIntegrator pi_;
  };
}
