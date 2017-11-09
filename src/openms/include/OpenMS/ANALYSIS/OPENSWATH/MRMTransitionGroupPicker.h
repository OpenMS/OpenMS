// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMTRANSITIONGROUPPICKER_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMTRANSITIONGROUPPICKER_H

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

// Cross-correlation
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>

#include <numeric>

//#define DEBUG_TRANSITIONGROUPPICKER

namespace OpenMS
{

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
    ~MRMTransitionGroupPicker();
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

      std::vector<MSChromatogram > picked_chroms_;
      std::vector<MSChromatogram > smoothed_chroms_;

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
        picker_.pickChromatogram(chromatogram, picked_chrom, smoothed_chrom);
        picked_chrom.sortByIntensity();
        picked_chroms_.push_back(picked_chrom);
        smoothed_chroms_.push_back(smoothed_chrom);
      }

      // Pick precursor chromatograms
      if (use_precursors_)
      {
        for (Size k = 0; k < transition_group.getPrecursorChromatograms().size(); k++)
        {
          SpectrumT picked_chrom, smoothed_chrom;
          SpectrumT& chromatogram = transition_group.getPrecursorChromatograms()[k];
          String native_id = chromatogram.getNativeID();

          picker_.pickChromatogram(chromatogram, picked_chrom, smoothed_chrom);
          picked_chrom.sortByIntensity();
          picked_chroms_.push_back(picked_chrom);
          smoothed_chroms_.push_back(smoothed_chrom);
        }
      }

      // Find features (peak groups) in this group of transitions.
      // While there are still peaks left, one will be picked and used to create
      // a feature. Whenever we run out of peaks, we will get -1 back as index
      // and terminate.
      int chr_idx, peak_idx, cnt = 0;
      std::vector<MRMFeature> features;
      while (true)
      {
        chr_idx = -1; peak_idx = -1;
        findLargestPeak(picked_chroms_, chr_idx, peak_idx);
        if (chr_idx == -1 && peak_idx == -1) break;

        // Compute a feature from the individual chromatograms and add non-zero features
        MRMFeature mrm_feature = createMRMFeature(transition_group, picked_chroms_, smoothed_chroms_, chr_idx, peak_idx);
        if (mrm_feature.getIntensity() > 0)
        {
          features.push_back(mrm_feature);
        }

        cnt++;
        if (stop_after_feature_ > 0 && cnt > stop_after_feature_) {break;}
        if (mrm_feature.getIntensity() > 0 && 
            mrm_feature.getIntensity() / (double)mrm_feature.getMetaValue("total_xic") < stop_after_intensity_ratio_)
        {
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
          if ((double)mrm_feature.getMetaValue("leftWidth") >=  (double)features[j].getMetaValue("leftWidth") && 
              (double)mrm_feature.getMetaValue("rightWidth") <= (double)features[j].getMetaValue("rightWidth") )
          { skip = true; }
        }
        if (mrm_feature.getIntensity() > 0 && !skip)
        {
          transition_group.addFeature(mrm_feature);
        }
      }

    }

    /// Create feature from a vector of chromatograms and a specified peak
    template <typename SpectrumT, typename TransitionT>
    MRMFeature createMRMFeature(MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                                std::vector<SpectrumT>& picked_chroms, std::vector<SpectrumT>& smoothed_chroms, const int chr_idx, const int peak_idx)
    {
      OPENMS_PRECONDITION(transition_group.isInternallyConsistent(), "Consistent state required")
      OPENMS_PRECONDITION(transition_group.chromatogramIdsMatch(), "Chromatogram native IDs need to match keys in transition group")

      MRMFeature mrmFeature;
      mrmFeature.setIntensity(0.0);
      double best_left = picked_chroms[chr_idx].getFloatDataArrays()[1][peak_idx];
      double best_right = picked_chroms[chr_idx].getFloatDataArrays()[2][peak_idx];
      double peak_apex = picked_chroms[chr_idx][peak_idx].getRT();
      LOG_DEBUG << "**** Creating MRMFeature for peak " << chr_idx << " " << peak_idx << " " <<
        picked_chroms[chr_idx][peak_idx] << " with borders " << best_left << " " <<
        best_right << " (" << best_right - best_left << ")" << std::endl;

      if (recalculate_peaks_)
      {
        // This may change best_left / best_right
        recalculatePeakBorders_(picked_chroms, best_left, best_right, recalculate_peaks_max_z_);
        if (peak_apex < best_left || peak_apex > best_right)
        {
          // apex fell out of range, lets correct it
          peak_apex = (best_left + best_right) / 2.0;
        }
      }
      picked_chroms[chr_idx][peak_idx].setIntensity(0.0);

      // Remove other, overlapping, picked peaks (in this and other
      // chromatograms) and then ensure that at least one peak is set to zero
      // (the currently best peak).
      remove_overlapping_features(picked_chroms, best_left, best_right);

      // Check for minimal peak width -> return empty feature (Intensity zero)
      if (min_peak_width_ > 0.0 && std::fabs(best_right - best_left) < min_peak_width_) 
      {
        return mrmFeature;
      }

      if (compute_peak_quality_)
      {
        String outlier = "none";
        double qual = computeQuality_(transition_group, picked_chroms, chr_idx, best_left, best_right, outlier);
        if (qual < min_qual_) 
        {
          return mrmFeature;
        }
        mrmFeature.setMetaValue("potentialOutlier", outlier);
        mrmFeature.setMetaValue("initialPeakQuality", qual);
        mrmFeature.setOverallQuality(qual);
      }

      // Prepare linear resampling of all the chromatograms, here creating the
      // empty master_peak_container with the same RT (m/z) values as the reference
      // chromatogram.
      SpectrumT master_peak_container;
      const SpectrumT& ref_chromatogram = selectChromHelper_(transition_group, picked_chroms[chr_idx].getNativeID());
      prepareMasterContainer_(ref_chromatogram, master_peak_container, best_left, best_right);

      // Iterate over initial transitions / chromatograms (note that we may
      // have a different number of picked chromatograms than total transitions
      // as not all are detecting transitions).
      double total_intensity = 0; double total_peak_apices = 0; double total_xic = 0;
      for (Size k = 0; k < transition_group.getTransitions().size(); k++)
      {
        const SpectrumT& chromatogram = selectChromHelper_(transition_group, transition_group.getTransitions()[k].getNativeID()); 
        if (transition_group.getTransitions()[k].isDetectingTransition())
        {
          for (typename SpectrumT::const_iterator it = chromatogram.begin(); it != chromatogram.end(); it++)
          {
            total_xic += it->getIntensity();
          }
        }

        // resample the current chromatogram
        const SpectrumT used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, best_left, best_right);
        // const SpectrumT& used_chromatogram = chromatogram; // instead of resampling

        Feature f;
        double quality = 0;
        f.setQuality(0, quality);
        f.setOverallQuality(quality);

        ConvexHull2D::PointArrayType hull_points;
        double intensity_sum(0.0), intensity_integral(0), rt_sum(0.0);
        double peak_apex_int = -1;
        calculatePeakApexInt_(used_chromatogram, best_left, best_right, hull_points, intensity_sum, intensity_integral, rt_sum, peak_apex_int, peak_apex);

        double background(0), avg_noise_level(0);
        if (background_subtraction_ != "none")
        {
          if ((background_subtraction_ == "smoothed_average" || background_subtraction_ == "smoothed_exact") && smoothed_chroms.size() <= k)
          {
            std::cerr << "Tried to calculate background estimation without any smoothed chromatograms" << std::endl;
            background =  0;
          }
          else if (background_subtraction_ == "smoothed_average")
          {
            calculateBgEstimationAverage_(smoothed_chroms[k], best_left, best_right, background, avg_noise_level);
          }
          else if (background_subtraction_ == "smoothed_exact")
          {
            calculateBgEstimationExact_(smoothed_chroms[k], best_left, best_right, peak_apex_int, background, avg_noise_level);
          }
          else if (background_subtraction_ == "original_average")
          {
            calculateBgEstimationAverage_(used_chromatogram, best_left, best_right, background, avg_noise_level);
          }
          else if (background_subtraction_ == "original_exact")
          {
            calculateBgEstimationExact_(used_chromatogram, best_left, best_right, peak_apex_int, background, avg_noise_level);
          }
          intensity_sum -= background;
          peak_apex_int -= avg_noise_level;
          if (intensity_sum < 0) {intensity_sum = 0;}
          if (peak_apex_int < 0) {peak_apex_int = 0;}
        }

        f.setRT(picked_chroms[chr_idx][peak_idx].getMZ());
        f.setIntensity(intensity_sum);
        ConvexHull2D hull;
        hull.setHullPoints(hull_points);
        f.getConvexHulls().push_back(hull);
        if (chromatogram.metaValueExists("product_mz"))
        {
          f.setMetaValue("MZ", chromatogram.getMetaValue("product_mz"));
          f.setMZ(chromatogram.getMetaValue("product_mz"));
        }
        else
        {
          LOG_WARN << "Please set meta value 'product_mz' on chromatogram to populate feature m/z value" << std::endl;
        }
        f.setMetaValue("native_id", chromatogram.getNativeID());
        f.setMetaValue("peak_apex_int", peak_apex_int);
        if (background_subtraction_ != "none")
        {
          f.setMetaValue("area_background_level", background);
          f.setMetaValue("noise_background_level", avg_noise_level);
        }

        if (transition_group.getTransitions()[k].isDetectingTransition())
        {
          total_intensity += intensity_sum;
          total_peak_apices += peak_apex_int;
        }
    
        if (compute_peak_shape_metrics_)
        { //for backwards compatibility with TOPP tests
          // Calculate peak shape metrics that will be used for later QC  
          PeakShapeMetrics_ peakShapeMetrics;        
          calculatePeakShapeMetrics_(used_chromatogram, 
            best_left, best_right, 
            peak_apex_int, peak_apex, avg_noise_level,
            peakShapeMetrics);    

          f.setMetaValue("width_at_5", peakShapeMetrics.width_at_5);    
          f.setMetaValue("width_at_10", peakShapeMetrics.width_at_10);
          f.setMetaValue("width_at_50", peakShapeMetrics.width_at_50);
          f.setMetaValue("start_time_at_10", peakShapeMetrics.start_time_at_10);
          f.setMetaValue("start_time_at_5", peakShapeMetrics.start_time_at_5);
          f.setMetaValue("end_time_at_10", peakShapeMetrics.end_time_at_10);
          f.setMetaValue("end_time_at_5", peakShapeMetrics.end_time_at_5);
          f.setMetaValue("total_width", peakShapeMetrics.total_width);
          f.setMetaValue("tailing_factor", peakShapeMetrics.tailing_factor);
          f.setMetaValue("asymmetry_factor", peakShapeMetrics.asymmetry_factor);
          f.setMetaValue("baseline_delta_2_height", peakShapeMetrics.baseline_delta_2_height);
          f.setMetaValue("slope_of_baseline", peakShapeMetrics.slope_of_baseline);
          f.setMetaValue("points_across_baseline", peakShapeMetrics.points_across_baseline);
          f.setMetaValue("points_across_half_height", peakShapeMetrics.points_across_half_height);   
        }       

        mrmFeature.addFeature(f, chromatogram.getNativeID()); //map index and feature
      }

      // Also pick the precursor chromatogram(s); note total_xic is not
      // extracted here, only for fragment traces
      for (Size k = 0; k < transition_group.getPrecursorChromatograms().size(); k++)
      {
        const SpectrumT& chromatogram = transition_group.getPrecursorChromatograms()[k];
        // resample the current chromatogram
        const SpectrumT used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, best_left, best_right);

        Feature f;
        double quality = 0;
        f.setQuality(0, quality);
        f.setOverallQuality(quality);

        ConvexHull2D::PointArrayType hull_points;
        double intensity_sum(0.0), rt_sum(0.0);
        double peak_apex_int = -1;
        double peak_apex_dist = std::fabs(used_chromatogram.begin()->getMZ() - peak_apex);
        // FEATURE : use RTBegin / MZBegin -> for this we need to know whether the template param is a real chromatogram or a spectrum!
        for (typename SpectrumT::const_iterator it = used_chromatogram.begin(); it != used_chromatogram.end(); it++)
        {
          if (it->getMZ() > best_left && it->getMZ() < best_right)
          {
            DPosition<2> p;
            p[0] = it->getMZ();
            p[1] = it->getIntensity();
            hull_points.push_back(p);
            if (std::fabs(it->getMZ() - peak_apex) <= peak_apex_dist)
            {
              peak_apex_int = p[1];
              peak_apex_dist = std::fabs(it->getMZ() - peak_apex);
            }
            rt_sum += it->getMZ();
            intensity_sum += it->getIntensity();
          }
        }

        if (chromatogram.metaValueExists("precursor_mz")) 
        {
          f.setMZ(chromatogram.getMetaValue("precursor_mz"));
          mrmFeature.setMZ(chromatogram.getMetaValue("precursor_mz"));
        }

        f.setRT(picked_chroms[chr_idx][peak_idx].getMZ());
        f.setIntensity(intensity_sum);
        ConvexHull2D hull;
        hull.setHullPoints(hull_points);
        f.getConvexHulls().push_back(hull);
        f.setMetaValue("native_id", chromatogram.getNativeID());
        f.setMetaValue("peak_apex_int", peak_apex_int);

        if (use_precursors_ && transition_group.getTransitions().empty())
        {
          total_intensity += intensity_sum;
        }

        mrmFeature.addPrecursorFeature(f, chromatogram.getNativeID());
      }

      mrmFeature.setRT(peak_apex);
      mrmFeature.setIntensity(total_intensity);
      mrmFeature.setMetaValue("PeptideRef", transition_group.getTransitionGroupID());
      mrmFeature.setMetaValue("leftWidth", best_left);
      mrmFeature.setMetaValue("rightWidth", best_right);
      mrmFeature.setMetaValue("total_xic", total_xic);
      mrmFeature.setMetaValue("peak_apices_sum", total_peak_apices);

      mrmFeature.ensureUniqueId();
      return mrmFeature;
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
      //std::cout << "Removing features for peak  between " << best_left << " " << best_right << std::endl;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getMZ() >= best_left && picked_chroms[k][i].getMZ() <= best_right)
          {
            //std::cout << "For Chrom " << k << " removing peak " << picked_chroms[k][i].getMZ() << " l/r : " << picked_chroms[k].getFloatDataArrays()[1][i] << " " <<
            //  picked_chroms[k].getFloatDataArrays()[2][i] << " with int " <<  picked_chroms[k][i].getIntensity() <<std::endl;
            picked_chroms[k][i].setIntensity(0.0);
          }
        }
      }

      // delete all seeds that overlap within the current seed
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getIntensity() <= 0.0) {continue; }

          double left = picked_chroms[k].getFloatDataArrays()[1][i];
          double right = picked_chroms[k].getFloatDataArrays()[2][i];
          if ((left > best_left && left < best_right)
             || (right > best_left && right < best_right))
          {
            //std::cout << "= For Chrom " << k << " removing contained peak " << picked_chroms[k][i].getMZ() << " l/r : " << picked_chroms[k].getFloatDataArrays()[1][i] << " " <<
            //  picked_chroms[k].getFloatDataArrays()[2][i] << " with int " <<  picked_chroms[k][i].getIntensity() <<std::endl;
            picked_chroms[k][i].setIntensity(0.0);
          }
        }
      }
    }

    /// Find largest peak in a vector of chromatograms
    void findLargestPeak(std::vector<MSChromatogram >& picked_chroms, int& chr_idx, int& peak_idx);
    
    /**
      @brief The background noise is estimated based on the peak boundaries 

      The average noise level is computed as the average intensity of the left and
      right peak borders. The background is then computed by multiplying
      the average noise level by all data points under the peak 
      (i.e., integrating the background under the peak).
      The integration strategy assumes that all points are equally spaced.

    */
    void calculateBgEstimationAverage_(const MSChromatogram& chromatogram,
                                  double best_left, double best_right, double & background, double & avg_noise_level);
    
    /**
      @brief The background noise is estimated based on the exact intensity 
        at the left and right peak borders.

      The average noise level is computed by interpolating the intensity at the retention time
      of the peak apex by calculating a line from the intensities at the left and right peak borders. 
      The background is then computed by summing the intensities at each point along
      the line between the intensities at the left and right peak borders 
      (i.e., integrating the background under the peak).
      The integration strategy assumes that all points are equally spaced.

    */
    void calculateBgEstimationExact_(const MSChromatogram& chromatogram,
                                  double best_left, double best_right, double peak_height, double & background, double & avg_noise_level);

    /**
    @brief Will use the chromatogram to get the maximum peak intensity

    The maximum peak intensity/height is calculated.  The convex hull points,
    intensity_sum, and rt_sum are also calculated.
    */
    void calculatePeakApexInt_(const MSChromatogram& chromatogram,
    double best_left, double best_right, 
    ConvexHull2D::PointArrayType & hull_points,
    double & intensity_sum, 
    double & intensity_integral,
    double & rt_sum,
    double & peak_apex_int,
    double & peak_apex_rt);

    // internal structure to represent various peak shape metrics
    struct PeakShapeMetrics_ {
      double width_at_5 = 0.0;
      double width_at_10 = 0.0;
      double width_at_50 = 0.0;
      double start_time_at_10 = 0.0;
      double start_time_at_5 = 0.0;
      double start_time_at_50 = 0.0;
      double end_time_at_10 = 0.0;
      double end_time_at_5 = 0.0;
      double end_time_at_50 = 0.0;
      double total_width = 0.0;
      double tailing_factor = 0.0; /**The tailing factor is a measure of peak tailing. 
        It is defined as the distance from the front slope of the peak to the back slope 
        divided by twice the distance from the center line of the peak to the front slope, 
        with all measurements made at 5% of the maximum peak height.  
        tailing_factor = Tf = W0.05/2a
          where W0.05 is peak width at 5% max peak height
          a = min width to peak maximum at 5% max peak height
          b = max width to peak maximum at 5% max peak height
          0.9 < Tf < 1.2
          front Tf < 0.9
          tailing Tf > 1.2*/       
      double asymmetry_factor = 0.0; /**The asymmetry factor is a measure of peak tailing. 
        It is defined as the distance from the center line of the peak to the back slope 
        divided by the distance from the center line of the peak to the front slope, 
        with all measurements made at 10% of the maximum peak height. 
        asymmetry_factor = As = b/a
		      where a is min width to peak maximum at 10% max peak height
		      b is max width to peak maximum at 10% max peak height */
      double baseline_delta_2_height = 0.0; /**The change in baseline divided by the height is
        a way of comparing the influence of the change of baseline on the peak height.*/
      double slope_of_baseline = 0.0; /**The slope of the baseline is a measure of slope change.
        It is approximated as the difference in baselines between the peak start and peak end.  */
      int points_across_baseline = 0;
      int points_across_half_height = 0;
    };

    /**
    @brief Calculates standard peak shape quality metrics

    Standard peak shape quality metrics are calculated for down stream QC/QA.
    */
    void calculatePeakShapeMetrics_(const MSChromatogram& chromatogram, 
    double best_left, double best_right, 
    double peak_height, double peak_apex_rt, double avg_noise_level,
    PeakShapeMetrics_ & peakShapeMetrics);
    
protected:

    /// Synchronize members with param class
    void updateMembers_();

    /// Assignment operator is protected for algorithm
    MRMTransitionGroupPicker& operator=(const MRMTransitionGroupPicker& rhs);

    /**
      @brief Select matching precursor or fragment ion chromatogram
    */
    template <typename SpectrumT, typename TransitionT>
    const SpectrumT& selectChromHelper_(MRMTransitionGroup<SpectrumT, TransitionT>& transition_group, String native_id)
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
    double computeQuality_(MRMTransitionGroup<SpectrumT, TransitionT>& transition_group,
                           std::vector<SpectrumT>& picked_chroms, const int chr_idx,
                           const double best_left, const double best_right, String& outlier)
    {

      // Resample all chromatograms around the current estimated peak and
      // collect the raw intensities. For resampling, use a bit more on either
      // side to correctly identify shoulders etc.
      double resample_boundary = resample_boundary_; // sample 15 seconds more on each side
      SpectrumT master_peak_container;
      const SpectrumT& ref_chromatogram = selectChromHelper_(transition_group, picked_chroms[chr_idx].getNativeID());
      prepareMasterContainer_(ref_chromatogram, master_peak_container, best_left - resample_boundary, best_right + resample_boundary);
      std::vector<std::vector<double> > all_ints;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        const SpectrumT chromatogram = selectChromHelper_(transition_group, picked_chroms[k].getNativeID());
        const SpectrumT used_chromatogram = resampleChromatogram_(chromatogram, 
            master_peak_container, best_left - resample_boundary, best_right + resample_boundary);

        std::vector<double> int_here;
        for (Size i = 0; i < used_chromatogram.size(); i++)
        {
          int_here.push_back(used_chromatogram[i].getIntensity());
        }
        all_ints.push_back(int_here);
      }

      // Compute the cross-correlation for the collected intensities
      std::vector<double> mean_shapes;
      std::vector<double> mean_coel;
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
      std::vector<double> left_borders;
      std::vector<double> right_borders;
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
          if (picked_chroms[k][i].getMZ() >= best_left && picked_chroms[k][i].getMZ() <= best_right)
          {
            pfound++;
            if (picked_chroms[k][i].getIntensity() > max_int)
            {
              max_int = picked_chroms[k][i].getIntensity() > max_int;
              l_tmp = picked_chroms[k].getFloatDataArrays()[1][i];
              r_tmp = picked_chroms[k].getFloatDataArrays()[2][i];
            }
          }
        }

        if (l_tmp > 0.0) left_borders.push_back(l_tmp);
        if (r_tmp > 0.0) right_borders.push_back(r_tmp);

        if (pfound == 0) missing_peaks++;
        if (pfound > 1) multiple_peaks++;
      }

      // Check how many chromatograms had exactly one peak picked between our
      // current left/right borders -> this would be a sign of consistency.
      LOG_DEBUG << " Overall found missing : " << missing_peaks << " and multiple : " << multiple_peaks << std::endl;

      /// left_borders / right_borders might not have the same length since we might have peaks missing!!

      // Is there one transitions that is very different from the rest (e.g.
      // the same element has a bad shape and a bad coelution score) -> potential outlier
      if (min_index_shape == max_index_coel)
      {
        LOG_DEBUG << " element " << min_index_shape << " is a candidate for removal ... " << std::endl;
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

      LOG_DEBUG << " computed score  " << score << " (from " <<  shape_score << 
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
    void recalculatePeakBorders_(std::vector<SpectrumT>& picked_chroms, double& best_left, double& best_right, double max_z)
    {
      // 1. Collect all seeds that lie within the current seed 
      // - Per chromatogram only the most intense one counts, otherwise very
      // - low intense peaks can contribute disproportionally to the voting
      // - procedure.
      std::vector<double> left_borders;
      std::vector<double> right_borders;
      for (Size k = 0; k < picked_chroms.size(); k++)
      {
        double max_int = -1;
        double left = -1;
        double right = -1;
        for (Size i = 0; i < picked_chroms[k].size(); i++)
        {
          if (picked_chroms[k][i].getMZ() >= best_left && picked_chroms[k][i].getMZ() <= best_right)
          {
            if (picked_chroms[k].getFloatDataArrays()[0][i] > max_int)
            {
              max_int = picked_chroms[k].getFloatDataArrays()[0][i];
              left = picked_chroms[k].getFloatDataArrays()[1][i];
              right = picked_chroms[k].getFloatDataArrays()[2][i];
            }
          }
        }
        if (max_int > -1 )
        {
          left_borders.push_back(left);
          right_borders.push_back(right);
          LOG_DEBUG << " * " << k << " left boundary " << left_borders.back()   <<  " with int " << max_int << std::endl;
          LOG_DEBUG << " * " << k << " right boundary " << right_borders.back() <<  " with int " << max_int << std::endl;
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

      LOG_DEBUG << " - Recalculating right peak boundaries " << mean << " mean / best " 
                << best_right << " std " << stdev << " : "  << std::fabs(best_right - mean) / stdev 
                << " coefficient of variation" << std::endl;

      // Compare right borders of best transition with the mean
      if (std::fabs(best_right - mean) / stdev > max_z)
      {
        best_right = right_borders[right_borders.size() / 2]; // pseudo median
        LOG_DEBUG << " - Setting right boundary to  " << best_right << std::endl;
      }

      // Left borders
      mean = std::accumulate(left_borders.begin(), left_borders.end(), 0.0) / (double) left_borders.size();
      stdev = std::sqrt(std::inner_product(left_borders.begin(), left_borders.end(), left_borders.begin(), 0.0)
                        / left_borders.size() - mean * mean);
      std::sort(left_borders.begin(), left_borders.end());

      LOG_DEBUG << " - Recalculating left peak boundaries " << mean << " mean / best " 
                << best_left << " std " << stdev << " : "  << std::fabs(best_left - mean) / stdev 
                << " coefficient of variation" << std::endl;

      // Compare left borders of best transition with the mean
      if (std::fabs(best_left - mean)  / stdev > max_z)
      {
        best_left = left_borders[left_borders.size() / 2]; // pseudo median
        LOG_DEBUG << " - Setting left boundary to  " << best_left << std::endl;
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

      @param ref_chromatogram Reference chromatogram containing mz / RT values (possibly beyond the desired range)
      @param master_peak_container Output container to be populated
      @param left_boundary Left boundary of values the container should be populated with
      @param right_boundary Right boundary of values the container should be populated with

    */
    template <typename SpectrumT>
    void prepareMasterContainer_(const SpectrumT& ref_chromatogram,
                                 SpectrumT& master_peak_container, double left_boundary, double right_boundary)
    {
      OPENMS_PRECONDITION(master_peak_container.empty(), "Master peak container must be empty")

      // get the start / end point of this chromatogram => then add one more
      // point beyond the two boundaries to make the resampling accurate also
      // at the edge.
      typename SpectrumT::const_iterator begin = ref_chromatogram.begin();
      while (begin != ref_chromatogram.end() && begin->getMZ() < left_boundary) {begin++; }
      if (begin != ref_chromatogram.begin()) {begin--; }

      typename SpectrumT::const_iterator end = begin;
      while (end != ref_chromatogram.end() && end->getMZ() < right_boundary) {end++; }
      if (end != ref_chromatogram.end()) {end++; }

      // resize the master container and set the m/z values to the ones of the master container
      master_peak_container.resize(distance(begin, end)); // initialize to zero
      typename SpectrumT::iterator it = master_peak_container.begin();
      for (typename SpectrumT::const_iterator chrom_it = begin; chrom_it != end; chrom_it++, it++)
      {
        it->setMZ(chrom_it->getMZ());
      }
    }

    /**
      @brief Resample a container at the positions indicated by the master peak container

      @param chromatogram Container with the input data
      @param master_peak_container Container with the mz / RT values at which to resample
      @param left_boundary Left boundary of values the container should be resampled
      @param right_boundary Right boundary of values the container should be resampled

      @return A container which contains the data from the input chromatogram resampled at the positions of the master container
    */
    template <typename SpectrumT>
    SpectrumT resampleChromatogram_(const SpectrumT& chromatogram,
                                    const SpectrumT& master_peak_container, double left_boundary, double right_boundary)
    {
      // get the start / end point of this chromatogram => then add one more
      // point beyond the two boundaries to make the resampling accurate also
      // at the edge.
      typename SpectrumT::const_iterator begin = chromatogram.begin();
      while (begin != chromatogram.end() && begin->getMZ() < left_boundary) {begin++;}
      if (begin != chromatogram.begin()) {begin--;}

      typename SpectrumT::const_iterator end = begin;
      while (end != chromatogram.end() && end->getMZ() < right_boundary) {end++;}
      if (end != chromatogram.end()) {end++;}

      SpectrumT resampled_peak_container = master_peak_container; // copy the master container, which contains the RT values
      LinearResamplerAlign lresampler;
      lresampler.raster(begin, end, resampled_peak_container.begin(), resampled_peak_container.end());

      return resampled_peak_container;
    }

    //@}

    // Members
    String background_subtraction_;
    bool recalculate_peaks_;
    bool use_precursors_;
    bool compute_peak_quality_;
    bool compute_peak_shape_metrics_;
    double min_qual_;

    int stop_after_feature_;
    double stop_after_intensity_ratio_;
    double min_peak_width_;
    double recalculate_peaks_max_z_;
    double resample_boundary_;

    PeakPickerMRM picker_;
  };
}

#endif //  OPENMS_ANALYSIS_OPENSWATH_MRMTRANSITIONGROUPPICKER_H

