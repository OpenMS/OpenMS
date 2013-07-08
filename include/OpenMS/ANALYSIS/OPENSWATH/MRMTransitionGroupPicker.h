// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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


#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>

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

    // this is the type in which we store the chromatograms for this analysis
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; 

    //@{
    /// Constructor
    MRMTransitionGroupPicker();

    /// Destructor
    ~MRMTransitionGroupPicker();
    //@}

    /**
      @brief Pick a group of chromatograms belonging to the same peptide

      Will identify peaks in a set of chromatograms that belong to the same
      peptide. The chromatograms are given inthe MRMTransitionGroup container
      which also contains the mapping of the chromatograms to their metadata.

      The resulting features are added added to the MRMTransitionGroup. Each feature contains the following meta-data:

      - PeptideRef
      - leftWidth
      - rightWidth
      - total_xic
      - peak_apices_sum

    */
    template <typename SpectrumT, typename TransitionT>
    void pickTransitionGroup(MRMTransitionGroup<SpectrumT, TransitionT> & transition_group)
    {
      std::vector<RichPeakChromatogram> picked_chroms_;
      std::vector<RichPeakChromatogram> smoothed_chroms_;

      PeakPickerMRM picker;
      picker.setParameters(param_.copy("PeakPickerMRM:", true));

      // Pick chromatograms
      for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
      {
        RichPeakChromatogram& chromatogram = transition_group.getChromatograms()[k];
        if (!chromatogram.isSorted()) 
        { 
          chromatogram.sortByPosition(); 
        }

        RichPeakChromatogram picked_chrom, smoothed_chrom;
        picker.pickChromatogram(chromatogram, smoothed_chrom, picked_chrom);
        picked_chrom.sortByIntensity(); // we could do without that
        picked_chroms_.push_back(picked_chrom);
        smoothed_chroms_.push_back(smoothed_chrom);
      }

      // Find features (peak groups) in this group of transitions.
      // While there are still peaks left, one will be picked and used to create
      // a feature. Whenever we run out of peaks, we will get -1 back as index
      // and terminate.
      int chr_idx, peak_idx, cnt = 0;
      while (true)
      {
        chr_idx = -1; peak_idx = -1;
        findLargestPeak(picked_chroms_, chr_idx, peak_idx);
        if (chr_idx == -1 && peak_idx == -1) break;

        // get feature, prevent non-extended zero features to be added
        MRMFeature mrm_feature = createMRMFeature(transition_group, picked_chroms_, smoothed_chroms_, chr_idx, peak_idx);
        if (mrm_feature.getIntensity() > 0)
        {
          transition_group.addFeature(mrm_feature);
        }

        cnt++;
        if ((stop_after_feature_ > 0 && cnt > stop_after_feature_) && mrm_feature.getIntensity() / (double)mrm_feature.getMetaValue("total_xic") < stop_after_intensity_ratio_)
        {
          break;
        }
      }
    }

    /// Create feature from a vector of chromatograms and a specified peak
    template <typename SpectrumT, typename TransitionT>
    MRMFeature createMRMFeature(MRMTransitionGroup<SpectrumT, TransitionT> & transition_group,
      std::vector<SpectrumT> & picked_chroms, std::vector<SpectrumT> & smoothed_chroms_, int & chr_idx, int & peak_idx)
    {
      MRMFeature mrmFeature;
      mrmFeature.setIntensity(0.0);
      double best_left = picked_chroms[chr_idx].getFloatDataArrays()[1][peak_idx];
      double best_right = picked_chroms[chr_idx].getFloatDataArrays()[2][peak_idx];
      const double peak_apex = picked_chroms[chr_idx][peak_idx].getRT();

      // Remove other, overlapping, picked peaks (in this and other
      // chromatograms) and then ensure that at least one peak is set to zero
      // (the currently best peak).
      remove_overlapping_features(picked_chroms, best_left, best_right);
      picked_chroms[chr_idx][peak_idx].setIntensity(0.0);

      // Check for minimal peak width
      if( min_peak_width_ > 0.0 && std::fabs(best_right-best_left) < min_peak_width_) {return mrmFeature;}

      // Prepare linear resampling of all the chromatograms, here creating the
      // empty master_peak_container with the same RT (m/z) values as the reference
      // chromatogram.
      SpectrumT master_peak_container;
      const SpectrumT & ref_chromatogram = transition_group.getChromatograms()[chr_idx];
      prepareMasterContainer_(ref_chromatogram, master_peak_container, best_left, best_right);

      double total_intensity = 0; double total_peak_apices = 0; double total_xic = 0;
      for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
      {
        const SpectrumT& chromatogram = transition_group.getChromatograms()[k];
        for (typename SpectrumT::const_iterator it = chromatogram.begin(); it != chromatogram.end(); it++)
        {
          total_xic += it->getIntensity();
        }

        // resample the current chromatogram
        const SpectrumT used_chromatogram = resampleChromatogram_(chromatogram, master_peak_container, best_left, best_right);
        // const SpectrumT& used_chromatogram = chromatogram; // instead of resampling

        Feature f;
        double quality = 0;
        f.setQuality(0, quality);
        f.setOverallQuality(quality);

        ConvexHull2D::PointArrayType hull_points;
        DoubleReal intensity_sum(0.0), rt_sum(0.0);
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

        if (background_subtraction_ != "none")
        {
          double background = 0;
          // we use the smoothed chromatogram here to have a more accurate estimatation of the noise at the flanks of the peak
          if (background_subtraction_ == "smoothed")
          {
            if (smoothed_chroms_.size() <= k)
            {
              std::cerr << "Tried to calculate background estimation without any smoothed chromatograms" << std::endl;
              background =  0;
            }
            else
            {
              background = calculateBgEstimation_(smoothed_chroms_[k], best_left, best_right);
            }
          }
          else if (background_subtraction_ == "original")
          {
            background = calculateBgEstimation_(used_chromatogram, best_left, best_right);
          }
          intensity_sum -= background;
          if (intensity_sum < 0)
          {
            std::cerr << "Warning: Intensity was below 0 after background subtraction: " << intensity_sum << ". Setting it to 0." << std::endl;
            intensity_sum = 0;
          }
        }

        f.setRT(picked_chroms[chr_idx][peak_idx].getMZ());
        f.setMZ(chromatogram.getMetaValue("product_mz"));
        f.setIntensity(intensity_sum);
        ConvexHull2D hull;
        hull.setHullPoints(hull_points);
        f.getConvexHulls().push_back(hull);
        f.setMetaValue("MZ", chromatogram.getMetaValue("product_mz"));
        f.setMetaValue("native_id", chromatogram.getNativeID());
        f.setMetaValue("peak_apex_int", peak_apex_int);
        //f.setMetaValue("leftWidth", best_left);
        //f.setMetaValue("rightWidth", best_right);

        total_intensity += intensity_sum;
        total_peak_apices += peak_apex_int;
        mrmFeature.addFeature(f, chromatogram.getNativeID()); //map index and feature
      }
      mrmFeature.setRT(picked_chroms[chr_idx][peak_idx].getMZ());
      mrmFeature.setIntensity(total_intensity);
      mrmFeature.setMetaValue("PeptideRef", transition_group.getTransitionGroupID());
      mrmFeature.setMetaValue("leftWidth", best_left);
      mrmFeature.setMetaValue("rightWidth", best_right);
      mrmFeature.setMetaValue("total_xic", total_xic);
      mrmFeature.setMetaValue("peak_apices_sum", total_peak_apices);

      return mrmFeature;
    }

    // maybe private, but we have tests
    /**
      @brief Remove overlaping features.
      
      Remove features that are within the current seed (between best_left and
      best_right) or overlap with it. An overlapping feature is defined as a
      feature that has either of its borders within the border of the current
      peak 
      
      Directly adjacent features are allowed, e.g. they can share one
      border.

    */
    template <typename SpectrumT>
    void remove_overlapping_features(std::vector<SpectrumT> & picked_chroms, double best_left, double best_right)
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
    void findLargestPeak(std::vector<RichPeakChromatogram>& picked_chroms, int& chr_idx, int& peak_idx);

protected:

    /// Synchronize members with param class
    void updateMembers_();

    /// Assignment operator is protected for algorithm
    MRMTransitionGroupPicker& operator=(const MRMTransitionGroupPicker& rhs);

    /// @name Resampling methods
    //@{
    /// create an empty master peak container that has the correct mz / RT values set
    template <typename SpectrumT>
    void prepareMasterContainer_(const SpectrumT& ref_chromatogram,
      SpectrumT & master_peak_container, double best_left, double best_right)
    {
      // search for begin / end of the reference chromatogram (and add one more point)
      typename SpectrumT::const_iterator begin = ref_chromatogram.begin();
      while (begin != ref_chromatogram.end() && begin->getMZ() < best_left) {begin++; }
      if (begin != ref_chromatogram.begin()) {begin--; }

      typename SpectrumT::const_iterator end = begin;
      while (end != ref_chromatogram.end() && end->getMZ() < best_right) {end++; }
      if (end != ref_chromatogram.end()) {end++; }

      // resize the master container and set the m/z values to the ones of the master container
      master_peak_container.resize(distance(begin, end));
      typename SpectrumT::iterator it = master_peak_container.begin();
      for (typename SpectrumT::const_iterator chrom_it = begin; chrom_it != end; chrom_it++, it++)
      {
        it->setMZ(chrom_it->getMZ());
      }
    }

    /// use the master container from above to resample a chromatogram at those points stored in the master container
    template <typename SpectrumT>
    SpectrumT resampleChromatogram_(const SpectrumT & chromatogram,
      SpectrumT & master_peak_container, double best_left, double best_right)
    {
      // get the start / end point of this chromatogram => go one past
      // best_left / best_right to make the resampling accurate also at the
      // edge.
      typename SpectrumT::const_iterator begin = chromatogram.begin();
      while (begin != chromatogram.end() && begin->getMZ() < best_left) {begin++; }
      if (begin != chromatogram.begin()) {begin--; }

      typename SpectrumT::const_iterator end = begin;
      while (end != chromatogram.end() && end->getMZ() < best_right) {end++; }
      if (end != chromatogram.end()) {end++; }

      SpectrumT resampled_peak_container = master_peak_container; // copy the master container, which contains the RT values
      LinearResamplerAlign lresampler;
      lresampler.raster(begin, end, resampled_peak_container.begin(), resampled_peak_container.end());

      return resampled_peak_container;
    }

    //@}

    /**
      @brief Will use the smoothed chromatograms to estimate the background noise and then subtract it

      The background is estimated by averaging the noise on either side of the
      peak and then subtracting that from the total intensity.
    */
    double calculateBgEstimation_(const RichPeakChromatogram& smoothed_chromat, double best_left, double best_right);

    // Members
    String background_subtraction_;

    int stop_after_feature_;
    DoubleReal stop_after_intensity_ratio_;
    DoubleReal min_peak_width_;
  };
}

#endif
