// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>



#include <vector>
#include <functional>


namespace OpenMS
{
  class OPENMS_DLLAPI FragmentIndexTagGenerator
  {
  private:
    MSSpectrum spectrum_;
    std::vector<bool> selected_peaks_;                   // The peaks we actually want to use
    uint32_t n;                                          // the number of globally selected peaks;




  public:

    class MultiPeak
    {
    public:
      MultiPeak();

      MultiPeak(float peak, float score);

      /// Copy
      MultiPeak(const MultiPeak& other);
      /// Assignment
      MultiPeak& operator=(const MultiPeak& other);
      /// Destructor
      virtual ~MultiPeak() = default;

      float getPeak() const;
      float getScore() const;
      const std::string& getFollowUpPeaksAa() const;
      const std::vector<float>& getFollowUpPeaks() const;

      void addFollowUpPeak(float distance, const std::string& AA);
      void addScore(float score);

    protected:
      float peak_;
      float score_;
      std::string follow_up_peaks_AA;       //TODO:: Only for debugging
      std::vector<float> follow_up_peaks;
    };

    class MultiFragment
    {
    public:
      MultiFragment();

      MultiFragment(UInt32 peptide_idx,
                    float fragment_mz,
                    const std::vector<float>& follow_up);

      MultiFragment(UInt32 peptide_idx, float fragment_mz, const MultiPeak& multiPeak);

      MultiFragment(const MultiFragment& other);

      /// Assignment operator
      MultiFragment& operator=(const MultiFragment& other);

      ///Destructor
      virtual ~MultiFragment() = default;

      /// ValueSwappable
      void swap(MultiFragment& other);


      UInt32 getPeptideIdx() const;
      float getFragmentMz() const;
      //const std::string& getFollowUpPeaksAa() const;
      const std::vector<float>& getFollowUpPeaks() const;


    protected:
      UInt32 peptide_idx_;
      float fragment_mz_;
      std::vector<float> follow_up_peaks_;
    };


    /**
     * @brief A struct class containing the intensity and the idx of a peak.
     * Used for simple global and local selection
     */
    struct IdxAndIntensity{
      uint32_t idx_;
      float intensity_;
      IdxAndIntensity(uint32_t idx, float intensity)
      {
        idx_ = idx;
        intensity_ = intensity;
      }
    };
  public:
    /// Constructor
    explicit FragmentIndexTagGenerator(const MSSpectrum& spectrum);

    /// copy constructor
    FragmentIndexTagGenerator(const FragmentIndexTagGenerator& cp);

    /// assignemnt operator
    FragmentIndexTagGenerator& operator=(const FragmentIndexTagGenerator& source);

    /// Destructor
    ~FragmentIndexTagGenerator();

    /// setter
    void setMSSpectrum(const MSSpectrum &spectrum);


    /**@brief top N peaks are selected according to their intensities over the entire m/z range of a spectrum where N is
     *related to a precursor ion mass
     *
     */
    void globalSelection();

    /**@brief peaks are selected by sliding a window of 70 Da (window increment, 35 Da)
     * when fewer than two peaks are selected in any window during the global selection
     *
     */
    void localSelection();


    /**
     * @brief Skipes dag generation and creates all MultiPeaks directly from the Spectrum!
     * Assumes that the spectrum only has charge 1 peaks.
     * @param[out] all_multi_peaks ouput vector containing all the MultiPeaks
     * @param depth The depth of the MultiPeak (= # of follow up peaks)
     */
    void generateAllMultiPeaksFast(std::vector<MultiPeak>& all_multi_peaks,
                                   size_t depth,
                                   float tolerance,
                                   bool tolerance_unit);


    /**
     * @brief Given a THEORETICAL SPECTRA compute all multifragments. The spectra must contain meta data about the iontypes
     * @param multi_frags : output
     * @param depth : Number of follow up peaks/fragments one peak should include
     * @param peptide_idx : The index of the peptide, currently processed
     * @param frag_min_mz : The min frag mz that is included
     * @param frag_max_mz
     */
    void generateAllMultiFragments(std::vector<MultiFragment>& multi_frags, size_t depth, size_t peptide_idx, float frag_min_mz, float frag_max_mz);

  private:
    /**
     * @brief Recursive chain which builds all peaks
     * @param all_multi_peaks[out] When recursion_step reaches 1 and we found a multi_peak it is placed here
     * @param multi_peak COPY of the current multi_peak
     * @param recursion_step
     */
    void generateAllMultiPeaksFastRecursion(std::vector<MultiPeak>& all_multi_peaks, FragmentIndexTagGenerator::MultiPeak& multi_peak, size_t current_peak_idx, size_t recursion_step, float tolerance);

  };
}