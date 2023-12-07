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
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>


#include <vector>
#include <functional>


namespace OpenMS
{
  class OPENMS_DLLAPI TagGenerator
  {
  private:
    MSSpectrum spectrum_;
    std::vector<bool> selected_peaks_;                   // The peaks we actually want to use
    uint32_t n;                                          // the number of globally selected peaks;
    std::vector<std::shared_ptr<TagGeneratorNode>> dag_; // directed acyclic graph containing all peaks



  public:

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

    class MultiPeak
    {
    public:
      MultiPeak();

      MultiPeak(Peak1D peak, float score);

      /// Copy
      MultiPeak(const MultiPeak& other);
      /// Assignment
      MultiPeak& operator=(const MultiPeak& other);
      /// Destructor
      virtual ~MultiPeak() = default;

      [[nodiscard]] const Peak1D& getPeak() const;
      float getScore() const;
      const std::string& getFollowUpPeaksAa() const;
      const std::vector<float>& getFollowUpPeaks() const;

      void addFollowUpPeak(float distance, const std::string& AA);
      void addScore(float score);

    protected:
      Peak1D peak_;
      float score_;
      std::string follow_up_peaks_AA;
      std::vector<float> follow_up_peaks;
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

    /// Constructor
    TagGenerator(const MSSpectrum& spectrum);

    /// copy constructor
    TagGenerator(const TagGenerator& cp);

    /// assignemnt operator
    TagGenerator& operator=(const TagGenerator& source);

    /// Destructor
    ~TagGenerator();

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
     * @brief Gets from the peaks of the spectrum all possible nodes for the dag(directed acyclic graph)
     * @param max_charge : The maximal charge a peak should have. Runtime scales with max_charge, but is not the bottleneck
     */
    void generateAllNodes(uint32_t max_charge);

    /**
     * @brief Generates a directed acyclic graph from all selected peaks. Nodes are connected
     * if they have a distance equal to the monoisotopic mass of any Amino Acid
     * @param fragment_tolerance : The allowed tolerance between peak distance and any monoisotopic mass
     */
    void generateDirectedAcyclicGraph(float fragment_tolerance);

    /**
     * @brief Given the dag, generates all possible Multi Peaks and calculates their score.
     * This method is to be used on EXPERIMENTAL SPECTRA
     * @param quad_peaks : vector the Multipeaks are stored in
     * @param depth : The number of follow up peaks one peak should include
     */
    void generateAllMultiPeaks(std::vector<MultiPeak>& quad_peaks, size_t depth);

    /**
     * @brief Given a THEORETICAL SPECTRA compute all multifragments. The spectra must contain meta data about the iontypes
     * @param multi_frags : output
     * @param depth : Number of follow up peaks/fragments one peak should include
     * @param peptide_idx : The index of the peptide, currently processed
     * @param frag_min_mz : The min frag mz that is included
     * @param frag_max_mz
     */
    void generateAllMultiFragments(std::vector<MultiFragment>& multi_frags, size_t depth, size_t peptide_idx, float frag_min_mz, float frag_max_mz);

  };
}