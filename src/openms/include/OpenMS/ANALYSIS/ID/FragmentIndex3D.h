// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGeneratorNode.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/MultiFragment.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <functional>
#include <utility>
#include <vector>


namespace OpenMS
{
  /** @brief Generates from a set of Fasta files a 3D-datastructure which stores all theoretical masses of all
   * b and y ions from all peptides generated from the Fasta file. The datastructure is build such that on one axis
   * the fragments are sorted by their own mass and the axis by the mass of their precursor/protein.
   * In addition a third axis is introduced where the Fragments are sorted based on neighbor infomation.
   * The FI has two options: Bottom-up and Top Down. In later digestion is skiped and the fragments have a direct
   * reference to the mass of the proteins instead of digested peptides.
   */
  class OPENMS_DLLAPI FragmentIndex3D : public DefaultParamHandler
  {
  public:
    // STRUCTS
    struct Peptide {
      UInt32 protein_idx_;            ///< The index in fasta entries vector
      float precursor_mz;                  ///< mz of the peptide
    };

    /**
     * @brief Every potential Peptide/Protein has such an struct. Inside the number of peaks-to-Fragment hits are safed
     */
    struct SpectrumMatch
    {
      uint32_t num_matched_{};       ///< Number of peaks-fragment hits
      uint16_t precursor_charge_{};  ///< The precursor_charged used for the performed search
      int16_t isotope_error_{};      /// < The isotope_error used for the performed search
      size_t peptide_idx_{};         ///< The idx this struct belongs to
    };

    struct SpectrumMatchesTopN
    {
      std::vector<SpectrumMatch> hits_;     ///< The preliminary candidates
      uint32_t matched_peaks_{};      ///< The number of matched peaks TODO: statistic needed?
      uint32_t scored_candidates_{};  ///< The number of scored candidates

      SpectrumMatchesTopN() = default;

      /**
       * @brief Appends the a SpectrumMatchesTopN to another one. Add the number of all matched peaks up. Same for number of scored candidates
       * The
       * @param other The appended struct
       * @return The struct after the attachment
       */
      SpectrumMatchesTopN& operator+=(const SpectrumMatchesTopN& other)
      {
        this->matched_peaks_ += other.matched_peaks_;
        this->scored_candidates_ += other.scored_candidates_;
        this->hits_.insert(this->hits_.end(), other.hits_.begin(), other.hits_.end());
        return *this;
      }
    };

    struct Hit
    {
      Hit(UInt32 peptide_idx, float fragment_mz) :
          peptide_idx(peptide_idx),
          fragment_mz(fragment_mz)
      {}
      UInt32 peptide_idx; // index in database
      float fragment_mz;
    };


    /**
     * Builds the Index Database in a multi-dimensional/multi-leveled tree structure
     * @param fasta_entries
     */
    void build(const std::vector<FASTAFile::FASTAEntry> &fasta_entries);

    std::pair<size_t, size_t > getPeptidesInPrecursorRange(float precursor_mass, std::pair<float, float> window);

    /**
     * Takes a MultiPeak as input and recursively searches for a hit in the DB
     * @param hits output
     * @param peak the query Multipeak
     * @param peptide_idx_range The range of all possible Peptides
     * @param window The window in which we want to search, enabeling finding of modified peptides
     */
    void query(std::vector<Hit>& hits, const MultiPeak& peak, std::pair<size_t, size_t> peptide_idx_range, std::pair<float, float> window);

    /// DefaultConstructor
    FragmentIndex3D() ;

    /// Default destructor
    //FragmentIndex;

  protected:

    void updateMembers_() override;

  private:

    void generate_peptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

    /**
     * @brief Check if potential hit is in the range of [query-tolerance + window, query + tolerance + window]
     * @param hit entry in the DB which could be a potential hit
     * @param query the current query
     * @param tolerance
     * @param window
     * @return
     */
    static bool inRange(float hit, float query, float tolerance, std::pair<float, float> window);
    /**
     * @brief compares two vectors and checks if each single element of the hit vector is in range of its corresponding query counterpart
     * @param hit
     * @param query
     * @param tolerance
     * @return
     */
    static bool inRangeFollowUpPeaks(std::vector<float> hit, std::vector<float> query, float tolerance);



    void recursiveQuery(std::vector<Hit>& hits,
                        const MultiPeak& peak,
                        std::pair<size_t, size_t> peptide_idx_range,
                        std::pair<float, float> window,  // The window for the Fragment mass ONLY!!!
                        size_t recursion_step,
                        size_t current_slice,
                        float fragment_tolerance);

    bool is_build_;

    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;

    float fragment_min_mz_;  ///< smallest fragment mz
    float fragment_max_mz_;  ///< largest fragment mz

    float precursor_mz_tolerance_;
    bool precursor_mz_tolerance_unit_ppm_{true};
    float fragment_mz_tolerance_;
    bool fragment_mz_tolerance_unit_ppm_{true};

    uint16_t depth_; // The depth of the database (e.q. Depth 3. We include the next three peaks on the right. The database is then (3+2) Dimensional)
    std::vector<Peptide> fi_peptides_;   ///< vector of all (digested) peptides
    std::vector<MultiFragment> fi_fragments_;
    std::vector<std::vector<float>> follow_up_peaks_buckets_min_mz;
    size_t bucketsize_;       ///< number of fragments per outer node
    std::vector<float> bucket_min_mz_;  ///< vector of the smalles fragment mz of each bucket

    // Search Related:
    uint16_t min_matched_peaks_;  ///< PSM with less hits are discarded
    int16_t min_precursor_charge_; ///< minimal possible precursor charge (usually always 1)
    uint16_t max_precursor_charge_; ///< maximal possible precursor charge
    uint16_t max_fragment_charge_;  ///< The maximal possible charge of the fragments
    uint32_t max_processed_hits_;   ///< The amount of PSM that will be used. the rest is filtered out
    float open_precursor_window_lower; ///< Defines the lower bound of the precursor-mass range
    float open_precursor_window_upper; ///< Defines the upper bound of the precursor-mass range
    float open_fragment_window_lower;  ///< Defines the lower bound of the fragment-mass range
    float open_fragment_window_upper;  ///< Defines the upper bound of the fragment-mass range
  };

}