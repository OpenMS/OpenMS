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
  /** @brief Generates from a set of Fasta files a 2D-datastructure which stores all theoretical masses of all
   * b and y ions from all peptides generated from the Fasta file. The datastructure is build such that on one axis
   * the fragments are sorted by their own mass and the axis by the mass of their precursor/protein.
   * The FI has two options: Bottom-up and Top Down. In later digestion is skiped and the fragments have a direct
   * reference to the mass of the proteins instead of digested peptides.
   */
  class OPENMS_DLLAPI FragmentIndex : public DefaultParamHandler
  {
  public:



    /** @brief Peptide with all important infos needed for the FI-structure
     */
    struct Peptide {

      // We need a constructor in order to emplace back
      Peptide(UInt32 protein_idx, UInt32 modification_idx, std::pair<uint16_t , uint16_t> sequence, float precursor_mz):
          protein_idx(protein_idx),
        modification_idx_(modification_idx),
        sequence_(sequence),
        precursor_mz_(precursor_mz)
        {}

        UInt32 protein_idx;            ///< The index in fasta entries vector
        UInt32 modification_idx_;        ///< The modification index which needed to reconstruct the modification
        std::pair<uint16_t , uint16_t> sequence_; ///< The substring of the protein at position protein_idx
        float precursor_mz_;                  ///< mz of the peptide
    };

    /**
     * @brief Match between a query peak and an entry in the DB
     */
    struct SpectrumMatch
    {
      uint32_t num_matched_{};      ///< Number of peaks-fragment hits
      uint16_t precursor_charge_{};  ///< The precursor_charged used for the performed search
      int16_t isotope_error_{};      /// < The isotope_error used for the performed search
      size_t peptide_idx_{};         ///< The idx this struct belongs to
    };


    /**
     * @brief container for SpectrumMatch. Also keeps count of total number of candidates and total number of matches.
     */
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
    /// DefaultConstructor
    FragmentIndex();

    /// Default destructor
    ~FragmentIndex() override = default;

    // returns true if already built, false otherwise
    bool isBuild() const;

    const std::vector<Peptide>& getPeptides() const;

#ifdef DEBUG_FRAGMENT_INDEX
    /**
     * @brief (only for debugging) Manually add a peptide to the DB. The peptide can be modified
     * IMPORTANT: add peptide after the other peptides are generated and before the FragmentIndex is build
     * @param peptide
     */
    void addSpecialPeptide(AASequence& peptide, Size source_idx);
#endif

    /** @brief Given a set of Fasta files, builds the Fragment Index datastructure (FID). First all fragments are sorted
     * by their own mass. Next they are placed in buckets. The min-fragment mass is stored for each bucket, whereupon
     * the fragments are sorted within the buckets by their originating precursor mass.
     *
     * @param fasta_entries
     */
    void build(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);

    /** @brief Delete fragment index. Sets is_build=false*/
    void clear();


    /** Return index range of all possible Peptides/Proteins, such that a vector can be created fitting that range (safe some memory)
     * @param precursor_mass The mono-charged precursor mass (M+H)
     * @param window Defines the lower and upper bound for the precusor mass. For closed search it only contains the tolerance. In case of open search
     *                  it contains both tolerance and open-search-window
     * @return a pair of indexes defining all possible peptides which the current peak could hit
     */
    std::pair<size_t, size_t> getPeptidesInPrecursorRange(float precursor_mass, const std::pair<float, float>& window);

    /**
     * A match between a single query peak and a database fragment
     */
    struct Hit
    {
      Hit(UInt32 peptide_idx, float fragment_mz) :
        peptide_idx(peptide_idx),
        fragment_mz(fragment_mz)
      {}
      UInt32 peptide_idx; // index in database
      float fragment_mz;
    };

    /**@brief Queries one peak
     * @param peak The queried peak
     * @param peptide_idx_range The range of precursors/peptides the peptide could potentially belongs to
     * @param peak_charge The charge of the peak. Is used to calculate the mass from the mz
     * @return a vector of Hits(matching peptide_idx_range and matching fragment_mz_) containing the idx of the hitted peptide and the mass of the hit
     */
    std::vector<Hit> query(const Peak1D& peak, const std::pair<size_t, size_t>& peptide_idx_range, const uint16_t peak_charge);

    /**
     * @brief: queries one complete experimental spectra against the Database. Loops over all precursor charges
     * Starts at min_precursor_charge and iteratively goes to max_precursor_charge. We query all peaks multiple times with all the
     * different precursor charges and corresponding precursor masses
     * @param spectrum experimental spectrum
     * @param sms[out] The n best Spectrum matches
     */
    void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms);

protected:

  /**@brief One entry in the fragment index
   */
  struct Fragment
  {
      Fragment(UInt32 peptide_idx, float fragment_mz):
          peptide_idx_(peptide_idx),
          fragment_mz_(fragment_mz)
      {}
      UInt32 peptide_idx_; // 32 bit in sage
      float fragment_mz_;
  };

    bool is_build_{false};              ///< true, if the database has been populated with fragments

    void updateMembers_() override;

     /**@brief Generates all peptides from given fasta entries. If Bottom-up is set to false
     * skips digestion. If set to true the Digestion enzyme can be set in the parameters.
     * Additionally introduces fixed and variable modifications for restrictive PSM search.
     *
     * @param fasta_entries
     */
    void generatePeptides(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

    std::vector<Peptide> fi_peptides_;   ///< vector of all (digested) peptides
    std::vector<Fragment> fi_fragments_; ///< vector of all theoretical fragments (b- and y- ions)

    float fragment_min_mz_;  ///< smallest fragment mz
    float fragment_max_mz_;  ///< largest fragment mz    
    size_t bucketsize_;       ///< number of fragments per outer node
    std::vector<float> bucket_min_mz_;  ///< vector of the smalles fragment mz of each bucket
    float precursor_mz_tolerance_;
    bool precursor_mz_tolerance_unit_ppm_{true};
    float fragment_mz_tolerance_;
    bool fragment_mz_tolerance_unit_ppm_{true};    
private:


    /**
     * @brief queries exactly one peak, with a set range of potential peptides, isotope error and precursor charge Hits are transfered into the a PSM
     * Technically an adapter between query(...) and openSearch(...)/searchDifferentPrecursorRanges(...)
     * @param[out] candidates The n best Spectrum matches
     * @param peak The queried peak
     * @param candidates_range The range of precursors/peptides the peptide could potentially belongs to
     * @param isotope_error The applied isotope_error
     * @param precursor_charge The applied precursor charge
     */
    void queryPeak(SpectrumMatchesTopN& candidates,
                   const Peak1D& peak,
                   const std::pair<size_t, size_t>& candidates_range,
                   const int16_t isotope_error,
                   const uint16_t precursor_charge);
    /**
     * @brief If closed search loops over all isotope errors. For each iteration loop over all peaks with queryPeak.
     * @brief If open search applies a precursor-mass window
     * @param spectrum experimental query-spectrum
     * @param precursor_mass The mass of the precursor (mz * charge)
     * @param[out] sms The Top m SpectrumMatches
     * @param charge Applied charge
     */
    void searchDifferentPrecursorRanges(const MSSpectrum& spectrum, float precursor_mass, SpectrumMatchesTopN& sms, uint16_t charge);

    /** @brief places the k-largest elements in the front of the input array. Inside of the k-largest elements and outside the elements are not sorted
     *
     */
    void trimHits(SpectrumMatchesTopN& init_hits) const;

    //since we work with TheoreticalSpectrumGenerator, we must transfer some of those member variables
    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;

    // SpectrumGenerator independend member variables
    std::string digestion_enzyme_;

    size_t missed_cleavages_; ///< number of missed cleavages
    float peptide_min_mass_;
    float peptide_max_mass_;
    size_t peptide_min_length_;
    size_t peptide_max_length_;
  
    StringList modifications_fixed_;    ///< Modification that are one all peptides
    StringList modifications_variable_; ///< Variable Modification -> all possible comibnations are created
    size_t max_variable_mods_per_peptide_;



    // Search Related member variables

    uint16_t min_matched_peaks_;  ///< PSM with less hits are discarded
    int16_t min_isotope_error_;   ///< Minimal possible isotope error
    int16_t max_isotope_error_;   ///< Maximal possible isotope error (both only used for closed search)
    uint16_t min_precursor_charge_; ///< minimal possible precursor charge (usually always 1)
    uint16_t max_precursor_charge_; ///< maximal possible precursor charge
    uint16_t max_fragment_charge_;  ///< The maximal possible charge of the fragments
    uint32_t max_processed_hits_;   ///< The amount of PSM that will be used. the rest is filtered out
    bool open_search;               ///< true if a unrestrictive open search for potential PTM is performed
    float open_precursor_window_lower_; ///< Defines the lower bound of the precursor-mass range
    float open_precursor_window_upper_; ///< Defines the upper bound of the precursor-mass range


  };

}
