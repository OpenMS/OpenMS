// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTagGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <functional>
#include <vector>


namespace OpenMS
{
  /**@brief This class is used to get the most promising hits within a FID given a spectrum.
   * For unrestrictive PTM search open-search can be used.
   * Also methods for scoring the best candidates are supported
   *
   */
  class OPENMS_DLLAPI FragmentIndexScorer : public DefaultParamHandler
  {
  private:
    std::string fragmentation_;
    FragmentIndex db_;            ///< The FragmentIndex database
    uint16_t min_matched_peaks_;  ///< PSM with less hits are discarded
    int16_t min_isotope_error_;   ///< Minimal possible isotope error
    int16_t max_isotope_error_;   ///< Maximal possible isotope error (both only used for closed search)
    uint16_t min_precursor_charge_; ///< minimal possible precursor charge (usually always 1)
    uint16_t max_precursor_charge_; ///< maximal possible precursor charge
    uint16_t max_fragment_charge_;
    uint32_t max_processed_hits_;   ///< The amount of PSM that will be used. the rest is filtered out
    bool open_search;               ///< true if a unrestrictive open search for potential PTM is performed
    float open_fragment_window;    // TODO: Should be dependent of prec_window and also adaptable!!! so I guess i could delete this later
    float open_precursor_window;   ///< (Only for open_search) The window in which the precursor-mz could be located, due to PTM
    std::string fragmentation_method;

    /**
     * @brief simple minHeapify for entry idx inside the first size elements of an vector
     * @tparam T The class the input vector is filled with
     * @tparam A A member variable within class T, which is used for performing heapify
     * @param slice The input vector
     * @param idx The element inside the vector on which heapify is performed with
     * @param size The size of the heap within the vector
     * @param access The function to access A within T
     */
    template<class T, class A>
    void minHeapify(std::vector<T>& slice, size_t idx, uint32_t size, A (*access) (T));

    /**
     * @brief Using minHeapify build a complete min-Heap
     * @tparam T
     * @tparam A
     * @param slice
     * @param size
     * @param access
     */
    template<class T, class A>
    void buildKMinHeap(std::vector<T>& slice, uint32_t size, A (*access)(T));

  public:

    const FragmentIndex& getDB() const;

    /**For some reason templated functions can not be called in different classes? So this is the work around
     * A function only for TagGeneration
     *
     * @param slice
     * @param size
     * @param access
     */
    void buildKMinHeapforTag(std::vector<TagGenerator::IdxAndIntensity>& slice, uint32_t size);

    void testHeapify();

    /**
     * @brief Every potential Peptide/Protein has such an struct. Inside the number of peaks-to-Fragment hits are safed
     */
    struct SpectrumMatch
    {
      uint32_t num_matched_{};       ///< Number of peaks-fragment hits
      uint32_t precursor_charge_{};  ///< The precursor_charged used for the performed search
      size_t peptide_idx_{};         ///< The idx this struct belongs to
      int16_t isotope_error_{};      /// < The isotope_error used for the performed search
    };

    /**
     * @brief all SpectrumMatchesTopN of a spectrum, the number of matched peaks and the number of scored candidates
     */
    struct SpectrumMatchesTopN
    {
      std::vector<SpectrumMatch> hits_;     ///< The preliminary candidates
      uint32_t matched_peaks_{};      ///< The number of matched peaks TODO: statistic needed?
      uint32_t scored_candidates_{};  ///< The number of scored candidates

      SpectrumMatchesTopN() = default;

      SpectrumMatchesTopN& operator+=(const SpectrumMatchesTopN& other)
      {
        this->matched_peaks_ += other.matched_peaks_;
        this->scored_candidates_ += other.scored_candidates_;
        this->hits_.insert(this->hits_.end(), other.hits_.begin(), other.hits_.end());
        return *this;
      }
    };

    /// Constructor
    FragmentIndexScorer();

    void buildDB(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);  // Only builds standard

    void extractHits(SpectrumMatchesTopN& candidates ,
                     const std::vector<FragmentIndex::Hit>& hits,
                     uint32_t charge,
                     int16_t isotope_error,
                     std::pair<size_t , size_t > peptide_range);

    /// The "scoring" function
    void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms);

    void queryPeak(SpectrumMatchesTopN& candidates, const Peak1D& peak, std::pair<size_t, size_t> candidates_range, int16_t isotope_error, uint16_t precursor_charge);
    
    /// Closed Scoring. NO open window for PTM search
    void closedSearch(const MSSpectrum& spectrum, float mz, SpectrumMatchesTopN& sms, uint16_t charge);

    /** @brief The idea is to search witch an wider precursor tolerance window and actively adjust the fragment window
     *
     */
    void openSearch(const MSSpectrum& spectrum, float precursor_mass, SpectrumMatchesTopN& sms, uint16_t charge);

/*
     * @brief Scoring for the MultiDim FragmentIndex!
     * @param spectrum
     * @param SpectrumMatchesTopN
    void multiDimScoring(const MSSpectrum& spectrum, SpectrumMatchesTopN& SpectrumMatchesTopN);
*/

  protected:
    void updateMembers_() override;

    /** @brief places the k-largest elements in the front of the input array. Inside of the k-largest elements and outside the elements are not sorted
     * Brazenly stolen from sage.
     */
    void trimHits(SpectrumMatchesTopN& init_hits);

  };
}
