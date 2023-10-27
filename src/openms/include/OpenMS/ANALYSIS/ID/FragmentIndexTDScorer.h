//
// Created by trapho on 10/12/23.
//

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTD.h>
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>

#include <vector>
#include <functional>


namespace OpenMS
{
  /**@brief This class is used to get the most promising hits within a FID given a spectrum.
   * For unrestrictive PTM search open-search can be used.
   * Also methods for scoring the best candidates are supported
   *
   */
  class OPENMS_DLLAPI FragmentIndexTDScorer : public DefaultParamHandler
  {
  private:
    FragmentIndexTD db_;       ///< The FragmentIndex(TD) database
    uint16_t min_matched_peaks_;  ///< PSM with less hits are discarded
    int16_t min_isotope_error_;   ///< Minimal possible isotope error
    int16_t max_isotope_error_;   ///< Maximal possible isotope error (both only used for closed search)
    uint16_t min_precursor_charge_; ///< minimal possible precursor charge (usually always 1)
    uint16_t max_precursor_charge_; ///< maximal possible precursor charge
    uint32_t max_processed_hits_;   ///< The amount of PSM that will be used. the rest is filtered out
    bool open_search;               ///< true if a unrestrictive open search for potential PTM is performed
    double open_fragment_window;    // TODO: Should be dependent of prec_window and also adaptable!!! so I guess i could delete this later
    double open_precursor_window;   ///< (Only for open_search) The window in which the precursor-mz could be located, due to PTM



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

    /**For some reason templated functions can not be called in different classes? So this is the work around
     * A function only for TagGeneration
     *
     * @param slice
     * @param size
     * @param access
     */
    void buildKMinHeapforTag(std::vector<TagGenerator::IdxAndIntensity>& slice, uint32_t size);

    void testHeapify();

    /// getter
    const FragmentIndexTD& getDb() const;

    /**
     * @brief Every potential Peptide/Protein has such an struct. Inside the number of peaks-to-Fragment hits are safed
     */
    struct PreHits{
      uint32_t num_matched_; ///< Number of peaks-fragment hits
      size_t peptide_idx_;   ///< The idx this struct belongs to
      uint32_t precursor_charge_;  ///< The precursor_charged used for the performed search
      int16_t isotope_error_;      /// < The isotope_error used for the performed search

      PreHits(){
        num_matched_ = 0;
        peptide_idx_ = 0;
        precursor_charge_ = 0;
        isotope_error_ = 0;
      }
    };

    /**
     * @brief A container for the PreHits
     */
    struct InitHits{
      uint32_t matched_peaks_;
      uint32_t scored_candidates_;
      std::vector<PreHits> hits_;

      InitHits(){
        matched_peaks_ = 0;
        scored_candidates_ = 0;

      }

      InitHits& operator+=(const InitHits& other){
        this-> matched_peaks_ += other.matched_peaks_;
        this->scored_candidates_ += other.scored_candidates_;
        this->hits_.insert(this->hits_.end(), other.hits_.begin(), other.hits_.end());
        return *this;
      }
    };

    /// Constructor
    FragmentIndexTDScorer();

    /// FID setter and builder
    void setDB(FragmentIndexTD& db);
    void buildDB(const std::vector<FASTAFile::FASTAEntry> & fasta_entries);

    /// The "scoring" function
    void simpleScoring(MSSpectrum& spectrum, InitHits& initHits);

    /// Closed Scoring. NO open window for PTM search
    void closedScoring(MSSpectrum& spectrum, double mz, InitHits& initHits, uint32_t charge);

    /** @brief The idea is to search witch an wider precursor tolerance window and actively adjust the fragment window
     *
     */
    void openScoring(MSSpectrum& spectrum, double precursor_mass, InitHits& initHits, uint32_t charge);



  protected:
    void updateMembers_() override;

    /** @brief places the k-largest elements in the front of the input array. Inside of the k-largest elements and outside the elements are not sorted
     * Brazenly stolen from sage.
     */
    void trimHits(InitHits& init_hits);




  };
}
