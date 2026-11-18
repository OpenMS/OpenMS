// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <memory>

namespace OpenMS
{
  /**
  @brief FLASHDeconv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
  From MSSpectrum, this class outputs DeconvolvedSpectrum.
  Deconvolution takes three steps:
   i) decharging and select candidate masses - speed up via binning
   ii) collecting isotopes from the candidate masses and deisotope - peak groups are defined here
   iii) scoring and filter out low scoring masses (i.e., peak groups)
  @ingroup Topdown
*/

  class OPENMS_DLLAPI SpectralDeconvolution : public DefaultParamHandler
  {
  public:
    typedef FLASHHelperClasses::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHHelperClasses::LogMzPeak LogMzPeak;

    /// default constructor
    SpectralDeconvolution();

    /// copy constructor
    SpectralDeconvolution(const SpectralDeconvolution&) = default;

    /// move constructor
    SpectralDeconvolution(SpectralDeconvolution&& other) = default;

    /// assignment operator
    SpectralDeconvolution& operator=(const SpectralDeconvolution& fd) = default;

    /// move assignment operator
    SpectralDeconvolution& operator=(SpectralDeconvolution&& fd) = default;

    /// destructor
    ~SpectralDeconvolution() = default;

    /**
      @brief main deconvolution function that generates the deconvolved target and dummy spectrum based on the original spectrum.
      @param spec the original spectrum
      @param scan_number scan number from input spectrum.
      @param precursor_peak_group precursor peak group
    */
    void performSpectrumDeconvolution(const MSSpectrum& spec, int scan_number, const PeakGroup& precursor_peak_group);

    /// return deconvolved spectrum
    DeconvolvedSpectrum& getDeconvolvedSpectrum();

    /// get calculated averagine. Call after calculateAveragine is called.
    const PrecalculatedAveragine& getAveragine();

    /// set calculated averagine
    void setAveragine(const PrecalculatedAveragine& avg);

    /** @brief set targeted or excluded masses for targeted deconvolution. Masses are targeted or excluded in all ms levels.
     *  @param masses target masses to set
        @param exclude if set, masses are excluded.
     */
    void setTargetMasses(const std::vector<double>& masses, bool exclude = false);

    /** @brief precalculate averagine (for predefined mass bins) to speed up averagine generation
        @param use_RNA_averagine if set, averagine for RNA (nucleotides) is calculated
       */
    void calculateAveragine(bool use_RNA_averagine);

    /// when estimating tolerance, max_mass_dalton_tolerance_ should be large
    void setToleranceEstimation()
    {
      max_mass_dalton_tolerance_ = 100;
    }

    /// convert double to nominal mass
    static int getNominalMass(double mass);

    /** calculate cosine between two vectors a and b with additional parameters for fast calculation
     * @param a vector a
     * @param a_start non zero start index of a
     * @param a_end non zero end index of a (exclusive)
     * @param b vector b
     * @param offset element index offset between a and b
     * @param min_iso_len minimum isotope size. If isotope size is less than this, return 0
     */
    static float getCosine(const std::vector<float>& a, int a_start, int a_end, const IsotopeDistribution& b, int offset, int min_iso_len);


    /** @brief Examine intensity distribution over isotope indices. Also determines the most plausible isotope index or, monoisotopic mono_mass
        @param mono_mass monoisotopic mass
        @param per_isotope_intensities vector of intensities associated with each isotope - aggregated through charges
        @param offset output offset between input monoisotopic mono_mass and determined monoisotopic mono_mass
        @param avg precalculated averagine
        @param iso_int_shift isotope shift in per_isotope_intensities.
        @param window_width isotope offset value range. If -1, set automatically.
        @param excluded_masses masses not considered in the calculation.
        @return calculated cosine similarity score
     */
    static float getIsotopeCosineAndIsoOffset(double mono_mass, const std::vector<float>& per_isotope_intensities, int& offset,
                                              const PrecalculatedAveragine& avg, int iso_int_shift,
                                              int window_width, const std::vector<double>& excluded_masses);

    /**
     *  set target dummy type for the SpectralDeconvolution run. All masses from the target SpectralDeconvolution run will have the target_decoy_type_.
     * @param target_decoy_type  This target_decoy_type_ specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
     * @param target_dspec_for_decoy_calcualtion target masses from normal deconvolution
     */
    void setTargetDecoyType(PeakGroup::TargetDecoyType target_decoy_type, const DeconvolvedSpectrum& target_dspec_for_decoy_calcualtion);

    /// minimum isotopologue count in a peak group
    const static int min_iso_size = 2;

  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconvAlgorithm attaches one immutable averagine to its per-thread spectrum workers (setSharedAveragine_)
    friend class FLASHDeconvAlgorithm;

    /**
      @brief Use @p avg, an immutable averagine shared by several instances (the spectrum workers of FLASHDeconvAlgorithm), instead of the
      own one. calculateAveragine and setAveragine return to the own averagine.
      @param avg the shared averagine (not null; not modified while shared)
      @throws Exception::InvalidParameter if @p avg is null
    */
    void setSharedAveragine_(std::shared_ptr<const PrecalculatedAveragine> avg);

    /// frees the pooled per-spectrum arrays; the decoy instances of an FDR run call it after each spectrum, so a worker keeps one set
    /// of pools instead of three
    void releasePools();

    /// the averagine in use: the shared one if set, otherwise the own one
    const PrecalculatedAveragine& activeAveragine_() const noexcept
    {
      return shared_avg_ ? *shared_avg_ : avg_;
    }

    /// FLASHDeconv parameters

    /// allowed isotope error in deconvolved mass to calculate qvalue
    int allowed_iso_error_ = -1;

    /// min charge and max charge subject to analysis, set by users
    int min_abs_charge_, max_abs_charge_;
    /// is positive mode
    bool is_positive_;
    /// mass ranges of deconvolution, set by users
    double min_mass_, max_mass_;
    /// current_max_charge_: controlled by precursor charge for MSn n>1; otherwise just max_abs_charge_
    int current_max_charge_;
    /// max mass is controlled by precursor mass for MSn n>1; otherwise just max_mass
    double current_max_mass_;
    /// max mass is max_mass for MS1 and 50 for MS2
    double current_min_mass_;
    /// isotope distance for noise decoy. This distance has been determined so i) it is close to the mass of 13C - 12C yet is not close to any correct isotope distances of different charges (from 1 to 10).
    double noise_iso_delta_ =  .9444;
    /// minimum number of peaks supporting a mass
    const static int min_support_peak_count_ = 2;
    /// tolerance in ppm for each MS level
    DoubleList tolerance_;
    /// bin multiplication factor (log mz * bin_mul_factors_ = bin number) - for fast convolution, binning is used
    DoubleList bin_mul_factors_;
    /// Isotope cosine threshold for each MS level
    DoubleList min_isotope_cosine_;
    /// snr threshold for each MS level
    DoubleList min_snr_;
    /// minimum Qscore of a deconvolved mass
    double min_qscore_ = .2;

    /// the peak group vector from normal run. This is used when dummy masses are generated.
    const DeconvolvedSpectrum* target_dspec_for_decoy_calculation_;

    /// PeakGroup::TargetDecoyType values
    PeakGroup::TargetDecoyType target_decoy_type_ = PeakGroup::TargetDecoyType::target;

    /// precalculated averagine distributions for fast averagine generation
    FLASHHelperClasses::PrecalculatedAveragine avg_;
    /// an immutable averagine shared with other instances (setSharedAveragine_); used instead of avg_ when set
    std::shared_ptr<const PrecalculatedAveragine> shared_avg_;

    /// mass bins that are targeted for FLASHIda global targeting mode
    boost::dynamic_bitset<> target_mass_bins_;
    std::vector<double> target_mono_masses_;

    /// mass bins that are excluded for FLASHIda global targeting mode
    std::vector<double> excluded_masses_;

    /// mass bins that are previously deconvolved and excluded for decoy mass generation
    boost::dynamic_bitset<> excluded_mass_bins_for_decoy_runs_;
    std::vector<double> excluded_peak_masses_for_decoy_runs_;
    std::vector<double> excluded_masses_for_decoy_runs_;

    /// Stores log mz peaks
    std::vector<LogMzPeak> log_mz_peaks_;
    /// selected_peak_groups_ stores the deconvolved mass peak groups
    DeconvolvedSpectrum deconvolved_spectrum_;
    /// binned_log_masses_ stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> binned_log_masses_;
    /// log m/z of every peak of the original spectrum of deconvolved_spectrum_, computed once per spectrum (FLASHHelperClasses::getLogMz)
    std::vector<double> original_log_mzs_;
    /// nearest-peak lookup of the original spectrum of deconvolved_spectrum_, built once per spectrum (PeakGroup::recruitAllPeaksInSpectrum)
    PeakGroup::NearestPeakIndex original_nearest_peaks_;

    /// binned_log_mz_peaks_ stores the binned log mz peaks
    boost::dynamic_bitset<> binned_log_mz_peaks_;

    /// This stores the "universal pattern"
    std::vector<double> universal_pattern_;
    /// This stores the patterns for harmonic reduction
    Matrix<double> harmonic_pattern_matrix_;

    /// This stores the "universal pattern" in binned dimension
    std::vector<int> binned_universal_pattern_;
    /// This stores the patterns for harmonic reduction in binned dimension
    Matrix<int> binned_harmonic_patterns;

    /// minimum mass and mz values representing the first bin of massBin and mzBin, respectively: to save memory space
    double mass_bin_min_value_;
    double mz_bin_min_value_;

    /// current ms Level
    uint ms_level_;

    /// isotope dalton distance (set per spectrum; initialised for setTargetMasses before the first deconvolution)
    double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;

    /// for precursor targetting.
    int target_precursor_charge_ = 0;
    double target_precursor_mz_ = 0;
    /// mass of the isolated precursor ion (precursor_mass): gives the precursor charge of each MSn spectrum when precursor_charge is 0
    double target_precursor_mass_ = 0;

    /// this is additional mass tolerance in Da to get more high signal-to-ratio peaks in this candidate peakgroup finding
    double max_mass_dalton_tolerance_ = .16;

    /// a mass bin selected by filterMassBins_ and the range of the charge indices (charge - 1) it was selected for
    struct MassBinChargeRange_
    {
      Size bin;
      int min_charge_index;
      int max_charge_index;
    };

    /// the state per mass bin of updateMassBins_, kept from spectrum to spectrum so that it is not allocated and filled for every spectrum
    /// (see updateCandidateMassBins_); a copy starts empty
    struct MassBinScratch_
    {
      /// the rest of the state of a mass bin in the consecutive-charge chain search, valid only where prev_charges holds a tag of the current spectrum
      struct ChainState
      {
        float prev_intensity;
        float prev_residual;
        unsigned short support_peak_count;
        short chain_offset;
      };
      /// charge index of the last peak per mass bin, as a tag: charge_tag_base + charge index in the current spectrum; lower values: none
      std::vector<unsigned short> prev_charges;
      std::vector<ChainState> chain_states;
      /// the tag of charge index 0 in the next spectrum; 0: prev_charges must be reset first
      int charge_tag_base = 0;
      /// intensity per mass bin; all zero between spectra (when mass_intensities_clean)
      std::vector<float> mass_intensities;
      /// the mass bins of mass_intensities written in the current spectrum (to zero them again)
      std::vector<Size> written_mass_bins;
      bool mass_intensities_clean = true;

      MassBinScratch_() = default;
      MassBinScratch_(const MassBinScratch_&) {}
      MassBinScratch_(MassBinScratch_&&) = default;
      MassBinScratch_& operator=(const MassBinScratch_&) { return *this; }
      MassBinScratch_& operator=(MassBinScratch_&&) = default;
    };
    MassBinScratch_ mass_bin_scratch_;

    /// Pools for the per-spectrum arrays (sized by the m/z bins, several MB, or by the peaks; the mass-bin state is mass_bin_scratch_) and
    /// the candidate peak group vector:
    /// an instance reuses their capacity from spectrum to spectrum instead of allocating (and faulting in) fresh pages every time. Only
    /// the capacity carries over; every use refills them first, so the values are those of freshly constructed vectors.
    std::vector<float> mz_intensities_pool_;
    std::vector<double> exact_mz_positions_pool_;
    std::vector<Size> mz_bin_index_reverse_pool_;
    std::vector<Size> peak_bin_numbers_pool_;
    std::vector<PeakGroup> candidate_pool_;

    /// whether getCandidatePeakGroups_ scores each candidate as soon as it is built and keeps only those that pass (set per spectrum by
    /// generatePeakGroupsFromSpectrum_: inside a parallel region, where the scoring loop runs serially anyway)
    bool fused_scoring_ = false;
    /// detection limit of the current spectrum (the intensity of its weakest positive peak), set by generatePeakGroupsFromSpectrum_
    float detection_limit_ = 0;

    /** @brief static function that converts bin to value
        @param bin bin number
        @param min_value minimum value (corresponding to bin number = 0)
        @param bin_mul_factor bin multiplication factor: bin_value = (min_value + bin_number/ bin_mul_factors_)
        @return value corresponding to bin
     */
    static double getBinValue_(Size bin, double min_value, double bin_mul_factor);

    /** @brief static function that converts value to bin
        @param value value
        @param min_value minimum value (corresponding to bin number = 0)
        @param bin_mul_factor bin multiplication factor: bin_number = (bin_value * bin_mul_factors_ - min_value)
        @return bin corresponding to value
     */
    static Size getBinNumber_(double value, double min_value, double bin_mul_factor);

    /// generate log mz peaks from the input spectrum
    void updateLogMzPeaks_();

    /** @brief generate mz bins and intensity per mz bin from log mz peaks
        @param bin_number number of mz bins
        @param binned_log_mz_peak_intensities intensity per mz bin
     */
    void binLogMzPeaks_(Size bin_number, std::vector<float>& binned_log_mz_peak_intensities);

    /// get mass value for input mass bin: used for debugging
    double getMassFromMassBin_(Size mass_bin, double bin_mul_factor) const;

    /// get mz value for input mz bin: used for debugging
    double getMzFromMzBin_(Size mass_bin, double bin_mul_factor) const;

    /// Generate peak groups from the input spectrum
    void generatePeakGroupsFromSpectrum_();

    /// filter out overlapping masses
    static void removeOverlappingPeakGroups_(DeconvolvedSpectrum& dspec, double tol, PeakGroup::TargetDecoyType target_decoy_type = PeakGroup::TargetDecoyType::target);

    /** @brief Update binned_log_masses_. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses. This function does not perform deisotoping
        @param mz_intensities per mz bin intensity
        @return the selected mass bins in increasing order, each with its charge range
     */
    std::vector<MassBinChargeRange_> updateMassBins_(const std::vector<float>& mz_intensities);

    /** @brief Subfunction of updateMassBins_.
        @param mass_intensities per mass bin intensity
        @param chain_redirects (m/z bin x current max charge + charge index, mass bin) for the peaks updateCandidateMassBins_ added to the mass bin
               of their consecutive-charge chain instead of their own (rounded) mass bin, sorted here
        @return the selected mass bins in increasing order, each with its charge range
     */
    std::vector<MassBinChargeRange_> filterMassBins_(const std::vector<float>& mass_intensities, std::vector<std::pair<Size, long>>& chain_redirects);

    /** @brief Subfunction of updateMassBins_. It select candidate masses and update binned_log_masses_ using the universal pattern, eliminate possible harmonic masses
        The mass bin intensities go to mass_bin_scratch_.mass_intensities (the written bins to mass_bin_scratch_.written_mass_bins).
        @param mz_intensities mz bin intensities
        @param chain_redirects output: (m/z bin x current max charge + charge index, mass bin) for the peaks added to the mass bin of their
               consecutive-charge chain instead of their own (rounded) mass bin
     */
    void updateCandidateMassBins_(const std::vector<float>& mz_intensities, std::vector<std::pair<Size, long>>& chain_redirects);

    /** @brief For the selected mass bins, select the peaks from the original spectrum. Also isotopic peaks are clustered in this function.
        @param mass_bin_charge_ranges the selected mass bins in increasing order, each with its charge range (filterMassBins_)
     */
    void getCandidatePeakGroups_(const std::vector<MassBinChargeRange_>& mass_bin_charge_ranges);

    /// Make the universal pattern.
    void setFilters_();

    /// function for peak group scoring and filtering
    void scoreAndFilterPeakGroups_();

    /// isotope/charge scoring of one candidate peak group (the first filter of scoreAndFilterPeakGroups_); true if it is kept
    bool scorePeakGroup_(PeakGroup& peak_group);

    /// filter out charge error masses
    void removeChargeErrorPeakGroups_(DeconvolvedSpectrum& dspec, const PeakGroup::TargetDecoyType& target_decoy_type) const;

    /// filter out excluded masses
    void removeExcludedMasses_(DeconvolvedSpectrum& dspec, std::vector<double> excluded_masses, double tol) const;

    void setTargetPrecursorCharge_();

    bool isPeakGroupInExcludedMassForDecoyRuns_(const PeakGroup& peak_group, double tol, int offset = 0) const;
  };
} // namespace OpenMS