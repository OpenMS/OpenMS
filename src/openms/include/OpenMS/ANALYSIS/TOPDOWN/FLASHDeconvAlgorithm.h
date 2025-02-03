// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <boost/dynamic_bitset.hpp>
#include <iostream>

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

  class OPENMS_DLLAPI FLASHDeconvAlgorithm : public DefaultParamHandler
  {
  public:
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    FLASHDeconvAlgorithm();

    /// copy constructor
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm&) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& other) = default;

    /// assignment operator
    FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& fd) = default;

    /// move assignment operator
    FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& fd) = default;

    /// destructor
    ~FLASHDeconvAlgorithm() = default;

    /**
      @brief main deconvolution function that generates the deconvolved target and dummy spectrum based on the original spectrum.
      @param spec the original spectrum
      @param survey_scans the survey scans to assign precursor mass to the deconvolved spectrum.
      @param scan_number scan number from input spectrum.
      @param precursor_map_for_FLASHIda deconvolved precursor information from FLASHIda
    */
    void performSpectrumDeconvolution(const MSSpectrum& spec, const std::vector<DeconvolvedSpectrum>& survey_scans, int scan_number,
                                      const std::map<int, std::vector<std::vector<float>>>& precursor_map_for_FLASHIda);

    /// return deconvolved spectrum
    DeconvolvedSpectrum& getDeconvolvedSpectrum();

    /// get calculated averagine. Call after calculateAveragine is called.
    const PrecalculatedAveragine& getAveragine();

    /// set calculated averagine
    void setAveragine(const PrecalculatedAveragine& avg);

    /** @brief set targeted or excluded masses for targeted deconvolution. Masses are targeted or excluded in all ms levels.
        
        @param masses Masses to be targeted or excluded (@p exclude=true)
        @param exclude if set, masses are excluded.
     */
    void setTargetMasses(const std::vector<double>& masses, bool exclude = false);

    /** @brief precalculate averagine (for predefined mass bins) to speed up averagine generation
        @param use_RNA_averagine if set, averagine for RNA (nucleotides) is calculated
     */
    void calculateAveragine(bool use_RNA_averagine);

    /// convert double to nominal mass
    static int getNominalMass(double mass);

    /** calculate cosine between two vectors a and b with additional parameters for fast calculation
     * @param a vector a
     * @param a_start non zero start index of a
     * @param a_end non zero end index of a (exclusive)
     * @param b vector b
     * @param b_size size of b
     * @param offset element index offset between a and b
     * @param min_iso_size minimum isotope size. If isotope size is less than this, return 0
     */
    static float getCosine(const std::vector<float>& a, int a_start, int a_end, const IsotopeDistribution& b, int b_size, int offset, int min_iso_size);


    /** @brief Examine intensity distribution over isotope indices. Also determines the most plausible isotope index or, monoisotopic mono_mass
        @param mono_mass monoisotopic mass
        @param per_isotope_intensities vector of intensities associated with each isotope - aggregated through charges
        @param offset output offset between input monoisotopic mono_mass and determined monoisotopic mono_mass
        @param avg precalculated averagine
        @param iso_int_shift isotope shift in per_isotope_intensities.
        @param window_width isotope offset value range. If -1, set automatically.
        @param allowed_iso_error_for_second_best_cos allowed isotope error to calculate the second best cos. If target_dummy_type is not PeakGroup::TargetDummyType::target, the second best cosine and
       its corresponding offset will be output
        @param target_dummy_type  This target_dummy_type specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
        @return calculated cosine similar score
     */
    static float getIsotopeCosineAndDetermineIsotopeIndex(double mono_mass, const std::vector<float>& per_isotope_intensities, int& offset, const PrecalculatedAveragine& avg, int iso_int_shift = 0,
                                                          int window_width = -1, int allowed_iso_error_for_second_best_cos = 0,
                                                          PeakGroup::TargetDummyType target_dummy_type = PeakGroup::TargetDummyType::target);

    /**
     * add m/zs in input DeconvolvedSpectrum into exclusion list. The exclusion list is used to generate noise dummy masses.
     * @param dspec input DeconvolvedSpectrum
     * @param excluded_mzs mz exclusion list to be updated
     */
    static void addMZsToExcludsionList(const DeconvolvedSpectrum& dspec, std::unordered_set<double>& excluded_mzs);

    /**
     *  set target dummy type for the FLASHDeconvAlgorithm run. All masses from the target FLASHDeconvAlgorithm run will have the target_dummy_type_.
     * @param target_dummy_type  This target_dummy_type_ specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
     * @param target_dspec_for_dummy_calcualtion target masses from normal deconvolution
     */
    void setTargetDummyType(PeakGroup::TargetDummyType target_dummy_type, DeconvolvedSpectrum& target_dspec_for_dummy_calcualtion);

  protected:
    void updateMembers_() override;

  private:
    /// FLASHDeconv parameters

    /// minimum isotopologue count in a peak group
    const static int min_iso_size_ = 2;

    /// allowed isotope error in deconvolved mass to calculate qvalue
    int allowed_iso_error_ = 1;

    /// range of RT subject to analysis (in seconds)
    double min_rt_, max_rt_;
    /// range of mz subject to analysis
    double min_mz_, max_mz_;
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
    /// peak intensity threshold subject to analysis
    double intensity_threshold_;
    /// minimum number of peaks supporting a mass minus one
    const static int min_support_peak_count_ = 2;
    /// tolerance in ppm for each MS level
    DoubleList tolerance_;
    /// bin multiplication factor (log mz * bin_mul_factors_ = bin number) - for fast convolution, binning is used
    DoubleList bin_mul_factors_;
    /// cosine threshold between observed and theoretical isotope patterns for each MS level
    DoubleList min_isotope_cosine_;
    /// the deconvolved spectrum from normal run. This is used when dummy masses are generated.
    DeconvolvedSpectrum* target_dspec_for_dummy_calcualtion_;

    /// PeakGroup::TargetDummyType values
    PeakGroup::TargetDummyType target_dummy_type_ = PeakGroup::TargetDummyType::target;

    /// precalculated averagine distributions for fast averagine generation
    FLASHDeconvHelperStructs::PrecalculatedAveragine avg_;

    /// mass bins that are targeted for FLASHIda global targeting mode
    boost::dynamic_bitset<> target_mass_bins_;
    std::vector<double> target_mono_masses_;

    /// mass bins that are excluded for FLASHIda global targeting mode
    std::vector<double> excluded_masses_;

    /// mass bins that are previsouly deconvolved and excluded for dummy mass generation
    boost::dynamic_bitset<> previously_deconved_mass_bins_for_dummy_;
    std::vector<double> previously_deconved_mono_masses_for_dummy_;
    std::unordered_set<double> excluded_peak_mzs_;

    /// Stores log mz peaks
    std::vector<LogMzPeak> log_mz_peaks_;
    /// selected_peak_groups_ stores the deconvolved mass peak groups
    DeconvolvedSpectrum deconvolved_spectrum_;
    /// mass_bins_ stores the selected bins for this spectrum + overlapped spectrum (previous a few spectra).
    boost::dynamic_bitset<> mass_bins_;
    /// mz_bins_ stores the binned log mz peaks
    boost::dynamic_bitset<> mz_bins_;

    /// This stores the "universal pattern"
    std::vector<double> filter_;
    /// This stores the patterns for harmonic reduction
    Matrix<double> harmonic_filter_matrix_;

    /// isotope dalton distance
    double iso_da_distance_;

    /// This stores the "universal pattern" in binned dimension
    std::vector<int> bin_offsets_;
    /// This stores the patterns for harmonic reduction in binned dimension
    Matrix<int> harmonic_bin_offset_matrix_;

    /// minimum mass and mz values representing the first bin of massBin and mzBin, respectively: to save memory space
    double mass_bin_min_value_;
    double mz_bin_min_value_;

    /// current ms Level
    uint ms_level_;

    /// default precursor isolation window size.
    double isolation_window_size_;

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
        @param mz_bin_intensities intensity per mz bin
     */
    void updateMzBins_(Size bin_number, std::vector<float>& mz_bin_intensities);


    /// get mass value for input mass bin
    double getMassFromMassBin_(Size mass_bin, double bin_mul_factor) const;

    /// get mz value for input mz bin
    double getMzFromMzBin_(Size mass_bin, double bin_mul_factor) const;

    /// Generate peak groups from the input spectrum
    void generatePeakGroupsFromSpectrum_();

    /** @brief Update mass_bins_. It select candidate mass bins using the universal pattern, eliminate possible harmonic masses. This function does not perform deisotoping
        @param mz_intensities per mz bin intensity
        @return a matrix containing charge ranges for all found masses
     */
    Matrix<int> updateMassBins_(const std::vector<float>& mz_intensities);

    /** @brief Subfunction of updateMassBins_.
        @param mass_intensities per mass bin intensity
        @return a matrix containing charge ranges for all found masses
     */
    Matrix<int> filterMassBins_(const std::vector<float>& mass_intensities);

    /** @brief Subfunction of updateMassBins_. It select candidate masses and update mass_bins_ using the universal pattern, eliminate possible harmonic masses
        @param mass_intensities mass bin intensities which are updated in this function
        @param mz_intensities mz bin intensities
     */
    void updateCandidateMassBins_(std::vector<float>& mass_intensities, const std::vector<float>& mz_intensities);

    /** @brief For selected masses in mass_bins_, select the peaks from the original spectrum. Also isotopic peaks are clustered in this function.
        @param per_mass_abs_charge_ranges charge range per mass
     */
    void getCandidatePeakGroups_(const Matrix<int>& per_mass_abs_charge_ranges);

    /// Make the universal pattern.
    void setFilters_();

    /// function for peak group scoring and filtering
    void scoreAndFilterPeakGroups_();

    /// filter out charge error masses
    void removeChargeErrorPeakGroups_(DeconvolvedSpectrum& dspec) const;

    /// filter out overlapping masses
    void removeOverlappingPeakGroups_(DeconvolvedSpectrum& dspec, double tol);

    /**
    @brief register the precursor peak as well as the precursor peak group (or mass) if possible for MSn (n>1) spectrum.
    Given a precursor peak (found in the original MS n-1 Spectrum) the masses containing the precursor peak are searched.
    If multiple masses are detected, the one with the best Qscore is selected. For the selected mass, its corresponding peak group (along with precursor peak) is registered.
    If no such mass exists, only the precursor peak is registered.
    @param survey_scans the candidate precursor spectra - the user may allow search of previous N survey scans.
    @param precursor_map_for_real_time_acquisition this contains the deconvolved mass information from FLASHIda runs.
    */
    bool registerPrecursor_(const std::vector<DeconvolvedSpectrum>& survey_scans, const std::map<int, std::vector<std::vector<float>>>& precursor_map_for_real_time_acquisition);
  };
} // namespace OpenMS
