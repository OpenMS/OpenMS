// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{
  /**
@brief  Class describing a deconvolved mass.
   A mass contains multiple (LogMz) peaks of different charges and isotope indices.
   PeakGroup is the set of such peaks representing a single monoisotopic mass.
   PeakGroup also contains features that define the quality of it. It is used by Qscore calculation.
   DeconvolvedSpectrum consists of PeakGroups.
@ingroup Topdown
*/

  class OPENMS_DLLAPI PeakGroup
  {
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;

  public:
    /// target dummy type of PeakGroup. This specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3).
    enum TargetDummyType
    {
      target = 0,
      charge_dummy,
      noise_dummy,
      isotope_dummy
    };


    /// default constructor
    PeakGroup() = default;

    /**
           @brief Constructor specifying charge range
           @param min_abs_charge min Charge
           @param max_abs_charge max Charge
           @param is_positive whether MS is positive mode
      */
    explicit PeakGroup(int min_abs_charge, int max_abs_charge, bool is_positive);

    /// default destructor
    ~PeakGroup() = default;

    /// copy constructor
    PeakGroup(const PeakGroup&) = default;

    /// move constructor
    PeakGroup(PeakGroup&& other) = default;

    /// comparison operators
    bool operator<(const PeakGroup& a) const;

    bool operator>(const PeakGroup& a) const;

    bool operator==(const PeakGroup& a) const;

    /// assignment operator
    PeakGroup& operator=(const PeakGroup& t) = default;

    /**
           @brief add monoisotopic indices of peaks by offset and discard negative isotope peaks. Total intensity is also updated
      */
    void updateMonoMassAndIsotopeIntensities();

    /**
           @brief Update Qscore. Cosine and SNRs are also updated.
           @param noisy_peaks noisy peaks to calculate Qscore
           @param avg precalculated averagine
           @param min_cos the peak groups with cosine score less than this will have Qscore 0.
           @param allowed_iso_error this set the allowed isotope error in dummy mass generation.
           @return returns isotope offset after isotope cosine calculation
      */
    int updateQscore(std::vector<LogMzPeak>& noisy_peaks, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos, int allowed_iso_error = 1);

    /**
     * @brief given a monoisotopic mass, recruit raw peaks from the raw input spectrum and add to this peakGroup. This is a bit time-consuming and is done for only a small number of selected
     * high-quality peakgroups.
     * @param spec raw spectrum
     * @param tol mass tolerance
     * @param avg precalculated averagine
     * @param mono_mass monoisotopic mass
     * @param excluded_peak_mzs mzs that will be included - only for dummy generation
     * @return returns the noisy peaks for this peakgroup - i.e., the raw peaks within the range of this peakGroup that are not matched to any istope of this peakGroup mass.
     */
    std::vector<LogMzPeak> recruitAllPeaksInSpectrum(const MSSpectrum& spec, double tol, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double mono_mass,
                                                     const std::unordered_set<double>& excluded_peak_mzs);

    /// determine is an mz is a signal of this peakgroup. Input tol is ppm tolerance (e.g., 10.0 for 10ppm tolerance). Assume logMzPeaks are sorted.
    bool isSignalMZ(double mz, double tol) const;

    /// set scan number
    void setScanNumber(int scan_number);

    /// set per abs_charge isotope cosine
    void setChargeIsotopeCosine(int abs_charge, float cos);

    /// set min_abs_charge and max_abs_charge charge range
    void setAbsChargeRange(int min_abs_charge, int max_abs_charge);

    /// set isotope cosine score
    void setIsotopeCosine(float cos);

    /// set representative max_snr_abs_charge
    void setRepAbsCharge(int max_snr_abs_charge);

    /// set monoisotopic mass
    void setMonoisotopicMass(double mono_mass);

    /// set Q score - for FLASHIda log file parsing
    void Qscore(float qscore);

    /// set charge score - for FLASHIda log file parsing
    void setChargeScore(float charge_score);

    /// set average mass ppm error
    void setAvgPPMError(float error);

    /// set SNR manually - for FLASHIda log file parsing
    void setSNR(float snr);

    /// set charge SNR manually - for FLASHIda log file parsing
    void setChargeSNR(int abs_charge, float c_snr);

    /// set if it is targeted
    void setTargeted();

    /// get scan number
    int getScanNumber() const;

    /// get monoisotopic mass
    double getMonoMass() const;

    /// get intensity
    float getIntensity() const;

    /// get per abs_charge SNR
    float getChargeSNR(int abs_charge) const;

    /// get per abs_charge isotope cosine
    float getChargeIsotopeCosine(int abs_charge) const;

    /// get per abs_charge intenstiy
    float getChargeIntensity(int abs_charge) const;

    /// get mz range that results in max Qscore
    std::tuple<double, double> getRepMzRange() const;

    /// get mz range of the charge
    std::tuple<double, double> getMzRange(int abs_charge) const;

    /// get charge range - the actual charge values
    std::tuple<int, int> getAbsChargeRange() const;

    /// get per isotope intensities
    const std::vector<float>& getIsotopeIntensities() const;

    /// get isotopic cosine score
    float getIsotopeCosine() const;

    /// get representative charge
    int getRepAbsCharge() const;

    /// get Q score
    float getQscore() const;

    /// get total SNR
    float getSNR() const;

    /// get charge score
    float getChargeScore() const;

    /// get average mass ppm error;
    float getAvgPPMError() const;

    /// get average mass ppm error;
    float getAvgDaError() const;

    /// get if it is positive mode
    bool isPositive() const;

    /// get if it is targeted
    bool isTargeted() const;

    /// get the target dummy type of this
    PeakGroup::TargetDummyType getTargetDummyType() const;

    /// for this PeakGroup, specify the target dummy type.
    void setTargetDummyType(PeakGroup::TargetDummyType index);

    /**
     * Get q values for different target_dummy_type. For charge, noise, isotope dummy types, q values corresponding to the type will be returned. For target (default), the final q value is calculated
     * by summing the q values of all dummy types and returned.
     * @param  target_dummy_type  This target_dummy_type_ specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
     * @return Q value of the peakGroup
     */
    float getQvalue(PeakGroup::TargetDummyType target_dummy_type = PeakGroup::TargetDummyType::target) const;

    /**
     * Set peakGroup q-value for different TargetDummyType. Q values are stored per TargetDummyType and later used for final q value calculation.
     * @param q The q-value
     * @param target_dummy_type  This target_dummy_type_ specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
     */
    void setQvalue(float q, PeakGroup::TargetDummyType target_dummy_type);

    /// set distance between consecutive isotopes
    void setIsotopeDaDistance(double d);

    /// get distance between consecutive isotopes
    double getIsotopeDaDistance() const;

    /// get minimum neagative isotope index
    int getMinNegativeIsotopeIndex() const;
    /// set index of this peak group
    void setIndex(uint i);

    /// get index of this peak group
    uint getIndex() const;

    /// iterators for the signal LogMz peaks in this PeakGroup
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator begin() const noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator end() const noexcept;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator begin() noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator end() noexcept;

    const FLASHDeconvHelperStructs::LogMzPeak& operator[](Size i) const;

    /// iterators for the noisy LogMz peaks in this PeakGroup
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator getNoisePeakBegin() const noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator getNoisePeakEnd() const noexcept;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator getNoisePeakBegin() noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator getNoisePeakEnd() noexcept;

    /// vector operators for the LogMzPeaks in this PeakGroup
    void push_back(const FLASHDeconvHelperStructs::LogMzPeak& pg);
    Size size() const noexcept;

    void reserve(Size n);
    bool empty() const;
    void swap(std::vector<FLASHDeconvHelperStructs::LogMzPeak>& x);
    void sort();

  private:
    /// update chargefit score and also update per charge intensities here.
    void updateChargeFitScoreAndChargeIntensities_();
    /// update avg ppm error
    void updateAvgPPMError_();
    /// update avg Da error
    void updateAvgDaError_();
    /// get ppm error of a logMzPeak
    float getAbsPPMError_(const LogMzPeak& p) const;
    /// get Da error of a logMzPeak from the closest isotope
    float getAbsDaError_(LogMzPeak& p) const;
    /// using signal and total (signal + noise) power, update SNR value
    void updateSNR_();
    /// clear peaks
    void clear_();
    /// update per charge intensities, noise power, and squared intensities. used for SNR estimation
    void updatePerChargeInformation_(const std::vector<LogMzPeak>& noisy_peaks);
    /// update the charge range using the calculated per charge information
    void updateChargeRange_(std::vector<LogMzPeak>& noisy_peaks);
    /// update per charge cosine values
    void updatePerChargeCos_(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg);

    /**
     * calculate noisy peak power. The goal of this function is to group noisy peaks that are possibly from the same molecule and sum their intensities before calculate power
     * @param noisy_peaks noisy peaks to calculate power
     * @param signal_peaks signal peaks - they may make a part of noisy isotopes
     * @return calculated noise power
     */
    float getNoisePeakPower_(const std::vector<LogMzPeak>& noisy_peaks, const std::vector<LogMzPeak>& signal_peaks) const;

    /// log Mz peaks
    std::vector<FLASHDeconvHelperStructs::LogMzPeak> logMzpeaks_;
    /// negative isotope index peaks
    std::vector<FLASHDeconvHelperStructs::LogMzPeak> negative_iso_peaks_;
    /// per charge SNR, isotope cosine, and intensity vectors
    std::vector<float> per_charge_sum_signal_squared_;
    std::vector<float> per_charge_noise_pwr_;
    std::vector<float> per_charge_cos_;
    std::vector<float> per_charge_int_;
    std::vector<float> per_charge_snr_;
    /// per isotope intensity.
    std::vector<float> per_isotope_int_;
    /// charge range
    int min_abs_charge_ = 0, max_abs_charge_ = -1;
    /// peak group index
    uint index_ = 0;
    /// scan number
    int scan_number_ = 0;
    /// is positive or not
    bool is_positive_ = false;
    /// if this peak group has been targeted
    bool is_targeted_ = false;
    /// information on the deconvolved mass
    double monoisotopic_mass_ = -1.0;
    float intensity_ = 0; // total intensity
    /// index to specify if this peak_group is a target (0), an isotope dummy (1), a noise (2), or a charge dummy (3)
    PeakGroup::TargetDummyType target_dummy_type_ = target;
    /// up to which negative isotope index should be considered. By considereing negative istoopes, one can reduce isotope index error.
    int min_negative_isotope_index_ = -1;
    /// distance between consecutive isotopes. Can be different for dummys
    double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;
    /// scoring variables
    int max_snr_abs_charge_ = -1;
    float isotope_cosine_score_ = 0;
    float charge_score_ = 0;
    float qscore_ = .0f;
    float avg_ppm_error_ = 0;
    float avg_da_error_ = 0;
    float snr_ = 0;
    /// q values with different dummy types
    std::map<PeakGroup::TargetDummyType, float> qvalue_;
  };
} // namespace OpenMS