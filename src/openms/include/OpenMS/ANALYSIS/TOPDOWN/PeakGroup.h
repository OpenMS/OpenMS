//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
    void updateMonomassAndIsotopeIntensities();

    /**
           @brief Update isotope cosine sore and qscore. Mono mass is also updated one last time. SNR, per charge SNR, and avg errors are updated here.
           @param avg precalculated averagine
           @param min_cos the peak groups with cosine score less than this will have Qscore 0.
           @return returns isotope offset after isotope cosine calculation
      */
    int updateIsotopeCosineSNRAvgErrorAndQscore(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos);

    /**
     * @brief given a monoisotopic mass, recruit raw peaks from the raw input spectrum and add to this peakGroup. This is a bit time-consuming and is done for only a small number of selected
     * high-quality peakgroups.
     * @param spec raw spectrum
     * @param tol mass tolerance
     * @param avg precalculated averagine
     * @param mono_mass monoisotopic mass
     * @param excluded_peak_mzs mzs that will be included - only for dummy generation
     * @param charge_offset charge offset from peaks to recruited peaks
     * @param charge_multiple charge multiplication factor for recruited peaks
     * @param mz_off mz offset for recruited peaks
     * @return returns the noisy peaks for this peakgroup - i.e., the raw peaks within the range of this peakGroup that are not matched to any istope of this peakGroup mass.
     */
    std::vector<LogMzPeak> recruitAllPeaksInSpectrum(const MSSpectrum& spec, double tol, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double mono_mass,
                                                     const std::unordered_set<double>& excluded_peak_mzs, int charge_offset = 0, double charge_multiple = 1.0, double mz_off = .0);

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

    /// set representative max_qscore_charge
    void setRepAbsCharge(int max_qscore_charge);

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
     * set peakGroup q value for different TargetDummyType. Q values are stored per TargetDummyType and later used for final q value calculation.
     * @param  target_dummy_type  This target_dummy_type_ specifies if a PeakGroup is a target (0), charge dummy (1), noise dummy (2), or isotope dummy (3)
     */
    void setQvalue(float q, PeakGroup::TargetDummyType target_dummy_type);

    /// set distance between consecutive isotopes
    void setIsotopeDaDistance(double d);

    /// get distance between consecutive isotopes
    double getIsotopeDaDistance() const;

    /// set index of this peak group
    void setIndex(uint i);

    /// update per charge cosine values
    void updatePerChargeCos(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg);

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
    void shrink_to_fit();
    void sort();

  private:
    /// set per abs_charge signal power
    void setChargePowers_(int abs_charge, float signal_pwr, float noise_pwr, float intensity);
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

    void clear_();

    /// log Mz peaks
    std::vector<FLASHDeconvHelperStructs::LogMzPeak> logMzpeaks_;

    /// per charge SNR, isotope cosine, and intensity vectors
    std::vector<float> per_charge_signal_pwr_;
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
    int scan_number_;
    /// is positive or not
    bool is_positive_;
    /// if this peak group has been targeted
    bool is_targeted_ = false;
    /// information on the deconvolved mass
    double monoisotopic_mass_ = -1.0;
    float intensity_; // total intensity
    /// index to specify if this peak_group is a target (0), an isotope dummy (1), a noise (2), or a charge dummy (3)
    PeakGroup::TargetDummyType target_dummy_type_ = target;

    /// distance between consecutive isotopes. Can be different for dummys
    double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;
    /// scoring variables
    int max_qscore_abs_charge_ = -1;
    float isotope_cosine_score_ = 0;
    float charge_score_;
    float qscore_ = .0f;
    float avg_ppm_error_ = 0;
    float avg_da_error_ = 0;
    float snr_ = 0;
    /// q values with different dummy types
    std::map<PeakGroup::TargetDummyType, float> qvalue_;
  };
} // namespace OpenMS