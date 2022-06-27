//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

namespace OpenMS
{
  /**
@brief  Class describing a deconvolved mass.
   A mass contains multiple peaks of different charges and isotope indices.
   PeakGroup is the set of such peaks representing a single monoisotopic mass.
   PeakGroup also contains features that define the quality of it. It is used by QScore calculation
   DeconvolvedSpectrum consists of PeakGroups.
@ingroup Topdown
*/

  class OPENMS_DLLAPI PeakGroup
  {
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
  public:
    /// default constructor
    PeakGroup() = default;

    /**
           @brief Constructor specifying charge range
           @param min_abs_charge min Charge
           @param max_abs_charge max Charge
           @param is_positive whether MS is positive mode
      */
    explicit PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive);

    /// default destructor
    ~PeakGroup();

    /// copy constructor
    PeakGroup(const PeakGroup& ) = default;

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
           @param offset isotope index offset
           @param max_isotope_index max isotopic index
      */
    void updateMassesAndIntensity(const int offset = 0,
                                  const int max_isotope_index = 0);


    double recruitAllPeaksInSpectrum(const MSSpectrum& spec, const double tol, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,  double mono_mass = 0);


    /// using signal and total (signal + noise) power, update SNR value
    void updateSNR();

    /// determine is an mz is a signal of this peakgroup. Input tol is ppm tolerance (e.g., 10.0 for 10ppm tolerance)
    bool isSignalMZ(const double mz, const double tol) const;

    /// set scan number
    void setScanNumber(const int scan_number);

    /// set per abs_charge isotope cosine
    void setChargeIsotopeCosine(const int abs_charge, const float cos);


    /// set mz range that results in max Qscore
    void setMaxQScoreMzRange(const double min, const double max);

    /// set min_abs_charge and max_abs_charge charge range
    void setAbsChargeRange(const int min_abs_charge, const int max_abs_charge);

    /// set isotope cosine score
    void setIsotopeCosine(const float cos);

    /// set representative max_qscore_charge
    void setRepAbsCharge(const int max_qscore_charge);

    /// set Q score - for FLASHIda log file parsing
    void setQScore(const float qscore);

    /// set charge score - for FLASHIda log file parsing
    void setChargeScore(const float charge_score);

    /// set average mass ppm error
    void setAvgPPMError(const float error);

    /// set SNR manually - for FLASHIda log file parsing
    void setSNR(const float snr);
    /// set charge SNR manually - for FLASHIda log file parsing
    void setChargeSNR(const int abs_charge, const float c_snr);

    /// set if it is targeted
    void setTargeted();

    /// get scan number
    int getScanNumber() const;

    /// get monoisotopic mass
    double getMonoMass() const;

    /// get intensity
    double getIntensity() const;

    /// get per abs_charge SNR
    float getChargeSNR(const int abs_charge) const;

    /// get per abs_charge isotope cosine
    float getChargeIsotopeCosine(const int abs_charge) const;

    /// get per abs_charge intenstiy
    float getChargeIntensity(const int abs_charge) const;

    /// get mz range that results in max QScore
    std::tuple<double, double> getMaxQScoreMzRange() const;

    /// get mz range of the charge
    std::tuple<double, double> getMzRange(int abs_charge) const;

    /// get charge range - the actual charge values
    std::tuple<int, int> getAbsChargeRange() const;

    /// get isotopic cosine score
    float getIsotopeCosine() const;

    /// get representative charge
    int getRepAbsCharge() const;

    /// get Q score
    float getQScore() const;

    /// get total SNR
    float getSNR() const;

    /// get charge score
    float getChargeScore() const;

    /// get average mass ppm error;
    float getAvgPPMError() const;

    /// get if it is positive mode
    bool isPositive() const;

    /// get if it is targeted
    bool isTargeted() const;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator begin() const noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator end() const noexcept;

    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator begin() noexcept;
    std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator end() noexcept;

    const FLASHDeconvHelperStructs::LogMzPeak& operator[](const Size i) const;

    void push_back (const FLASHDeconvHelperStructs::LogMzPeak& pg);
    Size size() const noexcept;
    void clear();
    void reserve (Size n);
    bool empty() const;
    void swap (std::vector<FLASHDeconvHelperStructs::LogMzPeak>& x);
    void shrink_to_fit();
    void sort();

  private:

    /// set per abs_charge signal power
    void setChargeSignalPower_(const int abs_charge, const double pwr);
    /// set per abs_charge intensity
    void setChargeIntensity_(const int abs_charge, const float intensity);
    /// set per abs_charge total peak power
    void setChargePower_(const int abs_charge, const double pwr);

    void setChargeFitScore_();

    /// log Mz peaks
    std::vector<FLASHDeconvHelperStructs::LogMzPeak> logMzpeaks_;
    /// per charge SNR, isotope cosine, and intensity vectors
    std::vector<float> per_charge_signal_pwr_;
    std::vector<float> per_charge_pwr_;
    std::vector<float> per_charge_cos_;
    std::vector<float> per_charge_int_;
    std::vector<float> per_charge_snr_;

    /// mz range resulting in maximum Q score
    double max_qscore_mz_end_, max_qscore_mz_start_;
    /// charge range
    int min_abs_charge_ = 0, max_abs_charge_ = -1;
    /// scan number
    int scan_number_;
    /// is positive or not
    bool is_positive_;
    /// if this peak group has been targeted
    bool is_targeted_ = false;
    /// information on the deconvolved mass
    double monoisotopic_mass_ = -1.0;
    double intensity_;// total intensity


    /// scoring variables
    int max_qscore_abs_charge_;
    float isotope_cosine_score_;
    float charge_score_;
    float qscore_;
    float avg_ppm_error_;
    float snr_;
  };
}
