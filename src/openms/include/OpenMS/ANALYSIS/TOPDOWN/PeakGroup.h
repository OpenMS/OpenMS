//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
@brief  A mass contains multiple peaks of different charges and isotope indices. PeakGroup is the set of such peaks representing a single monoisotopic mass
@ingroup Topdown
*/

  class OPENMS_DLLAPI PeakGroup :
      public std::vector<FLASHDeconvHelperStructs::LogMzPeak>
  {
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
  public:
    /// default constructor
    PeakGroup() = default;

    /**
           @brief Constructor specifying charge range
           @param minCharge min Charge
           @param maxCharge max Charge
      */
    explicit PeakGroup(const int min_charge, const int max_charge);

    /// default destructor
    ~PeakGroup();

    /// copy constructor
    PeakGroup(const PeakGroup& ) = default;

    /// move constructor
    PeakGroup(PeakGroup&& other) = default;

    /// clear per charge vectors to save memory
    void clearChargeInfo();

    /// comparison operators
    bool operator<(const PeakGroup& a) const;

    bool operator>(const PeakGroup& a) const;

    bool operator==(const PeakGroup& a) const;

    /// assignment operator
    PeakGroup& operator = (const PeakGroup& t) = default;

    /**
           @brief adjust monoisotopic indices of peaks by offset and discard negative isotpe peaks. Total intensity is also updated
           @param offset isotope index offset
           @param maxIsoIndex max isotopic index
      */
    void updateMassesAndIntensity(const int offset = 0,
                                  const int max_isotope_index = 0);
    /// set scan number
    void setScanNumber(const int scan_number);
    /// set per charge SNR
    void setChargeSNR(const int charge, const float snr);
    /// set per charge isotope cosine
    void setChargeIsotopeCosine(const int charge, const float cos);
    /// set per charge intensity
    void setChargeIntensity(const int charge, const float intensity);
    /// set mz range that results in max Qscore
    void setMaxQScoreMzRange(const double min, const double max);
    /// set min and max charge range
    void setChargeRange(const int min, const int max);
    /// set isotope cosine score
    void setIsotopeCosine(const float cos);
    /// set representative charge
    void setRepCharge(const int charge);
    /// set Q score
    void setQScore(const float qscore);
    /// set total SNR
    void setSNR(const float snr);
    /// set charge score
    void setChargeScore(const float charge_score);

    /// get scan number
    int getScanNumber() const;
    /// get monoisotoopic mass
    double getMonoMass() const;
    /// get intensity
    double getIntensity() const;
    /// get per charge SNR
    float getChargeSNR(const int charge) const;
    /// get per charge isotope cosine
    float getChargeIsotopeCosine(const int charge) const;
    /// get per charge intenstiy
    float getChargeIntensity(const int charge) const;
    /// get mz range that results in max Qscore
    std::tuple<double, double> getMaxQScoreMzRange() const;
    /// get charge range
    std::tuple<int, int> getChargeRange() const;
    /// get isotopic cosine score
    float getIsotopeCosine() const;
    /// get representative chrage
    int getRepCharge() const;
    /// get Q score
    float getQScore() const;
    /// get total SNR
    float getSNR() const;
    /// get charge score
    float getChargeScore() const;

  private:
    /// per charge SNR, isotope cosine, and intensity vectors
    std::vector<float> per_charge_snr;
    std::vector<float> per_charge_cos;
    std::vector<float> per_charge_int;

    /// mz range resulting in maximum Q score
    double max_qscore_mz_end, max_qscore_mz_start;
    /// charge range
    int max_charge, min_charge;
    /// scan number
    int scan_number;

    /// information on the deconvouted mass
    double monoisotopic_mass;
    double intensity;// total intensity

    /// scoring variables
    int max_qscore_charge;
    float isotope_cosine_score;
    float charge_score;
    float qscore;
    float total_snr;
  };
}