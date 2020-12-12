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
    explicit PeakGroup(int minCharge, int maxCharge);

    /// default destructor
    ~PeakGroup();

    /// copy constructor
    PeakGroup(const PeakGroup &) = default;

    /// move constructor
    PeakGroup(PeakGroup &&other) = default;

    /// clear per charge vectors to save memory
    void clearChargeInfo();

    /// comparison operators
    bool operator<(const PeakGroup &a) const;

    bool operator>(const PeakGroup &a) const;

    bool operator==(const PeakGroup &a) const;

    /// assignment operator
    PeakGroup& operator = (const PeakGroup &t) = default;

    /**
           @brief adjust monoisotopic indices of peaks by offset and discard negative isotpe peaks. Total intensity is also updated
           @param offset isotope index offset
           @param maxIsoIndex max isotopic index
      */
    void updateMassesAndIntensity(int offset = 0,
                                  int maxIsoIndex = 0);
    /// set scan number
    void setScanNumber(int sn);
    /// set per charge SNR
    void setChargeSNR(int charge, float snr);
    /// set per charge isotope cosine
    void setChargeIsotopeCosine(int charge, float cos);
    /// set per charge intensity
    void setChargeIntensity(int charge, float intensity);
    /// set mz range that results in max Qscore
    void setMaxQScoreMzRange(double min, double max);
    /// set min and max charge range
    void setChargeRange(int min, int max);
    /// set isotope cosine score
    void setIsotopeCosine(float cos);
    /// set representative charge
    void setRepCharge(int c);
    /// set Q score
    void setQScore(float q);
    /// set total SNR
    void setSNR(float snr);
    /// set charge score
    void setChargeScore(float s);

    /// get scan number
    int getScanNumber() const;
    /// get monoisotoopic mass
    double getMonoMass() const;
    /// get intensity
    double getIntensity() const;
    /// get per charge SNR
    float getChargeSNR(int charge) const;
    /// get per charge isotope cosine
    float getChargeIsotopeCosine(int charge) const;
    /// get per charge intenstiy
    float getChargeIntensity(int charge) const;
    /// get mz range that results in max Qscore
    std::tuple<double, double> getMzxQScoreMzRange() const;
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
    std::vector<float> perChargeSNR;
    std::vector<float> perChargeCos;
    std::vector<float> perChargeInt;

    /// mz range resulting in maximum Q score
    double maxQScoreMzEnd, maxQScoreMzStart;
    /// charge range
    int maxCharge, minCharge;
    /// scan number
    int scanNumber;

    /// information on the deconvouted mass
    double monoisotopicMass;
    double intensity;// total intensity

    /// scoring variables
    int maxQScoreCharge;
    float isotopeCosineScore;
    float chargeScore;
    float qScore;
    float totalSNR;
  };
}