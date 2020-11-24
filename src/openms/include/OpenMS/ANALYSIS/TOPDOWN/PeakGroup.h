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
  // data structure for peak group. A mass contains multiple peaks of different charges and isotope indices
  class OPENMS_DLLAPI PeakGroup :
      public std::vector<FLASHDeconvHelperStructs::LogMzPeak>
  {
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
  public:
    PeakGroup() = default;

    explicit PeakGroup(int maxChargeRange);

    ~PeakGroup() = default;

    void clearChargeInfo();

    bool operator<(const PeakGroup &a) const;

    bool operator>(const PeakGroup &a) const;

    bool operator==(const PeakGroup &a) const;

    void updateMassesAndIntensity(int offset = 0,
                                  int maxIsoIndex = 0);

    void setScanNumber(int sn);

    void setMonoMass(double m);

    void setIntensity(double i);

    void setChargeSNR(int charge, float snr);

    void setChargeIsotopeCosine(int charge, float cos);

    void setChargeIntensity(int charge, float intensity);

    void setMaxQScoreMzRange(double min, double max);

    void setChargeRange(int min, int max);

    void setIsotopeCosine(float cos);

    void setRepCharge(int c);

    void setQScore(float q);

    void setSNR(float snr);

    int getScanNumber() const;

    double getMonoMass() const;

    double getIntensity() const;

    float getChargeSNR(int charge) const;

    float getChargeIsotopeCosine(int charge) const;

    float getChargeIntensity(int charge) const;

    std::tuple<double, double> getMzxQScoreMzRange() const;

    std::tuple<int, int> getChargeRange() const;

    float getIsotopeCosine() const;

    int getRepCharge() const;

    float getQScore() const;

    float getSNR() const;
    // information on scoring

  private:  // the spectrum from which a peak group is generated
    // the peaks contained in this peak group
    std::vector<float> perChargeSNR;
    std::vector<float> perChargeCos;
    std::vector<float> perChargeInt;

    double maxQScoreMzEnd, maxQScoreMzStart;
    int maxCharge, minCharge;

    // indexing information
    int scanNumber;

    // information on the deconvouted mass
    double monoisotopicMass;
    double intensity;// total intensity

    int maxQScoreCharge;
    float isotopeCosineScore;
    float qScore;
    float totalSNR;
  };
}