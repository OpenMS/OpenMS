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
  struct OPENMS_DLLAPI PeakGroup
  {
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;

    // the spectrum from which a peak group is generated
    MSSpectrum *spec;
    // the peaks contained in this peak group
    std::vector<LogMzPeak> peaks;
    // information for identity
    int scanNumber, specIndex;

    // information on the deconvouted mass
    double monoisotopicMass;
    double avgMass;
    double intensity;// total intensity
    Size massBinIndex;
    int massIndex;

    // information on scoring
    float isotopeCosineScore;
    float chargeCosineScore;
    int maxQScoreCharge;
    float qScore;
    float totalSNR;

    // other information appeared on the output file
    int maxCharge, minCharge;
    double maxQScoreMzEnd, maxQScoreMzStart;

    // necessary temp information for scoring.
    std::unordered_map<int, std::vector<float>> perChargeInfo; // charge -> SNR, ICos, SumInt

    ~PeakGroup();

    void push_back(LogMzPeak &p);

    void reserve(Size n);

    bool empty();

    void clearChargeInfo();

    bool operator<(const PeakGroup &a) const;

    bool operator>(const PeakGroup &a) const;

    bool operator==(const PeakGroup &a) const;

    void updateMassesAndIntensity(FLASHDeconvHelperStructs::PrecalculatedAveragine &averagines,
                                  int offset = 0,
                                  int maxIsoIndex = 0);
  };
}