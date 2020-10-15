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

#include <Eigen/Dense>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

#ifdef _OPENMP

  #include <omp.h> // for test..
#endif

namespace OpenMS
{
  class OPENMS_DLLAPI PeakGroupScoring
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    PeakGroupScoring(std::vector<PeakGroup> &peakGroups, Parameter &param);

    /// default destructor
    ~PeakGroupScoring();

    // examine intensity distribution over charges
    static double getChargeFitScore(double *perChargeIntensity, int range);

    // examine intensity distribution over iostope indices. Also determines the most plausible isotope index or, monoisotopic mass
    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           double *perIsotopeIntensities,
                                                           int perIsotopeIntensitiesSize,
                                                           int &offset,
                                                           FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

    //the main function of this class
    std::vector<PeakGroup> &scoreAndFilterPeakGroups(unsigned int &msLevel,
                                                     FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

  protected:
    std::vector<PeakGroup> &peakGroups;
    Parameter &param;

    //filter out overlapping masses
    void removeOverlappingPeakGroups(double tol);

    //filter out possible harmonics
    void removeHarmonicPeakGroups(double tol);

    //From peaks distributions over charge and isotope are calculated
    std::vector<int> updatePerChargeIsotopeIntensity(
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg);

    //Filter out masses with low isotope cosine scores
    void filterPeakGroupsByIsotopeCosine(int currentMaxMassCount);

    //Filter out masses with low QScores
    void filterPeakGroupsByQScore(int currentMaxMassCount);

    //Filter out masses with low intensities
    void filterPeakGroupsByIntensity(int currentMaxMassCount);

    //For MS1, check intensity ratio between charges.
    static bool checkChargeDistribution(double *perChargeIntensity, int range, int threshold);

    //cosine function
    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    //cosine function
    static double getCosine(const double *a, double *b, Size size);

    //correlation function
    static double getCorrelation(const double *a,
                                 int &aStart,
                                 int &aEnd,
                                 IsotopeDistribution &b,
                                 int &bSize,
                                 int offset);

    //cosine function for fast calculatoin
    static double getCosine(const double *a,
                            int &aStart,
                            int &aEnd,
                            IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);

  public:

  };


}
