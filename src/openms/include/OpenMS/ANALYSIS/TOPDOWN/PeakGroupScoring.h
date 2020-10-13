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

  #include <omp.h>
#endif
namespace OpenMS
{
  class OPENMS_DLLAPI PeakGroupScoring
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    //typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    PeakGroupScoring(std::vector<PeakGroup> &peakGroups, Parameter &param);

    /// default destructor
    ~PeakGroupScoring();

    static double getChargeFitScore(double *perChargeIntensity, int range);

    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           double *perIsotopeIntensities,
                                                           int perIsotopeIntensitiesSize,
                                                           int &offset,
                                                           FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

    std::vector<PeakGroup> &scoreAndFilterPeakGroups(unsigned int &msLevel,
                                                     FLASHDeconvHelperStructs::PrecalculatedAveragine &avg);

  protected:
    std::vector<PeakGroup> &peakGroups;
    Parameter &param;


    void removeOverlappingPeakGroups(double tol);

    void removeHarmonicPeakGroups(double tol);

    std::vector<int> updatePerChargeIsotopeIntensity(
        //        double **intensityGrid,
        //        double **intensityGrid2,
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg);

    void filterPeakGroupsByIsotopeCosine(int currentMaxMassCount);

    void filterPeakGroupsByQScore(int currentMaxMassCount);

    void filterPeakGroupsByIntensity(int currentMaxMassCount);

    double getPredictionScore(PeakGroup &pg, int charge); //

    //static double getAvgMassPpmError(PeakGroup &pg);


    static bool checkChargeDistribution(double *perChargeIntensity, int range, int threshold);


    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    static double getCosine(const double *a, double *b, Size size);

    static double getCorrelation(const double *a,
                                 int &aStart,
                                 int &aEnd,
                                 IsotopeDistribution &b,
                                 int &bSize,
                                 int offset);

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
