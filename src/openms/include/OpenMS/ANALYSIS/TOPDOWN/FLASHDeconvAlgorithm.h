// --------------------------------------------------------------------------
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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
//#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>

#include "boost/dynamic_bitset.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

namespace OpenMS
{
  /** NEED to be modified
  @brief @ref

  @htmlinclude

  @ingroup Topdown
*/
  /** check OPENMS c++ guide **/
  class OPENMS_DLLAPI FLASHDeconvAlgorithm
  {
public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    FLASHDeconvAlgorithm();

    /// default destructor
    ~FLASHDeconvAlgorithm();

    /// copy constructor
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm &);

    /// assignment operator
    FLASHDeconvAlgorithm &operator=(const FLASHDeconvAlgorithm &fd);

    static double getChargeFitScore(double *perChargeIntensity, int range);
    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           double *perIsotopeIntensities,
                                                           int perIsotopeIntensitiesSize,
                                                           PrecalcularedAveragine &averagines,
                                                           int &offset);
    static std::vector<PeakGroup> Deconvolution(MSExperiment &map,
                                                Parameter &param,
                                                PrecalcularedAveragine &averagines,
                                                int &specCntr,
                                                int &qspecCntr,
                                                int &massCntr);

  protected:
    int getNominalMass(double &m);

    static double getBinValue(Size bin, double minV, double binWidth);

    static Size getBinNumber(double v, double minV, double binWidth);

  private:

    static void printProgress(float progress);

    static std::vector<LogMzPeak> getLogMzPeaks(MSSpectrum &spec, const Parameter &param);

    static std::vector<PeakGroup> getPeakGroupsFromSpectrum(std::vector<LogMzPeak> &logMzPeaks,
                                                     double *filter,
                                                     double **harmonicFilter,
                                                     std::vector<std::vector<Size>> &prevMassBinVector,
                                                     std::vector<double> &prevMinBinLogMassVector,
                                                     PrecalcularedAveragine &averagines,
                                                     const Parameter &param,
                                                     int &specCntr);

    static boost::dynamic_bitset<> getUnionMassBin(boost::dynamic_bitset<> &massBins,
                                            double &massBinMinValue,
                                            std::vector<std::vector<Size>> &prevMassBinVector,
                                            std::vector<double> &prevMassBinMinValue,
                                            const Parameter &param);

    static std::vector<PeakGroup> getPeakGroupsWithMassBins(boost::dynamic_bitset<> &unionedMassBins,
                                                     std::vector<LogMzPeak> &logMzPeaks,
                                                     double &mzBinMinValue,
                                                     double &massBinMinValue,
                                                     float *sumLogIntensities,
                                                     long *binOffsets,
                                                     Byte **chargeRanges,
                                                     const Parameter &param);

    static boost::dynamic_bitset<> getMzBins(std::vector<LogMzPeak> &logMzPeaks,
                                      double &mzBinMinValue,
                                      Size &binNumber,
                                      double binWidth,
                                      float *intensities);

    static Byte **getMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                       double &massBinMinValue,
                       float *sumLogIntensities,
                       long *binOffsets,
                       long **hBinOffsets,
                       boost::dynamic_bitset<> &unionMassBins,
                       float *intensities,
                       const Parameter &param, double &minMass, double &maxMass);

    void
    printMasses(boost::dynamic_bitset<> &massBins, double &massBinMinValue, Byte *continuousChargePeakPairCount,
                const Parameter &param);

    static void getInitialMassBins(boost::dynamic_bitset<> &massBins,
                            boost::dynamic_bitset<> &mzBins,
                            boost::dynamic_bitset<> &isQualified,
                            float *signal,
                            long **hBinOffsets,
                            long *binOffsets,
                            float *intensities,
                            const Parameter &param);

    static Byte **getFinalMassBins(boost::dynamic_bitset<> &massBins, boost::dynamic_bitset<> &mzBins,
                            boost::dynamic_bitset<> &isQualified,
                            boost::dynamic_bitset<> &unionMassBins,
                            float *sumLogIntensities,
        // double &massBinMinValue,
                            long *binOffsets,
                            const Parameter &param,
                            long &binStart, long &binEnd);

    static std::vector<FLASHDeconvAlgorithm::PeakGroup> scoreAndFilterPeakGroups(std::vector<PeakGroup> &peakGroups,
                                                    PrecalcularedAveragine &averagines,
                                                    const Parameter &param);

    static void removeOverlappingPeakGroups(std::vector<PeakGroup> &pgs, double tol);


    static void updatePerChargeIsotopeIntensity(//int *perIsotopeMinCharge, int *perIsotopeMaxCharge,
        //int *perChargeMinIsotope, int *perChargeMaxIsotope,
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg,
        const Parameter &param);


    static bool checkSpanDistribution(int *mins, int *maxs, int range, int threshold);



    static bool checkChargeDistribution(double *perChargeIntensity,
                                 int range,
                                 int threshold);

    static double getCosine(double *a,
                     int &aStart,
                     int &aEnd,
                     IsotopeDistribution &b,
                     int &bSize,
                     double &bNorm,
                     int offset = 1);

    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    void filterPeakGroupsByIntensity(std::vector<PeakGroup> &peakGroups,
                                     std::vector<double> &intensities,
                                     const Parameter &param);

  };
}// namespace OpenMS
