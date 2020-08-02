//
// Created by Kyowon Jeong on 9/25/19.
//

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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{

  class OPENMS_DLLAPI SpectrumDeconvolution
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    //typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    SpectrumDeconvolution(MSSpectrum &spec, Parameter &param);

    /// default destructor
    ~SpectrumDeconvolution();

    bool empty();

    std::vector<PeakGroup>& getPeakGroupsFromSpectrum(std::vector<std::vector<Size>> &prevMassBinVector,
                                   std::vector<double> &prevMinBinLogMassVector,
                                                     FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                                     unsigned int msLevel);

    std::vector<LogMzPeak> logMzPeaks;

  protected:
    MSSpectrum &spec;
    Parameter &param;

    boost::dynamic_bitset<> massBinsForThisSpectrum;

    boost::dynamic_bitset<> massBins;
    boost::dynamic_bitset<> mzBins;

    std::vector<PeakGroup> peakGroups;

    long *binOffsets{};
    long **hBinOffsets{};

    double *filter{};
    double **harmonicFilter{};
    //PrecalculatedAveragine averagines;

    static double getBinValue(Size bin, double minV, double binWidth);

    static Size getBinNumber(double v, double minV, double binWidth);

    static void pringMasses(boost::dynamic_bitset<> &massBins, double &minMass, double binWidth);


    void updateLogMzPeaks(double chargeMass);

    void updateMzBins(double &mzBinMinValue, Size &binNumber, double binWidth,
                      float *mzBinIntensities
    );

    void unionPrevMassBins(double &massBinMinValue,
                           std::vector<std::vector<Size>> &prevMassBinVector,
                           std::vector<double> &prevMassBinMinValue,
                           UInt msLevel);

    Byte **updateMassBins_(boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                           float *massIntensities,
                           long &binStart, long &binEnd,
                           unsigned int &msLevel
    );

    Byte **updateMassBins(double &massBinMinValue,
                          double &mzBinMinValue,
                          float *massIntensities,
                          float *mzIntensities,
                          unsigned int &msLevel
    );

    boost::dynamic_bitset<> getCandidateMassBinsForThisSpectrum(float *massIntensitites, float *mzIntensities, double& mzMinValue, unsigned int &msLevel);

    void getCandidatePeakGroups(double &mzBinMinValue,
                                double &massBinMinValue,
                                float *sumLogIntensities,
                                Byte **chargeRanges,
                                FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                unsigned int& msLevel
    );

    void setFilters();
  };
}
