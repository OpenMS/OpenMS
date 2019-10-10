//
// Created by Kyowon Jeong on 9/25/19.
//

#ifndef OPENMS_HOST_SPECTRUMDECONVOLUTION_H
#define OPENMS_HOST_SPECTRUMDECONVOLUTION_H


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

namespace OpenMS
{

  class OPENMS_DLLAPI SpectrumDeconvolution
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    SpectrumDeconvolution(MSSpectrum &spec, Parameter &param);

    /// default destructor
    ~SpectrumDeconvolution();

    std::vector<PeakGroup>& getPeakGroupsFromSpectrum(std::vector<std::vector<Size>> &prevMassBinVector,
                                   std::vector<double> &prevMinBinLogMassVector,
                                                     FLASHDeconvHelperStructs::PrecalcularedAveragine &avg,
                                                     int msLevel);


  protected:
    MSSpectrum spec;
    Parameter param;
    std::vector<LogMzPeak> logMzPeaks;
    boost::dynamic_bitset<> massBinsForThisSpectrum;

    boost::dynamic_bitset<> massBins;
    boost::dynamic_bitset<> mzBins;

    std::vector<PeakGroup> peakGroups;

    long *binOffsets;
    long **hBinOffsets;

    double *filter;
    double **harmonicFilter;
    //PrecalcularedAveragine averagines;

    static double getBinValue(Size bin, double minV, double binWidth);

    static Size getBinNumber(double v, double minV, double binWidth);

    void updateLogMzPeaks();

    void updateMzBins(double &mzBinMinValue, Size &binNumber, double binWidth,
                      float *mzBinIntensities
    );

    void unionPrevMassBins(double &massBinMinValue,
                           std::vector<std::vector<Size>> &prevMassBinVector,
                           std::vector<double> &prevMassBinMinValue);

    Byte **updateMassBins_(boost::dynamic_bitset<> &candidateMassBinsForThisSpectrum,
                           float *massIntensities,
                           long &binStart, long &binEnd,
                           int &msLevel
    );

    Byte **updateMassBins(double &massBinMinValue,
                          float *massIntensities,
                          float *mzIntensities,
                          int &msLevel
    );

    boost::dynamic_bitset<> getCandidateMassBinsForThisSpectrum(float *massIntensitites, float *mzIntensities, int &msLevel);

    void getCandidatePeakGroups(double &mzBinMinValue,
                                double &massBinMinValue,
                                float *sumLogIntensities,
                                Byte **chargeRanges,
                                FLASHDeconvHelperStructs::PrecalcularedAveragine &avg
    );

    void setFilters();
  };
}

#endif //OPENMS_HOST_SPECTRUMDECONVOLUTION_H
