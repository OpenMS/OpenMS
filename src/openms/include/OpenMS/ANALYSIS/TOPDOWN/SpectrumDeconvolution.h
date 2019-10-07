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
#include <Eigen/Dense>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
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

    std::vector<PeakGroup> getPeakGroupsFromSpectrum(std::vector<std::vector<Size>> &prevMassBinVector,
                                   std::vector<double> &prevMinBinLogMassVector,
                                                     FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

    static double getChargeFitScore(double *perChargeIntensity, int range);

    static double getIsotopeCosineAndDetermineIsotopeIndex(double mass,
                                                           double *perIsotopeIntensities,
                                                           int perIsotopeIntensitiesSize,
                                                           int &offset,
                                                           FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

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

    Byte **updateMassBins_(boost::dynamic_bitset<> candidateMassBinsForThisSpectrum,
                           float *massIntensities,
                           long &binStart, long &binEnd
    );

    Byte **updateMassBins(double &massBinMinValue,
                          float *massIntensities,
                          float *mzIntensities
    );

    void removeOverlappingPeakGroups();

    static double getCosine(std::vector<double> &a, std::vector<double> &b, int off = 0);

    static double getCosine(double *a,
                            int &aStart,
                            int &aEnd,
                            IsotopeDistribution &b,
                            int &bSize,
                            double &bNorm,
                            int offset);

    void updatePerChargeIsotopeIntensity(
        double *perIsotopeIntensity,
        double *perChargeIntensity,
        PeakGroup &pg);

    static bool checkChargeDistribution(double *perChargeIntensity, int range, int threshold);

    void scoreAndFilterPeakGroups(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

    boost::dynamic_bitset<> getCandidateMassBinsForThisSpectrum(float *massIntensitites, float *mzIntensities);

    void getCandidatePeakGroups(double &mzBinMinValue,
                                double &massBinMinValue,
                                float *sumLogIntensities,
                                Byte **chargeRanges
    );


    void setFilters();
  };
}

#endif //OPENMS_HOST_SPECTRUMDECONVOLUTION_H
