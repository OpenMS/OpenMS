// Created by Kyowon Jeong on 9/25/19.
//

#ifndef OPENMS_HOST_MASSFEATURETRACE_H
#define OPENMS_HOST_MASSFEATURETRACE_H

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
//#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>

#include "boost/dynamic_bitset.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <Eigen/Dense>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

namespace OpenMS
{
  class OPENMS_DLLAPI MassFeatureTrace
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
   // typedef FLASHDeconvHelperStructs::DeconvolutedSpectrum DeconvolutedSpectrum;

    /// default constructor
    MassFeatureTrace(Parameter &param, Param &mtd_param, PrecalculatedAveragine &averagines);

    /// default destructor
    ~MassFeatureTrace();

    void addDeconvolutedSpectrum(DeconvolutedSpectrum &deconvolutedSpectrum);

    void findFeatures(int &featureCntr, int &featureIndex, std::fstream &fsf, std::fstream &fsp);

    static void writeHeader(std::fstream &fs);
    static void writePromexHeader(std::fstream &fs);

  protected:
    Parameter &param;
    Param &mtd_param;
    PrecalculatedAveragine &averagines;

    std::unordered_map<double, std::unordered_map<double, PeakGroup>> peakGroupMap; // rt , mono mass, peakgroup
  };
}



#endif //OPENMS_HOST_MASSFEATURETRACE_H
