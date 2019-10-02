// Created by Kyowon Jeong on 9/25/19.
//

#ifndef OPENMS_HOST_MASSFEATURETRACE_H
#define OPENMS_HOST_MASSFEATURETRACE_H

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
  class OPENMS_DLLAPI MassFeatureTrace
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    MassFeatureTrace();

    /// default destructor
    ~MassFeatureTrace();

  protected:

  };
}



#endif //OPENMS_HOST_MASSFEATURETRACE_H
