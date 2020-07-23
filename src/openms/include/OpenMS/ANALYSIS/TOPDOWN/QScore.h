//
// Created by Kyowon Jeong on 7/20/20.
//


#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CONCEPT/Types.h>
//#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include  <OpenMS/METADATA/Precursor.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

namespace OpenMS
{
  class PeakGroup;

  class OPENMS_DLLAPI QScore
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    static double getQScore(PeakGroup *pg, double intensity, int charge);
    //static double getQScore(PeakGroup *pg, Precursor peak);
  };
}
