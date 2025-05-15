// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/FAIMSHelper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>

namespace OpenMS
{
  std::set<double> FAIMSHelper::getCompensationVoltages(const PeakMap& exp)
  {
    std::set<double> CVs;

    // is this FAIMS data?
    if ((exp.getSpectra().empty()) ||
        (exp.getSpectra()[0].getDriftTimeUnit() != DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE))
    {
      return CVs;
    }
  
    for (auto it = exp.begin(); it != exp.end(); ++it)
    {
      CVs.insert(it->getDriftTime());
    }

    if (CVs.find(IMTypes::DRIFTTIME_NOT_SET) != CVs.end())
    {
      OPENMS_LOG_WARN << "Warning: FAIMS compensation voltage is missing for at least one spectrum!" << std::endl;
    }
  
    return CVs;
  }

}  //end namespace OpenMS
