// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/FAIMSHelper.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
  std::set<double> FAIMSHelper::getCompensationVoltages(const PeakMap& exp)
  {
    std::set<double> CVs;

    if (exp.getSpectra().empty())
    {
      return CVs;
    }

    // Determine if this dataset actually contains FAIMS spectra (scan all; do not rely on the first spectrum)
    const bool has_faims = std::any_of(exp.begin(), exp.end(), [](const MSSpectrum& s)
    {
      return s.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
    });

    if (!has_faims)
    {
      return CVs;
    }

    // Collect compensation voltages from spectra that actually carry FAIMS CVs
    for (const auto& spec : exp)
    {
      if (spec.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
      {
        CVs.insert(spec.getDriftTime());
      }
    }

    // Remove placeholder/sentinel values while warning the user
    if (CVs.erase(IMTypes::DRIFTTIME_NOT_SET) > 0)
    {
      OPENMS_LOG_WARN << "Warning: FAIMS compensation voltage is missing for at least one spectrum!" << std::endl;
    }

    return CVs;
  }

  PeptideIdentificationList FAIMSHelper::filterPeptidesByFAIMSCV(
    const PeptideIdentificationList& peptides,
    double target_cv,
    double cv_tolerance)
  {
    PeptideIdentificationList filtered;
    filtered.reserve(peptides.size()); // Optimistic reservation

    for (const auto& pep : peptides)
    {
      if (pep.metaValueExists(Constants::UserParam::FAIMS_CV))
      {
        double pep_cv = pep.getMetaValue(Constants::UserParam::FAIMS_CV);
        if (std::abs(pep_cv - target_cv) < cv_tolerance)
        {
          filtered.push_back(pep);
        }
      }
      else
      {
        // IDs without FAIMS_CV annotation - include in all groups (backward compatibility)
        filtered.push_back(pep);
      }
    }

    return filtered;
  }

}  //end namespace OpenMS
