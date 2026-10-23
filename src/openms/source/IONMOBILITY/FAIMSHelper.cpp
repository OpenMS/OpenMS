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
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
  namespace
  {
    /// Whether @p item belongs to the FAIMS CV group @p target_cv: its FAIMS_CV matches, or it has none (such items
    /// go to every group, for backward compatibility)
    bool inFAIMSCVGroup(const MetaInfoInterface& item, double target_cv, double cv_tolerance)
    {
      return !item.metaValueExists(Constants::UserParam::FAIMS_CV)
             || std::abs(static_cast<double>(item.getMetaValue(Constants::UserParam::FAIMS_CV)) - target_cv) < cv_tolerance;
    }
  }

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
      if (inFAIMSCVGroup(pep, target_cv, cv_tolerance)) filtered.push_back(pep);
    }
    return filtered;
  }

  FeatureMap FAIMSHelper::filterFeaturesByFAIMSCV(
    const FeatureMap& features,
    double target_cv,
    double cv_tolerance)
  {
    FeatureMap filtered;
    static_cast<DocumentIdentifier&>(filtered) = features;
    for (const Feature& f : features)
    {
      if (inFAIMSCVGroup(f, target_cv, cv_tolerance)) filtered.push_back(f);
    }
    return filtered;
  }

}  //end namespace OpenMS
