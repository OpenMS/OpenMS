// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/QCBase.h>

/**
 * @brief Detected Compounds as a Metabolomics QC metric
 *
 * Simple class to return a summary of detected compounds
 * from a featureXML file.
 *
 */

namespace OpenMS
{
  class OPENMS_DLLAPI FeatureSummary : public QCBase
  {
  public:
    /// Constructor
    FeatureSummary() = default;

    /// Destructor
    virtual ~FeatureSummary() = default;

    // stores feature summary values calculated by compute function
    struct OPENMS_DLLAPI Result {
      UInt feature_count = 0;
      float rt_shift_mean = 0;

      bool operator==(const Result& rhs) const;
    };

    /**
   @brief computes a summary of a featureXML file

   @param feature_map FeatureMap
   @return result object with summary values:
           number of detected compounds (detected_compounds),
           retention time shift mean (rt_shift_mean)
   **/
    Result compute(const FeatureMap& feature_map);

    const String& getName() const override;

    QCBase::Status requirements() const override;

  private:
    const String name_ = "Summary of features from featureXML file";
  };
} // namespace OpenMS
