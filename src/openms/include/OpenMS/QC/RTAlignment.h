// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  class FeatureMap;
  class PeptideIdentification;
  class TransformationDescription;

  /**
    @brief Take the original retention time before map alignment and use the alignment's trafoXML
           for calculation of the new alignment retention times.

    Sets meta values "rt_raw" and "rt_align" in PeptideIdentifications of the featureMap's PepIDs.
    It does <b>not</b> change the RT of the features.

    **/
  class OPENMS_DLLAPI RTAlignment : public QCBase
  {
  public:
    /// Constructor
    RTAlignment() = default;

    /// Destructor
    virtual ~RTAlignment() = default;

    /**
     @brief Calculates retention time after map alignment
            and sets meta values "rt_raw" and "rt_align" in all PepIDs (on features and all unassigned PepIDs)

     @param fm: FeatureMap to receive the new metavalues
     @param trafo: Transformation information to get needed data from
    **/
    void compute(FeatureMap& fm, const TransformationDescription& trafo) const;

    /**
    @brief Calculates retention time after map alignment
    and sets meta values "rt_raw" and "rt_align" in all PepIDs

    @param ids: PepIDs to receive the new metavalues
    @param trafo: Transformation information to get needed data from
    **/
    void compute(std::vector<PeptideIdentification>& ids, const TransformationDescription& trafo) const;

    /// returns the name of the metric
    const String& getName() const override;

    /// define the required input file: featureXML before map alignment (=POSTFDRFEAT), trafoXML after map alignment (=TRAFOALIGN)
    Status requirements() const override;

  private:
    /// name of the metric
    const String name_ = "RTAlignment";
  };
} // namespace OpenMS
