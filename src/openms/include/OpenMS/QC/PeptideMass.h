// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

namespace OpenMS
{
  class FeatureMap;

  /**
    @brief QC metric calculating theoretical mass of a peptide sequence

    Each PeptideHit in the FeatureMap will be annotated with its theoretical mass as metavalue 'mass'

    **/
  class OPENMS_DLLAPI PeptideMass : public QCBase
  {
  public:
    /// Constructor
    PeptideMass() = default;

    /// Destructor
    virtual ~PeptideMass() = default;

    /**
    @brief Sets the 'mass' metavalue to all PeptideHits by computing the theoretical mass

    @param features FeatureMap with PeptideHits
    **/
    void compute(FeatureMap& features);


    const String& getName() const override;

    Status requirements() const override;
  };

} // namespace OpenMS
