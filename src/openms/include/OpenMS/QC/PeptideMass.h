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
    @brief QC metric annotating the OBSERVED neutral mass of a peptide from the precursor m/z (precursor_mz - proton_mass)*charge, NOT the theoretical sequence mass

    The first PeptideHit of each PeptideIdentification in the FeatureMap is annotated with its observed
    neutral mass (derived from the precursor m/z) as metavalue 'mass'.

    **/
  class OPENMS_DLLAPI PeptideMass : public QCBase
  {
  public:
    /// Constructor
    PeptideMass() = default;

    /// Destructor
    virtual ~PeptideMass() = default;

    /**
    @brief Sets the 'mass' metavalue on the first PeptideHit of each ID to the observed neutral mass from the precursor m/z ((precursor_mz - proton_mass)*charge)

    @param[in,out] features FeatureMap with PeptideHits
    **/
    void compute(FeatureMap& features);


    const std::string& getName() const override;

    Status requirements() const override;
  };

} // namespace OpenMS
