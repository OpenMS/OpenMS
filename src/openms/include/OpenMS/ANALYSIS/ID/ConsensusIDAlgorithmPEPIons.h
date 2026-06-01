// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>

namespace OpenMS
{
  /**
    @brief Calculates a consensus from multiple ID runs based on PEPs and shared ions.

    @htmlinclude OpenMS_ConsensusIDAlgorithmPEPIons.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmPEPIons :
    public ConsensusIDAlgorithmSimilarity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmPEPIons();

  private:
    /// Fragment mass tolerance (for "PEPIons_")
    double mass_tolerance_;

    /// Min. number of shared fragments (for "PEPIons")
    Size min_shared_;

    /// Not implemented
    ConsensusIDAlgorithmPEPIons(const ConsensusIDAlgorithmPEPIons&);

    /// Not implemented
    ConsensusIDAlgorithmPEPIons& operator=(const ConsensusIDAlgorithmPEPIons&);

    /// Docu in base class
    void updateMembers_() override;

    /// Sequence similarity based on matching ions
    double getSimilarity_(AASequence seq1, AASequence seq2) override;

  };

} // namespace OpenMS

