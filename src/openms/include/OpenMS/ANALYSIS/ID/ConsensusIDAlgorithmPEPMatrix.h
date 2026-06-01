// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>
#include <OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h>

namespace OpenMS
{
  /**
    @brief Calculates a consensus from multiple ID runs based on PEPs and sequence similarities.

    @note The similarity scoring is based on an amino acid substitution matrix. Therefore only the raw amino acid sequences, without post-translational modifications (PTMs), can be considered for similarity scoring - PTMs are ignored during this step. However, PTMs on peptides are retained and separate results are produced for differently-modified peptides.

    @htmlinclude OpenMS_ConsensusIDAlgorithmPEPMatrix.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmPEPMatrix :
    public ConsensusIDAlgorithmSimilarity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmPEPMatrix();


  private:

    /// object for alignment score calculation
    NeedlemanWunsch alignment_;

    /// Not implemented
    ConsensusIDAlgorithmPEPMatrix(const ConsensusIDAlgorithmPEPMatrix&);

    /// Not implemented
    ConsensusIDAlgorithmPEPMatrix& operator=(const ConsensusIDAlgorithmPEPMatrix&);

    /// Sequence similarity based on substitution matrix (ignores PTMs)
    double getSimilarity_(AASequence seq1, AASequence seq2) override;

    // Docu in base class
    void updateMembers_() override;
  };

} // namespace OpenMS

