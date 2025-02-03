// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak2D.h>

#include <vector>

namespace OpenMS
{

  class PeptideHit;
  class ConsensusMap;

  /**
    @brief [experimental class] given a peptide quantitation, infer corresponding protein quantities

    Infers protein ratios from peptide ratios (currently using unique peptides only).
    Use the IDMapper class to add protein and peptide information to a
    quantitative ConsensusMap prior to this step.
  */
  class OPENMS_DLLAPI ProteinInference
  {

public:

    typedef Peak2D::IntensityType IntensityType;

    /// Constructor
    ProteinInference();

    /// copy constructor
    ProteinInference(const ProteinInference& cp);

    /// assignment operator
    ProteinInference& operator=(const ProteinInference& rhs);

    /**
        @brief given a peptide quantitation, infer corresponding protein quantities

        Infers protein ratios from peptide ratios (currently using unique peptides only).
        Use the IDMapper class to add protein and peptide information to a
        quantitative ConsensusMap prior to this step.

        @param consensus_map Peptide quantitation with ProteinIdentifications attached, where
                     Protein quantitation will be attached
        @param reference_map Index of (iTRAQ) reference channel within the consensus map

        @throws Exception::MissingInformation if Protein/PeptideIdentifications are missing
    */
    void infer(ConsensusMap& consensus_map, const UInt reference_map);


protected:

    void infer_(ConsensusMap& consensus_map,
                const size_t protein_idenfication_index,
                const UInt reference_map);

    bool sortByUnique_(std::vector<PeptideHit>& peptide_hits_local, const bool is_higher_score_better);

  }; // !class

} // !namespace

