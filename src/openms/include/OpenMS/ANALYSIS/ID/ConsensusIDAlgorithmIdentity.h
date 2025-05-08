// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>

namespace OpenMS
{
  /**
    @brief Abstract base class for ConsensusID algorithms that compare only identical sequences.

    Search engine scores are grouped by peptide sequence in apply_(). For each sequence, getAggregateScore_() is called to produce a consensus score from the list of search engine scores.

    All derived classes should implement getAggregateScore_(). They may re-implement preprocess_() if required.

    @htmlinclude OpenMS_ConsensusIDAlgorithmIdentity.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmIdentity :
    public ConsensusIDAlgorithm
  {
  protected:
    /// Default constructor
    ConsensusIDAlgorithmIdentity();

    /**
       @brief Preprocessing of peptide IDs (in the beginning of "apply_").

       Checks whether the score types are the same (warns if not) and whether the score orientations agree (error if not).

       @param ids Input/output peptide identifications

       @throw Exception::InvalidValue Score orientations do not agree
    */
    virtual void preprocess_(std::vector<PeptideIdentification>& ids);

    /**
       @brief Aggregate peptide scores into one final score (to be implemented by subclasses).

       @param scores List of scores for the same peptide by different search engines
       @param higher_better Whether higher or lower scores are better

       @return Final score for the respective peptide
    */
    virtual double getAggregateScore_(std::vector<double>& scores,
                                      bool higher_better) = 0;

  private:
    /// Not implemented
    ConsensusIDAlgorithmIdentity(const ConsensusIDAlgorithmIdentity&);

    /// Not implemented
    ConsensusIDAlgorithmIdentity& operator=(const ConsensusIDAlgorithmIdentity&);

    /// Consensus scoring
    void apply_(std::vector<PeptideIdentification>& ids,
        const std::map<String, String>& se_info,
        SequenceGrouping& results) override;
  };

} // namespace OpenMS

