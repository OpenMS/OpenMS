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
    @brief Abstract base class for ConsensusID algorithms that take peptide similarity into account.

    Similarity-based algorithms require posterior error probabilities (PEPs) as peptide scores, in order to combine scores and similarities into a consensus score for each peptide. See the following publication for the formula governing this calculation:

    Nahnsen <em>et al.</em>: <a href="https://doi.org/10.1021/pr2002879">Probabilistic consensus scoring improves tandem mass spectrometry peptide identification</a> (J. Proteome Res., 2011, PMID: 21644507).

    Derived classes should implement getSimilarity_(), which defines how similarity of two peptide sequences is quantified.

    @htmlinclude OpenMS_ConsensusIDAlgorithmSimilarity.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmSimilarity :
    public ConsensusIDAlgorithm
  {
  protected:
    /// Default constructor
    ConsensusIDAlgorithmSimilarity();

    /// Mapping: pair of peptide sequences -> sequence similarity
    typedef std::map<std::pair<AASequence, AASequence>, double> SimilarityCache;

    /// Cache for already computed sequence similarities
    SimilarityCache similarities_;

    /**
       @brief Sequence similarity calculation (to be implemented by subclasses).

       Implementations should use/update the cache of previously computed similarities.

       @return Similarity between two sequences in the range [0, 1]
    */
    virtual double getSimilarity_(AASequence seq1, AASequence seq2) = 0;

  private:
    /// Not implemented
    ConsensusIDAlgorithmSimilarity(const ConsensusIDAlgorithmSimilarity&);

    /// Not implemented
    ConsensusIDAlgorithmSimilarity& operator=(const ConsensusIDAlgorithmSimilarity&);

    /// Consensus scoring
    void apply_(std::vector<PeptideIdentification>& ids,
        const std::map<String, String>& se_info,
        SequenceGrouping& results) override;
  };

} // namespace OpenMS

