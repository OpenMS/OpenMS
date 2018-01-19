// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMSIMILARITY_H
#define OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMSIMILARITY_H

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithm.h>

namespace OpenMS
{
  /**
    @brief Abstract base class for ConsensusID algorithms that take peptide similarity into account.

    Similarity-based algorithms require posterior error probabilities (PEPs) as peptide scores, in order to combine scores and similarities into a consensus score for each peptide. See the following publication for the formula governing this calculation:

    Nahnsen <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr2002879">Probabilistic consensus scoring improves tandem mass spectrometry peptide identification</a> (J. Proteome Res., 2011, PMID: 21644507).

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
                        SequenceGrouping& results) override;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMSIMILARITY_H
