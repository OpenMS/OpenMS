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

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMIDENTITY_H
#define OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMIDENTITY_H

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
                        SequenceGrouping& results) override;
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMIDENTITY_H
