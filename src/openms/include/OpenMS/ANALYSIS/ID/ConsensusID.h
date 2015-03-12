// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Sven Nahnsen $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSID_H
#define OPENMS_ANALYSIS_ID_CONSENSUSID_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates a consensus ID from several ID runs

    This class combines several ID runs using one of several, available algorithms.

        @htmlinclude OpenMS_ConsensusID.parameters

        @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusID :
    public DefaultParamHandler
  {
public:
    ///Default constructor
    ConsensusID();

    /**
        @brief Calculates the consensus ID for a set of PeptideIdentification instances of the same spectrum

        @note Make sure that the score orientation (PeptideIdentification::isHigherScoreBetter()) is set properly!
    */
    void apply(std::vector<PeptideIdentification>& ids);

private:
    /// mapping: peptide sequence -> (charge, scores)
    typedef std::map<AASequence, std::pair<Int, std::vector<double> > > 
      SequenceGrouping;

    /// Not implemented
    ConsensusID(const ConsensusID&);

    /// Not implemented
    ConsensusID& operator=(const ConsensusID&);

    /// consensus based on ranks
    void rank_(std::vector<PeptideIdentification>& ids);

    /// consensus based on average score
    void average_(std::vector<PeptideIdentification>& ids);

    /// consensus based on PEPs and scoring matrix
    void PEPMatrix_(std::vector<PeptideIdentification>& ids);

    /// consensus based on PEPs and ion similarity
    void PEPIons_(std::vector<PeptideIdentification>& ids);

    /// consensus based on best score
    void best_(std::vector<PeptideIdentification>& ids);

    /// helper function for 'best_' and 'average_'
    void groupHits_(std::vector<PeptideIdentification>& ids,
                    SequenceGrouping& grouping);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSID_H
