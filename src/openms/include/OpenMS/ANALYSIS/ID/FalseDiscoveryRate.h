// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <boost/unordered_map.hpp>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates false discovery rates (FDR) from identifications

    Either two runs of forward and decoy database identification or one run containing both (with annotations) can be used to annotate each of the peptide hits with an FDR or q-value.

    q-values are basically only adjusted p-values, also ranging from 0 to 1, with lower values being preferable.
    When looking at the list of hits ordered by q-values, then a specific q-value of @em x means that @em x*100 percent of hits with a q-value <= @em x are expected to be false positives.

    @htmlinclude OpenMS_FalseDiscoveryRate.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FalseDiscoveryRate :
    public DefaultParamHandler
  {
public:
    ///Default constructor
    FalseDiscoveryRate();

    /**
       @brief Calculates the FDR of two runs, a forward run and a decoy run on peptide level

       @param fwd_ids forward peptide identifications
       @param rev_ids reverse peptide identifications
    */
    void apply(std::vector<PeptideIdentification>& fwd_ids, std::vector<PeptideIdentification>& rev_ids) const;

    /**
        @brief Calculates the FDR of one run from a concatenated sequence DB search

        @param id peptide identifications, containing target and decoy hits
    */
    void apply(std::vector<PeptideIdentification>& id) const;

    /**
       @brief Calculates the FDR of two runs, a forward run and decoy run on protein level

       @param fwd_ids forward protein identifications
       @param rev_ids reverse protein identifications
    */
    void apply(std::vector<ProteinIdentification>& fwd_ids, std::vector<ProteinIdentification>& rev_ids) const;

    /**
       @brief Calculate the FDR of one run from a concatenated sequence db search

       @param ids protein identifications, containing target and decoy hits
    */
    void apply(std::vector<ProteinIdentification>& ids) const;

    /**
       @brief Calculate FDR on the level of molecule-query matches (e.g. peptide-spectrum matches) for "general" identification data

       @param id_data Identification data
       @param score_key Key of the score to use for FDR calculation

       @return Key of the FDR score
    */
    IdentificationData::ScoreTypeRef applyToQueryMatches(IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref) const;

  private:
    /// Not implemented
    FalseDiscoveryRate(const FalseDiscoveryRate&);

    /// Not implemented
    FalseDiscoveryRate& operator=(const FalseDiscoveryRate&);

    /// calculates the FDR, given two vectors of scores
    void calculateFDRs_(std::map<double, double>& score_to_fdr, std::vector<double>& target_scores, std::vector<double>& decoy_scores, bool q_value, bool higher_score_better) const;

    /// Helper function for applyToQueryMatches()
    void handleQueryMatch_(
        IdentificationData::QueryMatchRef match_ref,
        IdentificationData::ScoreTypeRef score_ref,
        std::vector<double>& target_scores,
        std::vector<double>& decoy_scores,
        std::map<IdentificationData::IdentifiedMoleculeRef, bool>& molecule_to_decoy,
        std::map<IdentificationData::QueryMatchRef, double>& match_to_score) const;
  };

} // namespace OpenMS

