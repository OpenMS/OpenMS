// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <cfloat>
#include <map>
#include <vector>

namespace OpenMS
{
  class PeptideIdentification;
  class PeptideHit;

  /**
  * @brief This class serves as the library representation @ref TOPP_DatabaseSuitability
  *
  * This class holds the functionality of calculating the database suitability.
  * This can only be done if a combined deNovo+database identification search was performed.
  * Currently only Comet search is supported.
  *
  * Allows for multiple usage, because results are stored internally and can be returned
  * using getResults().
  *
  */
  class OPENMS_DLLAPI Suitability:
    public DefaultParamHandler
  {
  public:

    /// struct to store results
    struct SuitabilityData
    {
      Size num_top_novo = 0;
      Size num_top_db = 0;
      Size num_re_ranked = 0;
      Size num_interest = 0;
      double cut_off = DBL_MAX;
      double suitability = 0;
    };

    /// Constructor
    /// Settings are initialized with their default values:
    /// no_re_rank = false, novo_fract = 1, FDR = 0.01
    Suitability();

    /// Destructor
    ~Suitability() = default;

    /**
    * @brief Computes suitability of a database used to search a mzML
    *
    * Counts top deNovo and top database hits. The ratio of db hits vs
    * all hits yields the suitability.
    * To re-rank cases, where a de novo peptide scores just higher than
    * the database peptide, a decoy cut-off is calculated. This functionality
    * can be turned off. This will result in an underestimated suitability,
    * but it can solve problems like different search engines or to few decoy hits.
    *
    * Result is appended to the result member. This allows for multiple usage.
    *
    * @param pep_ids      vector containing pepIDs coming from a deNovo+database 
    *                     identification search (currently only Comet-support)
    * @throws             MissingInformation if decoy cut-off could not be calculated
    * @throws             MissingInformation if no target/decoy annotation is found
    * @throws             MissingInformation if no xcorr is found
    * @throws             Precondition if FDR wasn't calculated
    */
    void compute(std::vector<PeptideIdentification>& pep_ids);

    /// return results
    const std::vector<SuitabilityData>& getResults() const;

  private:
    /// result vector
    std::vector<SuitabilityData> results_;

    /**
    * @brief Calculates the xcorr difference between the top two hits marked as decoy
    *
    * Only searches the top ten hits for two decoys. If there aren't two decoys, DBL_MAX
    * is returned.
    *
    * @param pep_id     pepID from where the decoy difference will be calculated
    * @returns          xcorr difference
    * @throws           MissingInformation if no target/decoy annotation is found
    * @throws           MissingInformation if no xcorr is found
    */
    double getDecoyDiff_(const PeptideIdentification& pep_id);

    /**
    * @brief Calculates a xcorr cut-off based on decoy hits
    *
    * Decoy differences of all N pepIDs are calculated. The (1-novo_fract)*N highest
    * one is returned.
    * It is asssumed that this difference accounts for novo_fract of the re-ranking cases.
    *
    * @param pep_ids      vector containing the pepIDs
    * @param novo_fract   fraction of how many cases, where a de novo peptide scores just higher than the database peptide, will be re-rank
    * @returns            xcorr cut-off
    * @throws             MissingInformation if no more than 20 % of the pepIDs have two decoys in there top ten
    */
    double getDecoyCutOff_(const std::vector<PeptideIdentification>& pep_ids, double novo_fract);

    /**
    * @brief Tests if a PeptideHit is considered a deNovo hit
    *
    * To test this the function looks into the protein accessions.
    * If only the deNovo protein is found, 'true' is returned.
    * If one database protein is found, 'false' is returned.
    *
    * @param hit      PepHit in question
    * @returns        true/false
    */
    bool isNovoHit_(const PeptideHit& hit);

    /**
    * @brief Tests if a PeptideHit has a higher q-value the given FDR
    *
    * Q-value is searched at score and at meta-value level.
    *
    * @param hit            PepHit in question
    * @param FDR            FDR to check against
    * @returns              true/false
    */
    bool scoreHigherThanFDR_(const PeptideHit& hit, double FDR);
  };
}

