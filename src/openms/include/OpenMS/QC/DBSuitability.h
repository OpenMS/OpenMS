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
  * @brief This class holds the functionality of calculating the database suitability.
  *
  * To calculate the suitability of a database for a specific mzML for identification search, it
  * is vital to perform a combined deNovo+database identification search. Meaning that the database
  * should be appended with an additional entry derived from concatenated deNovo sequences from said mzML.
  * Currently only Comet search is supported.
  *
  * This class will calculate q-values by itself and will throw an error if any q-value calculation was done beforehand.
  *
  * The algorithm parameters can be set using setParams().
  *
  * Allows for multiple usage of the compute function. The result of each call is stored internally in a vector.
  * Therefore old results will not be overridden by a new call. This vector then can be returned using getResults().
  *
  * This class serves as the library representation of @ref TOPP_DatabaseSuitability
  */
  class OPENMS_DLLAPI DBSuitability:
    public DefaultParamHandler
  {
  public:
    /// struct to store results
    struct SuitabilityData
    {
      /// number of times the top hit is considered to be a deNovo hit
      Size num_top_novo = 0;

      /// number of times the top hit is considered to be a database hit
      Size num_top_db = 0;
      
      /// number of times a deNovo hit scored on top of a database hit
      Size num_interest = 0;

      /// number of times a deNovo hit scored on top of a database hit,
      /// but their score difference was small enough, that it was still counted as a database hit
      Size num_re_ranked = 0;

      /// the cut-off that was used to determine when a score difference was "small enough"
      /// this is normalized by mw
      double cut_off = DBL_MAX;

      /// the suitability of the database used for identification search, calculated with:
      ///               #db_hits / (#db_hits + #deNovo_hit)
      /// can reach from 0 -> the database was not at all suited to 1 -> the perfect database was used
      ///
      /// Preliminary tests have shown that databases of the right organism or close related organisms
      /// score around 0.9 to 0.95, organisms from the same class can still score around 0.8, organisms
      /// from the same phylum score around 0.5 to 0.6 and after that it quickly falls to suitabilities
      /// of 0.15 or even 0.05.
      /// Note that these test were only performed for one mzML and your results might differ.
      double suitability = 0;
    };

    /// Constructor
    /// Settings are initialized with their default values:
    /// no_rerank = false, reranking_cutoff_percentile = 1, FDR = 0.01
    DBSuitability();

    /// Destructor
    ~DBSuitability() = default;

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
    * Parameters can be set using the functionality of DefaultParamHandler.
    * Parameters are:
    *           no_rerank                   - re-ranking can be turned off with this
    *           reranking_cutoff_percentile - percentile that determines which cut-off will be returned
    *           FDR                         - q-value that should be filtered for
    *                                         Preliminary tests have shown that database suitability
    *                                         is rather stable across common FDR thresholds from 0 - 5 %
    *
    * Since q-values need to be calculated the identifications are taken by copy.
    *
    * Result is appended to the result member. This allows for multiple usage.
    *
    * @param pep_ids      vector containing pepIDs with target/decoy annotation coming from a deNovo+database 
    *                     identification search (currently only Comet-support) without FDR
    *                     vector is modified internally, and is thus copied
    * @throws             MissingInformation if no target/decoy annotation is found
    * @throws             MissingInformation if no xcorr is found
    * @throws             Precondition if a q-value is found in the input
    */
    void compute(std::vector<PeptideIdentification> pep_ids);

    /**
    * @brief Returns results calculated by this metric
    *
    * The returned vector contains one DBSuitabilityData object for each time compute was called.
    * Each of these objects contains the suitability information that was extracted from the
    * identifications used for the corresponding call of compute.
    *
    * @returns  DBSuitabilityData objects in a vector
    */
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
    * Decoy differences of all N pepIDs are calculated. The (1-reranking_cutoff_percentile)*N highest
    * one is returned.
    * It is asssumed that this difference accounts for 'reranking_cutoff_percentile' of the re-ranking cases.
    *
    * @param pep_ids                      vector containing the pepIDs
    * @param reranking_cutoff_percentile  percentile that determines which cut-off will be returned
    * @returns                            xcorr cut-off
    * @throws                             IllegalArgument if reranking_cutoff_percentile isn't in range [0,1]
    * @throws                             IllegalArgument if reranking_cutoff_percentile is too low for a decoy cut-off to be calculated
    * @throws                             MissingInformation if no more than 20 % of the peptide IDs have two decoys in their top ten peptide hits
    */
    double getDecoyCutOff_(const std::vector<PeptideIdentification>& pep_ids, double reranking_cutoff_percentile);

    /**
    * @brief Tests if a PeptideHit is considered a deNovo hit
    *
    * To test this the function looks into the protein accessions.
    * If only the deNovo protein is found, 'true' is returned.
    * If at least one database protein is found, 'false' is returned.
    *
    * @param hit      PepHit in question
    * @returns        true/false
    */
    bool isNovoHit_(const PeptideHit& hit);

    /**
    * @brief Tests if a PeptideHit has a lower q-value than the given FDR threshold, i.e. passes FDR
    *
    * Q-value is searched at score and at meta-value level.
    *
    * @param hit            PepHit in question
    * @param FDR            FDR threshold to check against
    * @returns              true/false
    */
    bool passesFDR_(const PeptideHit& hit, double FDR);
  };
}

