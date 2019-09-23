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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>

namespace OpenMS
{
  /** \brief Algorithm class that implements simple protein inference by aggregation of peptide scores.
   * It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
   * and if you want to use shared peptides ("use_shared_peptides").
   * First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
   * Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
   * - Modifications (modified sequence string)
   * - Charge states
   * The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
   * Possible aggregation methods that can be set via the parameter "aggregation_method" are:
   * - "maximum" (default)
   * - "sum"
   * - "product" (ignoring zeroes)
   * Annotation of the number of peptides used for aggregation can be disabled (see parameters).
   * Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
   */
  class OPENMS_DLLAPI BasicProteinInferenceAlgorithm :
    public DefaultParamHandler,
    public ProgressLogger
  {
    public:

    /**
     * @brief The aggregation method
     */
    enum class AggregationMethod
    {
      PROD, ///< aggregate by product (ignore zeroes)
      SUM, ///< aggregate by summing
      MAXIMUM ///< aggregate by maximum
    };

    /// Default constructor
    BasicProteinInferenceAlgorithm();

    /// main method of BasicProteinInferenceAlgorithm
    /// inputs are not const, since it will get annotated with results
    /// annotation of protein groups is currently only possible for a single protein ID run
    void run(std::vector<PeptideIdentification> &pep_ids, std::vector<ProteinIdentification> &prot_ids) const;
    void run(std::vector<PeptideIdentification> &pep_ids, ProteinIdentification &prot_id) const;

  private:

    /**
     * @brief Performs simple inference on one protein run.
     * @param acc_to_protein_hitP_and_count Maps Accessions to a pair of ProteinHit pointers
     *  and number of peptidoforms encountered @Todo could use member as hash to save strings
     * @param best_pep Maps (un)modified peptide sequence to a map from charge (0 when unconsidered) to the
     *  best PeptideHit pointer
     * @param prot_run The current run to process
     * @param pep_ids Peptides for the current run to process
     */
    void processRun_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      std::unordered_map<std::string, std::map<Int, PeptideHit*>>& best_pep,
      ProteinIdentification& prot_run,
      std::vector<PeptideIdentification>& pep_ids,
      Size min_peptides_per_protein) const;
  };
} //namespace OpenMS
