// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Abstract base class for all ConsensusID algorithms (that calculate a consensus from multiple ID runs).

    The main function is apply(), which aggregates several peptide identifications into one.

    Derived classes should implement apply_(), which takes a list of peptide identifications and produces a map of peptide sequences with accompanying scores (and charge states).
    Currently there are two derived classes, OpenMS::ConsensusIDAlgorithmIdentity and OpenMS::ConsensusIDAlgorithmSimilarity. They serve as abstract base classes for algorithms that score only
    identical peptide sequences together and algorithms that take similarities between peptides into account, respectively.

    See also the documentation of the TOPP tool, @ref TOPP_ConsensusID, for more information (e.g. on the @p filter: parameters).

    @htmlinclude OpenMS_ConsensusIDAlgorithm.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithm : public DefaultParamHandler
  {
  public:
    /**
        @brief Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature.

        Make sure that the score type (PeptideIdentification::getScoreType()) and the score orientation (PeptideIdentification::isHigherScoreBetter()) are set properly!

        @param ids Peptide identifications (input: more than one, output: one)
        @param number_of_runs Number of ID runs (default: size of "ids")
        @param se_info map from run identifiers to search engine infos to retain original search engine information
        @todo we could pass the score_types that we want to carry over in the map as well (right now it always takes main)
     */
    void apply(std::vector<PeptideIdentification>& ids, const std::map<String, String>& se_info, Size number_of_runs = 0);

    void apply(std::vector<PeptideIdentification>& ids, Size number_of_runs = 0);

    /// Virtual destructor
    ~ConsensusIDAlgorithm() override;

  protected:
    struct HitInfo {
      Int charge;
      std::vector<double> scores;
      std::vector<String> types;
      // in case too much information is stored, TD and evidence
      // could be re-annotated with PeptideIndexer later
      String target_decoy;
      std::set<PeptideEvidence> evidence;
      double final_score;
      double support;
      // TODO: we could gather spectrum_refs here as well,
      //  to support passing of spectrum_ref if ALL refs of a group are the same
      //  For now, we do it in the ConsensusID TOPP tool class in cases where we
      //  know that refs will be the same.
    };

    /// Mapping: peptide sequence -> (charge, scores)
    typedef std::map<AASequence, HitInfo> SequenceGrouping;

    /// Number of peptide hits considered per ID run (input parameter)
    Size considered_hits_;

    /// Number of ID runs
    Size number_of_runs_;

    /// Fraction of required support by other ID runs (input parameter)
    double min_support_;

    /// Count empty runs in "min_support" calculation? (input parameter)
    bool count_empty_;

    /// Keep old scores?
    bool keep_old_scores_;

    /// Default constructor
    ConsensusIDAlgorithm();

    /**
       @brief Consensus computation (to be implemented by subclasses).

       @param ids Peptide identifications (input)
       @param se_info mapping from run identifier to search engine to carry over infos to result
       @param results Algorithm results (output). For each peptide sequence, two scores are expected: the actual consensus score and the "support" value, in this order.
    */
    virtual void apply_(std::vector<PeptideIdentification>& ids, const std::map<String, String>& se_info, SequenceGrouping& results) = 0;

    /// Docu in base class
    void updateMembers_() override;

    /// Compare (and possibly update) charge state information
    void compareChargeStates_(Int& recorded_charge, Int new_charge, const AASequence& peptide);

  private:
    /// Not implemented
    ConsensusIDAlgorithm(const ConsensusIDAlgorithm&) = delete;

    /// Not implemented
    ConsensusIDAlgorithm& operator=(const ConsensusIDAlgorithm&) = delete;
  };

} // namespace OpenMS
