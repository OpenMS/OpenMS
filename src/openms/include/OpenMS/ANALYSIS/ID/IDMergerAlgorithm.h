// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <unordered_set>

namespace OpenMS
{

  //TODO add params for checking consistency (i.e. how strict to check)
  //TODO add another subclass that does score-aware merging? (i.e. only keep best per peptide[sequence])

  /**
   * @brief Creates a new Protein ID run into which other runs can be inserted.
   * Creates union of protein hits but concatenates PSMs. Checks
   * search engine consistency of all inserted runs. It differs from the IDMerger tool,
   * in that it is an algorithm class and it allows inserting multiple peptide hits per
   * peptide sequence (not only the first occurrence).
   *
   * @todo allow filtering for peptide sequence to supersede the IDMerger tool.
   *       Make it keep the best PSMs though.
   */
  class OPENMS_DLLAPI IDMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    explicit IDMergerAlgorithm (const String& runIdentifier = "merged");

    /// Insert (=move and clear) a run with its peptide IDs into the internal merged data structures,
    /// based on the initial mapping from fileorigins to new run
    void insertRuns(std::vector<ProteinIdentification>&& prots,
                    std::vector<PeptideIdentification>&& peps);
    void insertRuns(const std::vector<ProteinIdentification>& prots,
                    const std::vector<PeptideIdentification>& peps);

    //TODO add methods to just insert prots or just peps. Especially makes sense if you do re-indexing anyway,
    // then you do not need the proteins. But then we need origin information. Either externally in form of a
    // String or StringList (like the one from ProteinID.getPrimaryMSRunPath). Or by having the file annotated
    // at the PeptideID (with getBasename maybe?)
    // Current solution would be to clear the ProteinIdentification if you do not need the proteins and add all the
    // necessary information about origin(s) to this ProteinIdentification.

    /// Return the merged results and reset/clear all internal data
    void returnResultsAndClear(ProteinIdentification& prots,
                   std::vector<PeptideIdentification>& peps);

  private:

    /// Returns the new identifier. The initial identifier plus a timestamp.
    String getNewIdentifier_() const;

    /// Copies over search parameters
    static void copySearchParams_(const ProteinIdentification& from, ProteinIdentification& to);

    /// Checks consistency of search engines and settings across runs before merging.
    /// @param protRuns The runs to check (first = implicit reference)
    /// @param experiment_type allow some mismatches in case of other experiment types (e.g. SILAC)
    /// @return all same? TODO: a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const String& experiment_type) const;

    /// Same as above, if you want to use a specific reference
    /// @param protRuns The runs to check (first = reference)
    /// @param ref A possibly external protein run reference
    /// @param experiment_type allow some mismatches in case of other experiment types (e.g. SILAC)
    /// @return all same? TODO: a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const ProteinIdentification& ref,
        const String& experiment_type) const;

    /// moves and inserts protein IDs if not yet present
    /// then clears the input
    void insertProteinIDs_(
        std::vector<ProteinIdentification>&& old_protRuns
    );

    /// updates the references in pepIDs to the new protein ID run
    /// then moves the peptide IDs based on the
    /// mapping in
    void updateAndMovePepIDs_(
        std::vector<PeptideIdentification>&& pepIDs,
        const std::map<String, Size>& runID_to_runIdx,
        const std::vector<StringList>& originFiles,
        bool annotate_origin
    );


    void movePepIDsAndRefProteinsToResultFaster_(
        std::vector<PeptideIdentification>&& pepIDs,
        std::vector<ProteinIdentification>&& old_protRuns
    );

    /// the resulting new Protein IDs
    ProteinIdentification prot_result_;

    /// the resulting new Peptide IDs
    std::vector<PeptideIdentification> pep_result_;

    static size_t accessionHash_(const ProteinHit& p){
      return std::hash<String>()(p.getAccession());
    }
    static bool accessionEqual_(const ProteinHit& p1, const ProteinHit& p2){
      return p1.getAccession() == p2.getAccession();
    }
    using hash_type = std::size_t (*)(const ProteinHit&);
    using equal_type = bool (*)(const ProteinHit&, const ProteinHit&);
    std::unordered_set<ProteinHit, hash_type, equal_type> collected_protein_hits_;

    /// is the resulting protein ID already filled?
    bool filled_ = false;

    /// to keep track of the mzML origins of spectra
    std::map<String, Size> file_origin_to_idx_;

    /// the new identifier string
    String id_;
  };
} // namespace OpenMS
