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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <unordered_set>

namespace OpenMS
{

  //TODO add params for checking consistency (i.e. how strict to check)
  //TODO add another subclass that does score-aware merging? (i.e. only keep best per peptide[sequence])

  /**
   * @brief Algorithm for merging multiple protein and peptide identification runs.
   *
   * This class creates a new Protein ID run into which other runs can be inserted.
   * It performs the following operations:
   * - Creates a union of protein hits from all inserted runs
   * - Concatenates Peptide-Spectrum Matches (PSMs) from all runs
   * - Checks search engine consistency across all inserted runs
   * - Maintains references between peptide IDs and their corresponding protein IDs
   * 
   * The algorithm differs from the IDMerger tool in two key aspects:
   * 1. It is implemented as an algorithm class rather than a tool
   * 2. It allows inserting multiple peptide hits per peptide sequence (not only the first occurrence)
   * 
   * The class handles the complexity of merging identification data from different sources
   * while ensuring consistency and maintaining proper references between proteins and peptides.
   * It can be used in workflows where identification results from multiple files or runs
   * need to be combined into a single comprehensive result set.
   * 
   * The algorithm can optionally annotate the origin of each identification to maintain
   * traceability of the merged results back to their source files.
   * 
   * @see IDMerger
   * 
   * @todo Allow filtering for peptide sequence to supersede the IDMerger tool.
   *       Make it keep the best PSMs though.
   */
  class OPENMS_DLLAPI IDMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    /**
     * @brief Constructor for the IDMergerAlgorithm.
     * 
     * Initializes a new merger with the specified run identifier.
     * 
     * @param runIdentifier Base identifier for the merged run (default: "merged")
     * @param addTimeStampToID Whether to append a timestamp to the run identifier for uniqueness (default: true)
     */
    explicit IDMergerAlgorithm (const String& runIdentifier = "merged", bool addTimeStampToID = true);

    /**
     * @brief Insert runs using move semantics.
     * 
     * Inserts (moves and clears) protein and peptide identifications into the internal 
     * merged data structures. This version uses move semantics for better performance
     * when the source data is no longer needed.
     * 
     * @param prots Vector of protein identifications to be merged
     * @param peps Vector of peptide identifications to be merged
     */
    void insertRuns(std::vector<ProteinIdentification>&& prots,   
                    PeptideIdentificationList&& peps);
                    
    /**
     * @brief Insert runs using copy semantics.
     * 
     * Inserts (copies) protein and peptide identifications into the internal 
     * merged data structures. This version preserves the source data.
     * 
     * @param prots Vector of protein identifications to be merged
     * @param peps Vector of peptide identifications to be merged
     */
    void insertRuns(const std::vector<ProteinIdentification>& prots,
                    const PeptideIdentificationList& peps);

    //TODO add methods to just insert prots or just peps. Especially makes sense if you do re-indexing anyway,
    // then you do not need the proteins. But then we need origin information. Either externally in form of a
    // String or StringList (like the one from ProteinID.getPrimaryMSRunPath). Or by having the file annotated
    // at the PeptideID (with getBasename maybe?)
    // Current solution would be to clear the ProteinIdentification if you do not need the proteins and add all the
    // necessary information about origin(s) to this ProteinIdentification.

    /**
     * @brief Return the merged results and reset internal state.
     * 
     * Retrieves the merged protein and peptide identifications and clears all internal
     * data structures, preparing the algorithm instance for potential reuse.
     * 
     * This method should be called after all desired runs have been inserted to obtain
     * the final merged result.
     * 
     * @param prots [out] The merged protein identification containing the union of all protein hits
     * @param peps [out] The merged peptide identifications containing all PSMs from the inserted runs
     * 
     * @note After calling this method, the internal state is reset, and the algorithm
     *       can be reused for a new merging operation.
     */
    void returnResultsAndClear(ProteinIdentification& prots,
                   PeptideIdentificationList& peps);

  private:

    /**
     * @brief Generate a new identifier for the merged run.
     * 
     * Creates a new identifier by combining the base identifier with a timestamp
     * if requested.
     * 
     * @param addTimeStampToID Whether to append a timestamp to the identifier
     * @return The generated identifier string
     */
    String getNewIdentifier_(bool addTimeStampToID) const;

    /**
     * @brief Copy search parameters between protein identifications.
     * 
     * Transfers search parameters from one protein identification to another.
     * 
     * @param from Source protein identification
     * @param to Destination protein identification
     */
    static void copySearchParams_(const ProteinIdentification& from, ProteinIdentification& to);

    /**
     * @brief Check consistency of search engines and settings across runs.
     * 
     * Verifies that all runs have compatible search engine settings before merging.
     * Uses the first run as an implicit reference.
     * 
     * @param protRuns The runs to check (first = implicit reference)
     * @param experiment_type Experiment type to allow certain mismatches (e.g., "SILAC")
     * @return True if all runs are consistent, false otherwise
     * @throws BaseException for disagreeing settings
     * 
     * @todo Return a merged RunDescription about what to put in the new runs (e.g., for SILAC)
     */
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const String& experiment_type) const;

    /**
     * @brief Check consistency of search engines and settings against a reference.
     * 
     * Verifies that all runs have compatible search engine settings before merging,
     * using an explicitly provided reference run.
     * 
     * @param protRuns The runs to check
     * @param ref An external protein run to use as reference
     * @param experiment_type Experiment type to allow certain mismatches (e.g., "SILAC")
     * @return True if all runs are consistent with the reference, false otherwise
     * @throws BaseException for disagreeing settings
     * 
     * @todo Return a merged RunDescription about what to put in the new runs (e.g., for SILAC)
     */
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const ProteinIdentification& ref,
        const String& experiment_type) const;

    /**
     * @brief Insert protein identifications into the merged result.
     * 
     * Moves and inserts protein IDs if not yet present, then clears the input.
     * 
     * @param old_protRuns Vector of protein identifications to insert
     */
    void insertProteinIDs_(
        std::vector<ProteinIdentification>&& old_protRuns
    );

    /**
     * @brief Update peptide ID references and move them to the result.
     * 
     * Updates the references in peptide IDs to point to the new protein ID run,
     * then moves the peptide IDs based on the provided mapping.
     * 
     * @param pepIDs Vector of peptide identifications to update and move
     * @param runID_to_runIdx Mapping from run IDs to run indices
     * @param originFiles List of origin files for each run
     * @param annotate_origin Whether to annotate peptide IDs with their origin
     */
    void updateAndMovePepIDs_(
        PeptideIdentificationList&& pepIDs,
        const std::map<String, Size>& runID_to_runIdx,
        const std::vector<StringList>& originFiles,
        bool annotate_origin
    );

    /**
     * @brief Optimized method to move peptide IDs and reference proteins to result.
     * 
     * A faster implementation for moving peptide IDs and their referenced proteins
     * to the result data structures.
     * 
     * @param pepIDs Vector of peptide identifications to move
     * @param old_protRuns Vector of protein identifications to reference
     */
    void movePepIDsAndRefProteinsToResultFaster_(
        PeptideIdentificationList&& pepIDs,
        std::vector<ProteinIdentification>&& old_protRuns
    );

    /// The resulting merged protein identification
    ProteinIdentification prot_result_;

    /// The resulting merged peptide identifications
    PeptideIdentificationList pep_result_;

    /**
     * @brief Hash function for protein hits based on accession.
     * 
     * @param p Protein hit to hash
     * @return Hash value for the protein hit
     */
    static size_t accessionHash_(const ProteinHit& p){
      return std::hash<String>()(p.getAccession());
    }
    
    /**
     * @brief Equality function for protein hits based on accession.
     * 
     * @param p1 First protein hit to compare
     * @param p2 Second protein hit to compare
     * @return True if the accessions are equal, false otherwise
     */
    static bool accessionEqual_(const ProteinHit& p1, const ProteinHit& p2){
      return p1.getAccession() == p2.getAccession();
    }
    
    /// Type alias for the hash function
    using hash_type = std::size_t (*)(const ProteinHit&);
    
    /// Type alias for the equality function
    using equal_type = bool (*)(const ProteinHit&, const ProteinHit&);
    
    /// Set of collected protein hits using custom hash and equality functions
    std::unordered_set<ProteinHit, hash_type, equal_type> collected_protein_hits_;

    /// Flag indicating whether the resulting protein ID is already filled
    bool filled_ = false;

    /// Mapping to keep track of the mzML origins of spectra
    std::map<String, Size> file_origin_to_idx_;

    /// The new identifier string for the merged run
    String id_;

    /// Flag indicating whether the identifier should be fixed (i.e., not contain a timestamp)
    bool fixed_identifier_;
  };
} // namespace OpenMS
