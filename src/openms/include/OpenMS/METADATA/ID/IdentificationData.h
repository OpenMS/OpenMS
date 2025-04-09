// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/ProcessingStep.h>
#include <OpenMS/METADATA/ID/Observation.h>
#include <OpenMS/METADATA/ID/DBSearchParam.h>
#include <OpenMS/METADATA/ID/IdentifiedCompound.h>
#include <OpenMS/METADATA/ID/IdentifiedSequence.h>
#include <OpenMS/METADATA/ID/InputFile.h>
#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/ID/ParentMatch.h>
#include <OpenMS/METADATA/ID/ObservationMatch.h>
#include <OpenMS/METADATA/ID/ParentSequence.h>
#include <OpenMS/METADATA/ID/ParentGroup.h>
#include <OpenMS/METADATA/ID/ObservationMatchGroup.h>
#include <OpenMS/METADATA/ID/ScoreType.h>

#include <unordered_set>

namespace OpenMS
{
  /*!
    @brief Representation of spectrum identification results and associated data

    This class provides capabilities for storing spectrum identification results from different
    types of experiments/molecules (proteomics: peptides/proteins, metabolomics: small molecules, "nucleomics": RNA).

    The class design has the following goals:
    - Provide one structure for storing all relevant data for spectrum identification results.
    - Store data non-redundantly.
    - Ensure consistency (e.g. no conflicting information; no "dangling references").
    - Allow convenient and efficient querying.
    - Support different types of experiments, as mentioned above, in one common framework.

    The following important subordinate classes are provided to represent different types of data:
    <table>
    <tr><th>Class <th>Represents <th>Key <th>Proteomics example <th>Corresponding legacy class
    <tr><td>ProcessingStep <td>Information about a data processing step that was applied (e.g. input files, software used, parameters) <td>Combined information <td>Mascot search <td>ProteinIdentification
    <tr><td>Observation <td>A search query (with identifier, RT, m/z) from an input file, i.e. an MS2 spectrum or feature (for accurate mass search) <td>File/Identifier <td>MS2 spectrum <td>PeptideIdentification
    <tr><td>ParentSequence <td>An entry in a FASTA file with associated information (sequence, coverage, etc.) <td>Accession <td>Protein <td>ProteinHit
    <tr><td>IdentifiedPeptide/-Oligo/-Compound <td>An identified molecule of the respective type <td>Sequence (or identifier for a compound) <td>Peptide <td>PeptideHit
    <tr><td>ObservationMatch <td>A match between a query (Observation), identified molecule (Identified...), and optionally adduct <td>Combination of query/molecule/adduct references <td>Peptide-spectrum match (PSM) <td>PeptideIdentification/PeptideHit
    </table>

    To populate an IdentificationData instance with data, "register..." functions are used.
    These functions return "references" (implemented as iterators) that can be used to refer to stored data items and thus form connections.
    For example, a protein can be stored using registerParentSequence, which returns a corresponding reference.
    This reference can be used to build an IdentifiedPeptide object that references the protein.
    An identified peptide referencing a protein can only be registered if that protein has been registered already, to ensure data consistency.
    Given the identified peptide, information about the associated protein can be retrieved efficiently by simply dereferencing the reference.

    To ensure non-redundancy, many data types have a "key" (see table above) to which a uniqueness constraint applies.
    This means only one item of such a type with a given key can be stored in an IdentificationData object.
    If items with an existing key are registered subsequently, attempts are made to merge new information (e.g. additional scores) into the existing entry.
    The details of this merging are handled in the @p merge function in each data class.

    @warning This class is not thread-safe while being modified.

    @ingroup Metadata
  */

  /// Remove elements from a set (or ordered multi_index_container) if they fulfill a predicate (TODO: deprecate and use std::erase_if with C++20 adoption)
  template <typename ContainerType, typename PredicateType>
  static void removeFromSetIf_(ContainerType& container, PredicateType predicate)
  {
    for (auto it = container.begin(); it != container.end(); )
    {
      if (predicate(it))
      {
        it = container.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }

  class OPENMS_DLLAPI IdentificationData: public MetaInfoInterface
  {
  public:

    // to be able to add overloads and still find the inherited ones
    using MetaInfoInterface::setMetaValue;

    // type definitions:
    using MoleculeType = IdentificationDataInternal::MoleculeType;
    using MassType = IdentificationDataInternal::MassType;

    using InputFile = IdentificationDataInternal::InputFile;
    using InputFiles = IdentificationDataInternal::InputFiles;
    using InputFileRef = IdentificationDataInternal::InputFileRef;

    using ProcessingSoftware =
      IdentificationDataInternal::ProcessingSoftware;
    using ProcessingSoftwares =
      IdentificationDataInternal::ProcessingSoftwares;
    using ProcessingSoftwareRef =
      IdentificationDataInternal::ProcessingSoftwareRef;

    using ProcessingStep = IdentificationDataInternal::ProcessingStep;
    using ProcessingSteps = IdentificationDataInternal::ProcessingSteps;
    using ProcessingStepRef = IdentificationDataInternal::ProcessingStepRef;

    using DBSearchParam = IdentificationDataInternal::DBSearchParam;
    using DBSearchParams = IdentificationDataInternal::DBSearchParams;
    using SearchParamRef = IdentificationDataInternal::SearchParamRef;
    using DBSearchSteps = IdentificationDataInternal::DBSearchSteps;

    using ScoreType = IdentificationDataInternal::ScoreType;
    using ScoreTypes = IdentificationDataInternal::ScoreTypes;
    using ScoreTypeRef = IdentificationDataInternal::ScoreTypeRef;

    using ScoredProcessingResult =
      IdentificationDataInternal::ScoredProcessingResult;

    using AppliedProcessingStep =
      IdentificationDataInternal::AppliedProcessingStep;
    using AppliedProcessingSteps =
      IdentificationDataInternal::AppliedProcessingSteps;

    using Observation = IdentificationDataInternal::Observation;
    using Observations = IdentificationDataInternal::Observations;
    using ObservationRef = IdentificationDataInternal::ObservationRef;

    using ParentSequence = IdentificationDataInternal::ParentSequence;
    using ParentSequences = IdentificationDataInternal::ParentSequences;
    using ParentSequenceRef = IdentificationDataInternal::ParentSequenceRef;

    using ParentMatch = IdentificationDataInternal::ParentMatch;
    using ParentMatches = IdentificationDataInternal::ParentMatches;

    using IdentifiedPeptide = IdentificationDataInternal::IdentifiedPeptide;
    using IdentifiedPeptides = IdentificationDataInternal::IdentifiedPeptides;
    using IdentifiedPeptideRef =
      IdentificationDataInternal::IdentifiedPeptideRef;

    using IdentifiedCompound = IdentificationDataInternal::IdentifiedCompound;
    using IdentifiedCompounds = IdentificationDataInternal::IdentifiedCompounds;
    using IdentifiedCompoundRef =
      IdentificationDataInternal::IdentifiedCompoundRef;

    using IdentifiedOligo = IdentificationDataInternal::IdentifiedOligo;
    using IdentifiedOligos = IdentificationDataInternal::IdentifiedOligos;
    using IdentifiedOligoRef = IdentificationDataInternal::IdentifiedOligoRef;

    using IdentifiedMolecule = IdentificationDataInternal::IdentifiedMolecule;

    using PeakAnnotations = IdentificationDataInternal::PeakAnnotations;

    using Adducts = IdentificationDataInternal::Adducts;
    using AdductRef = IdentificationDataInternal::AdductRef;
    using AdductOpt = IdentificationDataInternal::AdductOpt;

    using ObservationMatch = IdentificationDataInternal::ObservationMatch;
    using ObservationMatches = IdentificationDataInternal::ObservationMatches;
    using ObservationMatchRef = IdentificationDataInternal::ObservationMatchRef;

    // @todo: allow multiple sets of groups, like with parent sequences
    // ("ParentGroupSets")?
    using ObservationMatchGroup = IdentificationDataInternal::ObservationMatchGroup;
    using ObservationMatchGroups = IdentificationDataInternal::ObservationMatchGroups;
    using MatchGroupRef = IdentificationDataInternal::MatchGroupRef;

    using ParentGroup = IdentificationDataInternal::ParentGroup;
    using ParentGroups =
      IdentificationDataInternal::ParentGroups;
    using ParentGroupRef = IdentificationDataInternal::ParentGroupRef;
    using ParentGroupSet =
      IdentificationDataInternal::ParentGroupSet;
    using ParentGroupSets =
      IdentificationDataInternal::ParentGroupSets;

    using AddressLookup = std::unordered_set<uintptr_t>;

    /// structure that maps references of corresponding objects after copying
    struct RefTranslator {
      std::map<InputFileRef, InputFileRef> input_file_refs;
      std::map<ScoreTypeRef, ScoreTypeRef> score_type_refs;
      std::map<ProcessingSoftwareRef, ProcessingSoftwareRef> processing_software_refs;
      std::map<SearchParamRef, SearchParamRef> search_param_refs;
      std::map<ProcessingStepRef, ProcessingStepRef> processing_step_refs;
      std::map<ObservationRef, ObservationRef> observation_refs;
      std::map<ParentSequenceRef, ParentSequenceRef> parent_sequence_refs;
      std::map<IdentifiedPeptideRef, IdentifiedPeptideRef> identified_peptide_refs;
      std::map<IdentifiedOligoRef, IdentifiedOligoRef> identified_oligo_refs;
      std::map<IdentifiedCompoundRef, IdentifiedCompoundRef> identified_compound_refs;
      std::map<AdductRef, AdductRef> adduct_refs;
      std::map<ObservationMatchRef, ObservationMatchRef> observation_match_refs;

      bool allow_missing = false;

      IdentifiedMolecule translate(IdentifiedMolecule old) const;

      ObservationMatchRef translate(ObservationMatchRef old) const;

    };

    /// Default constructor
    IdentificationData():
      current_step_ref_(processing_steps_.end()), no_checks_(false)
    {
    }

    /*!
      @brief Copy constructor

      Copy-constructing is expensive due to the necessary "rewiring" of references.
      Use the move constructor where possible.
    */
    IdentificationData(const IdentificationData& other);

    /*!
      @brief Copy assignment operator
    */
    IdentificationData& operator=(const IdentificationData& other);

    /*!
      @brief Move constructor
    */
    IdentificationData(IdentificationData&& other) noexcept;

    /*!
      @brief Move assignment operator
    */
    IdentificationData& operator=(IdentificationData&& other) noexcept;

    /*!
      @brief Register an input file

      @return Reference to the registered file
    */
    InputFileRef registerInputFile(const InputFile& file);

    /*!
      @brief Register data processing software

      @return Reference to the registered software
    */
    ProcessingSoftwareRef registerProcessingSoftware(
      const ProcessingSoftware& software);

    /*!
      @brief Register database search parameters

      @return Reference to the registered search parameters
    */
    SearchParamRef registerDBSearchParam(const DBSearchParam& param);

    /*!
      @brief Register a data processing step

      @return Reference to the registered processing step
    */
    ProcessingStepRef registerProcessingStep(const ProcessingStep&
                                                 step);

    /*!
      @brief Register a database search step with associated parameters

      @return Reference to the registered processing step
    */
    ProcessingStepRef registerProcessingStep(
      const ProcessingStep& step, SearchParamRef search_ref);

    /*!
      @brief Register a score type

      @return Reference to the registered score type
    */
    ScoreTypeRef registerScoreType(const ScoreType& score);

    /*!
      @brief Register an observation (e.g. MS2 spectrum or feature)

      @return Reference to the registered observation
    */
    ObservationRef registerObservation(const Observation& obs);

    /*!
      @brief Register a parent sequence (e.g. protein or intact RNA)

      @return Reference to the registered parent sequence
    */
    ParentSequenceRef registerParentSequence(const ParentSequence& parent);

    /// Register a grouping of parent sequences (e.g. protein inference result)
    void registerParentGroupSet(const ParentGroupSet& groups);

    /*!
      @brief Register an identified peptide

      @return Reference to the registered peptide
    */
    IdentifiedPeptideRef registerIdentifiedPeptide(const IdentifiedPeptide&
                                                   peptide);

    /*!
      @brief Register an identified compound (small molecule)

      @return Reference to the registered compound
    */
    IdentifiedCompoundRef registerIdentifiedCompound(const IdentifiedCompound&
                                                     compound);

    /*!
      @brief Register an identified RNA oligonucleotide

      @return Reference to the registered oligonucleotide
    */
    IdentifiedOligoRef registerIdentifiedOligo(const IdentifiedOligo& oligo);

    /*!
      @brief Register an adduct

      @return Reference to the registered adduct
    */
    AdductRef registerAdduct(const AdductInfo& adduct);

    /*!
      @brief Register an observation match (e.g. peptide-spectrum match)

      @return Reference to the registered observation match
    */
    ObservationMatchRef registerObservationMatch(const ObservationMatch& match);

    /*!
      @brief Register a group of observation matches that belong together

      @return Reference to the registered group of observation matches
    */
    MatchGroupRef registerObservationMatchGroup(const ObservationMatchGroup& group);

    /// Return the registered input files (immutable)
    const InputFiles& getInputFiles() const
    {
      return input_files_;
    }

    /// Return the registered data processing software (immutable)
    const ProcessingSoftwares& getProcessingSoftwares() const
    {
      return processing_softwares_;
    }

    /// Return the registered data processing steps (immutable)
    const ProcessingSteps& getProcessingSteps() const
    {
      return processing_steps_;
    }

    /// Return the registered database search parameters (immutable)
    const DBSearchParams& getDBSearchParams() const
    {
      return db_search_params_;
    }

    /// Return the registered database search steps (immutable)
    const DBSearchSteps& getDBSearchSteps() const
    {
      return db_search_steps_;
    }

    /// Return the registered score types (immutable)
    const ScoreTypes& getScoreTypes() const
    {
      return score_types_;
    }

    /// Return the registered observations (immutable)
    const Observations& getObservations() const
    {
      return observations_;
    }

    /// Return the registered parent sequences (immutable)
    const ParentSequences& getParentSequences() const
    {
      return parents_;
    }

    /// Return the registered parent sequence groupings (immutable)
    const ParentGroupSets& getParentGroupSets() const
    {
      return parent_groups_;
    }

    /// Return the registered identified peptides (immutable)
    const IdentifiedPeptides& getIdentifiedPeptides() const
    {
      return identified_peptides_;
    }

    /// Return the registered compounds (immutable)
    const IdentifiedCompounds& getIdentifiedCompounds() const
    {
      return identified_compounds_;
    }

    /// Return the registered identified oligonucleotides (immutable)
    const IdentifiedOligos& getIdentifiedOligos() const
    {
      return identified_oligos_;
    }

    /// Return the registered adducts (immutable)
    const Adducts& getAdducts() const
    {
      return adducts_;
    }

    /// Return the registered observation matches (immutable)
    const ObservationMatches& getObservationMatches() const
    {
      return observation_matches_;
    }

    /// Return the registered groups of observation matches (immutable)
    const ObservationMatchGroups& getObservationMatchGroups() const
    {
      return observation_match_groups_;
    }

    /// Add a score to an input match (e.g. PSM)
    void addScore(ObservationMatchRef match_ref, ScoreTypeRef score_ref,
                  double value);

    /*!
      @brief Set a data processing step that will apply to all subsequent "register..." calls.

      This step will be appended to the list of processing steps for all relevant elements that are registered subsequently (unless it is already the last entry in the list).
      If a score type without a software reference is registered, the software reference of this processing step will be applied.
      Effective until @ref clearCurrentProcessingStep() is called.
    */
    void setCurrentProcessingStep(ProcessingStepRef step_ref);

    /*!
      @brief Return the current processing step (set via @ref setCurrentProcessingStep()).

      If no current processing step has been set, @p processing_steps.end() is returned.
    */
    ProcessingStepRef getCurrentProcessingStep();

    /// Cancel the effect of @ref setCurrentProcessingStep().
    void clearCurrentProcessingStep();

    /*!
      @brief Return the best match for each observation, according to a given score type

      @param score_ref Score type to use
      @param require_score Exclude matches without score of this type, even if they are the only matches for their observations?
    */
    std::vector<ObservationMatchRef> getBestMatchPerObservation(ScoreTypeRef score_ref,
                                                                bool require_score = false) const;
    // @todo: this currently doesn't take molecule type into account - should it?

    /// Get range of matches (cf. @p equal_range) for a given observation
    std::pair<ObservationMatchRef, ObservationMatchRef> getMatchesForObservation(ObservationRef obs_ref) const;

    /*!
      @brief Helper function for filtering observation matches (e.g. PSMs) in IdentificationData

      If other parts are invalidated by filtering, the data structure is automatically cleaned up (IdentificationData::cleanup) to remove any invalidated references at the end of this operation.

      @param func Functor that returns true for container elements to be removed
    */
    template <typename PredicateType>
    void removeObservationMatchesIf(PredicateType&& func)
    {
      auto count = observation_matches_.size();
      removeFromSetIf_(observation_matches_, func);
      if (count != observation_matches_.size()) cleanup();
    }

    /*!
      @brief Helper function for filtering parent sequences (e.g. protein sequences) in IdentificationData

      If other parts are invalidated by filtering, the data structure is automatically cleaned up (IdentificationData::cleanup) to remove any invalidated references at the end of this operation.

      @param func Functor that returns true for container elements to be removed
    */
    template <typename PredicateType>
    void removeParentSequencesIf(PredicateType&& func)
    {
      auto count = parents_.size();
      removeFromSetIf_(parents_, func);
      if (count != parents_.size()) cleanup();
    }

    template <typename PredicateType>
    void applyToObservations(PredicateType&& func)
    {
      for (auto it = observations_.begin(); it != observations_.end(); ++it)
        observations_.modify(it, func);
    }

    /*!
      @brief Look up a score type by name.

      @return Reference to the score type, if found; otherwise @p getScoreTypes().end()
    */
    ScoreTypeRef findScoreType(const String& score_name) const;

    /// Calculate sequence coverages of parent sequences
    void calculateCoverages(bool check_molecule_length = false);

    /*!
      @brief Clean up the data structure after filtering parts of it.

      Make sure there are no invalid references or "orphan" data entries.

      @param require_observation_match Remove identified molecules, observations and adducts that aren't part of observation matches?
      @param require_identified_sequence Remove parent sequences (proteins/RNAs) that aren't referenced by identified peptides/oligonucleotides?
      @param require_parent_match Remove identified peptides/oligonucleotides that don't reference a parent sequence (protein/RNA)?
      @param require_parent_group Remove parent sequences that aren't part of parent sequence groups?
      @param require_match_group Remove input matches that aren't part of match groups?
    */
    void cleanup(bool require_observation_match = true,
                 bool require_identified_sequence = true,
                 bool require_parent_match = true,
                 bool require_parent_group = false,
                 bool require_match_group = false);

    /// Return whether the data structure is empty (no data)
    bool empty() const;

    /*!
      @brief Merge in data from another instance.

      Can be used to make a deep copy by calling merge() on an empty object.
      The returned translation table allows updating of references that are held externally.

      @param other Instance to merge in.

      @return Translation table for references (old -> new)
    */
    RefTranslator merge(const IdentificationData& other);

    /// Swap contents with a second instance
    void swap(IdentificationData& other);

    /// Clear all contents
    void clear();

    /*!
      Pick a score type for operations (e.g. filtering) on a container of scored processing results (e.g. input matches, identified peptides, ...).

      If @p all_elements is false, only the first element with a score will be considered (which is sufficient if all elements were processed in the same way).
      If @p all_elements is true, the score type supported by the highest number of elements will be chosen.

      If @p any_score is false, only the primary score from the most recent processing step (that assigned a score) is taken into account.
      If @p any_score is true, all score types assigned across all elements are considered (this implies @p all_elements = true).

      @param container Container with elements derived from @p ScoredProcessingResult
      @param all_elements Consider all elements?
      @param any_score Consider any score (or just primary/most recent ones)?

      @return Reference to the chosen score type (or @p getScoreTypes().end() if there were no scores)
    */
    template <class ScoredProcessingResults>
    ScoreTypeRef pickScoreType(const ScoredProcessingResults& container,
                               bool all_elements = false, bool any_score = false) const
    {
      std::map<ScoreTypeRef, Size> score_counts;

      if (any_score)
      {
        for (const auto& element : container)
        {
          for (const auto& step : element.steps_and_scores)
          {
            for (const auto& pair : step.scores)
            {
              score_counts[pair.first]++;
            }
          }
        }
      }
      else
      {
        for (const auto& element : container)
        {
          auto score_info = element.getMostRecentScore();
          if (std::get<2>(score_info)) // check success indicator
          {
            ScoreTypeRef score_ref = *std::get<1>(score_info); // unpack the option
            if (!all_elements) return score_ref;
            score_counts[score_ref]++; // elements are zero-initialized
          }
        }
      }
      if (score_counts.empty()) return score_types_.end();
      auto pos = max_element(score_counts.begin(), score_counts.end());
      // @TODO: break ties according to some criterion
      return pos->first;
    }

    /// Set a meta value on a stored observation match (e.g. PSM)
    void setMetaValue(const ObservationMatchRef ref, const String& key, const DataValue& value);

    /// Set a meta value on a stored observation
    void setMetaValue(const ObservationRef ref, const String& key, const DataValue& value);

    /// Set a meta value on a stored identified molecule (variant)
    void setMetaValue(const IdentifiedMolecule& var, const String& key, const DataValue& value);

    // @TODO: add overloads for other data types derived from MetaInfoInterface

    /// Remove a meta value (if it exists) from a stored observation match (e.g. PSM)
    /// @todo: return whether value existed? (requires changes in MetaInfo[Interface])
    void removeMetaValue(const ObservationMatchRef ref, const String& key);

  protected:

    // containers:
    InputFiles input_files_;
    ProcessingSoftwares processing_softwares_;
    ProcessingSteps processing_steps_;
    DBSearchParams db_search_params_;
    // @TODO: store SearchParamRef inside ProcessingStep? (may not be required
    // for many processing steps)
    DBSearchSteps db_search_steps_;
    ScoreTypes score_types_;
    Observations observations_;
    ParentSequences parents_;
    ParentGroupSets parent_groups_;
    IdentifiedPeptides identified_peptides_;
    IdentifiedCompounds identified_compounds_;
    IdentifiedOligos identified_oligos_;
    Adducts adducts_;
    ObservationMatches observation_matches_;
    ObservationMatchGroups observation_match_groups_;

    /// Reference to the current data processing step (see @ref setCurrentProcessingStep())
    ProcessingStepRef current_step_ref_;

    /*!
      @brief Suppress validity checks in @p register... calls?

      This is useful in situations where validity is already guaranteed (e.g. copying).
    */
    bool no_checks_;

    // look-up tables for fast checking of reference validity:
    AddressLookup observation_lookup_;
    AddressLookup parent_lookup_;
    // @TODO: just use one "identified_molecule_lookup_" for all molecule types?
    AddressLookup identified_peptide_lookup_;
    AddressLookup identified_compound_lookup_;
    AddressLookup identified_oligo_lookup_;
    AddressLookup observation_match_lookup_;

    /// Helper function to check if all score types are valid
    void checkScoreTypes_(const std::map<ScoreTypeRef, double>& scores) const;

    /// Helper function to check if all applied processing steps are valid
    void checkAppliedProcessingSteps_(const AppliedProcessingSteps&
                                      steps_and_scores) const;

    /// Helper function to check if all parent matches are valid
    void checkParentMatches_(const ParentMatches& matches,
                             MoleculeType expected_type) const;

    /*!
      @brief Helper function to merge scored processing results while updating references (to processing steps and score types)

      @param result Instance that gets updated
      @param other Instance to merge into @p result
      @param trans Mapping of corresponding references between @p other and @p result
    */
    void mergeScoredProcessingResults_(ScoredProcessingResult& result,
                                       const ScoredProcessingResult& other,
                                       const RefTranslator& trans);

    /*!
      @brief Helper functor for adding processing steps to elements in a @p boost::multi_index_container structure

      The validity of the processing step reference cannot be checked here!
    */
    template <typename ElementType>
    struct ModifyMultiIndexAddProcessingStep;

    /**
      @brief Helper functor for adding scores to elements in a @em boost::multi_index_container structure

      The validity of the score type reference cannot be checked here!
    */
    template <typename ElementType>
    struct ModifyMultiIndexAddScore;

    /**
      @brief Helper functor for removing invalid parent matches from elements in a @em boost::multi_index_container structure

      Used during filtering, to update parent matches after parents have been removed.
    */
    template <typename ElementType>
    struct ModifyMultiIndexRemoveParentMatches;

    /// Helper function for adding entries (derived from ScoredProcessingResult) to a @em boost::multi_index_container structure
    template <typename ContainerType, typename ElementType>
    typename ContainerType::iterator insertIntoMultiIndex_(ContainerType& container, const ElementType& element);

    /// Variant of insertIntoMultiIndex_() that also updates a look-up table of valid references (addresses)
    template <typename ContainerType, typename ElementType>
    typename ContainerType::iterator insertIntoMultiIndex_(
      ContainerType& container, const ElementType& element,
      AddressLookup& lookup);

  };

}
