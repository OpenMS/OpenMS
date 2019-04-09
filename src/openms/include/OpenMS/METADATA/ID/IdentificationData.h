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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/DataProcessingStep.h>
#include <OpenMS/METADATA/ID/DataQuery.h>
#include <OpenMS/METADATA/ID/DBSearchParam.h>
#include <OpenMS/METADATA/ID/IdentifiedCompound.h>
#include <OpenMS/METADATA/ID/IdentifiedSequence.h>
#include <OpenMS/METADATA/ID/MetaData.h>
#include <OpenMS/METADATA/ID/MoleculeParentMatch.h>
#include <OpenMS/METADATA/ID/MoleculeQueryMatch.h>
#include <OpenMS/METADATA/ID/ParentMolecule.h>
#include <OpenMS/METADATA/ID/ParentMoleculeGroup.h>
#include <OpenMS/METADATA/ID/QueryMatchGroup.h>
#include <OpenMS/METADATA/ID/ScoreType.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <boost/unordered_set.hpp>

namespace OpenMS
{
  /*!
    @brief Representation of spectrum identification results and associated data.

    This class provides capabilities for storing spectrum identification results from different types of experiments/molecules (proteomics: peptides/proteins, metabolomics: small molecules, "nucleomics": RNA).
    The class design has the following goals:
    - Provide one structure for storing all relevant data for spectrum identification results.
    - Store data non-redundantly.
    - Ensure consistency (e.g. no conflicting information; no "dangling references").
    - Allow convenient and efficient querying.
    - Support different types of experiments, as mentioned above, in one common framework.

    The following important subordinate classes are provided to represent different types of data:
    <table>
    <tr><th>Class <th>Represents <th>Key <th>Proteomics example <th>Corresponding legacy class
    <tr><td>DataProcessingStep <td>Information about a data processing step that was applied (e.g. input files, software used, parameters) <td>Combined information <td>Mascot search <td>ProteinIdentification
    <tr><td>DataQuery <td>A search query (with identifier, RT, m/z), i.e. an MS2 spectrum or feature (for accurate mass search) <td>Identifier <td>MS2 spectrum <td>PeptideIdentification
    <tr><td>ParentMolecule <td>An entry in a FASTA file with associated information (sequence, coverage, etc.) <td>Accession <td>Protein <td>ProteinHit
    <tr><td>IdentifiedPeptide/-Oligo/-Compound <td>An identified molecule of the respective type <td>Sequence (or identifier for a compound) <td>Peptide <td>PeptideHit
    <tr><td>MoleculeQueryMatch <td>A match between a query (DataQuery) and identified molecule (Identified...) <td>Combination of query and molecule references <td>Peptide-spectrum match (PSM) <td>PeptideIdentification/PeptideHit
    </table>

    To populate an IdentificationData instance with data, "register..." functions are used.
    These functions return "references" (implemented as iterators) that can be used to refer to stored data items and thus form connections.
    For example, a protein can be stored using registerParentMolecule, which returns a corresponding reference.
    This reference can be used to build an IdentifiedPeptide object that references the protein.
    An identified peptide referencing a protein can only be registered if that protein has been registered already, to ensure data consistency.
    Given the identified peptide, information about the associated protein can be retrieved efficiently by simply dereferencing the reference.

    To ensure non-redundancy, many data types have a "key" (see table above) to which a uniqueness constraint applies.
    This means only one item of such a type with a given key can be stored in an IdentificationData object.
    If items with an existing key are registered subsequently, attempts are made to merge new information (e.g. additional scores) into the existing entry.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI IdentificationData: public MetaInfoInterface
  {
  public:

    // type definitions:
    using MoleculeType = IdentificationDataInternal::MoleculeType;
    using MassType = IdentificationDataInternal::MassType;

    using InputFiles = IdentificationDataInternal::InputFiles;
    using InputFileRef = IdentificationDataInternal::InputFileRef;

    using DataProcessingSoftware =
      IdentificationDataInternal::DataProcessingSoftware;
    using DataProcessingSoftwares =
      IdentificationDataInternal::DataProcessingSoftwares;
    using ProcessingSoftwareRef =
      IdentificationDataInternal::ProcessingSoftwareRef;

    using DataProcessingStep = IdentificationDataInternal::DataProcessingStep;
    using DataProcessingSteps = IdentificationDataInternal::DataProcessingSteps;
    using ProcessingStepRef = IdentificationDataInternal::ProcessingStepRef;

    using DBSearchParam = IdentificationDataInternal::DBSearchParam;
    using DBSearchParams = IdentificationDataInternal::DBSearchParams;
    using SearchParamRef = IdentificationDataInternal::SearchParamRef;
    using DBSearchSteps = IdentificationDataInternal::DBSearchSteps;

    using ScoreType = IdentificationDataInternal::ScoreType;
    using ScoreTypes = IdentificationDataInternal::ScoreTypes;
    using ScoreTypeRef = IdentificationDataInternal::ScoreTypeRef;

    using AppliedProcessingStep =
      IdentificationDataInternal::AppliedProcessingStep;
    using AppliedProcessingSteps =
      IdentificationDataInternal::AppliedProcessingSteps;

    using DataQuery = IdentificationDataInternal::DataQuery;
    using DataQueries = IdentificationDataInternal::DataQueries;
    using DataQueryRef = IdentificationDataInternal::DataQueryRef;

    using ParentMolecule = IdentificationDataInternal::ParentMolecule;
    using ParentMolecules = IdentificationDataInternal::ParentMolecules;
    using ParentMoleculeRef = IdentificationDataInternal::ParentMoleculeRef;

    using MoleculeParentMatch = IdentificationDataInternal::MoleculeParentMatch;
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

    using PeakAnnotations = IdentificationDataInternal::PeakAnnotations;
    using IdentifiedMoleculeRef =
      IdentificationDataInternal::IdentifiedMoleculeRef;

    using MoleculeQueryMatch = IdentificationDataInternal::MoleculeQueryMatch;
    using MoleculeQueryMatches =
      IdentificationDataInternal::MoleculeQueryMatches;
    using QueryMatchRef = IdentificationDataInternal::QueryMatchRef;

    // @TODO: allow multiple sets of groups, like with parent molecules
    // ("ParentMoleculeGroupings")?
    using QueryMatchGroup = IdentificationDataInternal::QueryMatchGroup;
    using QueryMatchGroups = IdentificationDataInternal::QueryMatchGroups;
    using MatchGroupRef = IdentificationDataInternal::MatchGroupRef;

    using ParentMoleculeGroup = IdentificationDataInternal::ParentMoleculeGroup;
    using ParentMoleculeGroups =
      IdentificationDataInternal::ParentMoleculeGroups;
    using ParentGroupRef = IdentificationDataInternal::ParentGroupRef;
    using ParentMoleculeGrouping =
      IdentificationDataInternal::ParentMoleculeGrouping;
    using ParentMoleculeGroupings =
      IdentificationDataInternal::ParentMoleculeGroupings;

    using AddressLookup = boost::unordered_set<uintptr_t>;


    /// Default constructor
    IdentificationData():
      current_step_ref_(processing_steps_.end())
    {
    }

    // Copy constructor
    IdentificationData(const IdentificationData& other):
      MetaInfoInterface(other)
    {
      assign_(other);
    }

    // Assignment operator
    IdentificationData& operator=(const IdentificationData& other)     
    {
      MetaInfoInterface::operator=(other);
      assign_(other);
      return *this;
    }

    /// Move constructor
    IdentificationData(IdentificationData&& other):
      input_files_(std::move(other.input_files_)),
      processing_softwares_(std::move(other.processing_softwares_)),
      processing_steps_(std::move(other.processing_steps_)),
      db_search_params_(std::move(other.db_search_params_)),
      db_search_steps_(std::move(other.db_search_steps_)),
      score_types_(std::move(other.score_types_)),
      data_queries_(std::move(other.data_queries_)),
      parent_molecules_(std::move(other.parent_molecules_)),
      parent_molecule_groupings_(std::move(other.parent_molecule_groupings_)),
      identified_peptides_(std::move(other.identified_peptides_)),
      identified_compounds_(std::move(other.identified_compounds_)),
      identified_oligos_(std::move(other.identified_oligos_)),
      query_matches_(std::move(other.query_matches_)),
      query_match_groups_(std::move(other.query_match_groups_)),
      current_step_ref_(std::move(other.current_step_ref_)),
      // look-up tables:
      data_query_lookup_(std::move(other.data_query_lookup_)),
      parent_molecule_lookup_(std::move(other.parent_molecule_lookup_)),
      identified_peptide_lookup_(std::move(other.identified_peptide_lookup_)),
      identified_compound_lookup_(std::move(other.identified_compound_lookup_)),
      identified_oligo_lookup_(std::move(other.identified_oligo_lookup_)),
      query_match_lookup_(std::move(other.query_match_lookup_))
    {
    }

    /*!
      @brief Register an input file
      @return Reference to the registered file
    */
    InputFileRef registerInputFile(const String& file);

    /*!
      @brief Register data processing software
      @return Reference to the registered software
    */
    ProcessingSoftwareRef registerDataProcessingSoftware(
      const DataProcessingSoftware& software);

    /*!
      @brief Register database search parameters
      @return Reference to the registered search parameters
    */
    SearchParamRef registerDBSearchParam(const DBSearchParam& param);

    /*!
      @brief Register a data processing step
      @return Reference to the registered processing step
    */
    ProcessingStepRef registerDataProcessingStep(const DataProcessingStep&
                                                 step);

    /*!
      @brief Register a database search step with associated parameters
      @return Reference to the registered processing step
    */
    ProcessingStepRef registerDataProcessingStep(
      const DataProcessingStep& step, SearchParamRef search_ref);

    /*!
      @brief Register a score type
      @return Reference to the registered score type
    */
    ScoreTypeRef registerScoreType(const ScoreType& score);

    /*!
      @brief Register a data query (e.g. MS2 spectrum or feature)
      @return Reference to the registered data query
    */
    DataQueryRef registerDataQuery(const DataQuery& query);

    /*!
      @brief Register a parent molecule (e.g. protein or intact RNA)
      @return Reference to the registered parent molecule
    */
    ParentMoleculeRef registerParentMolecule(const ParentMolecule& parent);

    /// Register a grouping of parent molecules (e.g. protein inference result)
    void registerParentMoleculeGrouping(const ParentMoleculeGrouping& grouping);

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
      @brief Register a molecule-query match (e.g. peptide-spectrum match)
      @return Reference to the registered molecule-query match
    */
    QueryMatchRef registerMoleculeQueryMatch(const MoleculeQueryMatch& match);

    /*!
      @brief Register a group of associated molecule-query matches
      @return Reference to the registered group of matches
    */
    MatchGroupRef registerQueryMatchGroup(const QueryMatchGroup& group);

    /// Return the registered input files (immutable)
    const InputFiles& getInputFiles() const
    {
      return input_files_;
    }

    /// Return the registered data processing software (immutable)
    const DataProcessingSoftwares& getDataProcessingSoftwares() const
    {
      return processing_softwares_;
    }

    /// Return the registered data processing steps (immutable)
    const DataProcessingSteps& getDataProcessingSteps() const
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

    /// Return the registered data queries (immutable)
    const DataQueries& getDataQueries() const
    {
      return data_queries_;
    }

    /// Return the registered parent molecules (immutable)
    const ParentMolecules& getParentMolecules() const
    {
      return parent_molecules_;
    }

    /// Return the registered parent molecule groupings (immutable)
    const ParentMoleculeGroupings& getParentMoleculeGroupings() const
    {
      return parent_molecule_groupings_;
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

    /// Return the registered molecule-query matches (immutable)
    const MoleculeQueryMatches& getMoleculeQueryMatches() const
    {
      return query_matches_;
    }

    /// Return the registered groups of molecule-query matches (immutable)
    const QueryMatchGroups& getQueryMatchGroups() const
    {
      return query_match_groups_;
    }

    /// Add a score to a molecule-query match (e.g. PSM)
    void addScore(QueryMatchRef match_ref, ScoreTypeRef score_ref,
                  double value);

    /*!
      @brief Set a data processing step that will apply to all subsequent "register..." calls.

      This step will be appended to the list of processing steps for all relevant elements that are registered subsequently (unless it is already the last entry in the list).
      If a score type without a software reference is registered, the software reference of this processing step will be applied.

      Effective until @ref clearCurrentProcessingStep() is called.
     */
    void setCurrentProcessingStep(ProcessingStepRef step_ref);

    /*!
      Return the current processing step (set via @ref setCurrentProcessingStep()).

      If no current processing step has been set, @p processing_steps.end() is returned.
    */
    ProcessingStepRef getCurrentProcessingStep();

    /// Cancel the effect of @ref setCurrentProcessingStep().
    void clearCurrentProcessingStep();

    /// Return the best match for each data query, according to a given score type
    // @TODO: this currently doesn't take molecule type into account - shoud it?
    std::vector<QueryMatchRef> getBestMatchPerQuery(ScoreTypeRef
                                                    score_ref) const;

    /*!
      @brief Look up a score type by name
      @return A pair: 1. Reference to the score type, if found; 2. Boolean indicating success or failure
    */
    std::pair<ScoreTypeRef, bool> findScoreType(const String& score_name) const;

    /// Calculate sequence coverages of parent molecules
    void calculateCoverages(bool check_molecule_length = false);

    /*!
      @brief Clean up the data structure after filtering parts of it

      Make sure there are no invalid references or "orphan" data entries.

      @param require_query_match Remove identified molecules and data queries that aren't part of molecule-query matches?
      @param require_identified_sequence Remove parent molecules (proteins/RNAs) that aren't referenced by identified peptides/oligonucleotides?
      @param require_parent_match Remove identified peptides/oligonucleotides that don't reference a parent molecule (protein/RNA)?
      @param require_parent_group Remove parent molecules that aren't part of parent molecule groups?
      @param require_match_group Remove molecule-query matches that aren't part of match groups?
    */
    void cleanup(bool require_query_match = true,
                 bool require_identified_sequence = true,
                 bool require_parent_match = true,
                 bool require_parent_group = false,
                 bool require_match_group = false);

    /// Helper function to compare two scores
    static bool isBetterScore(double first, double second, bool higher_better)
    {
      if (higher_better) return first > second;
      return first < second;
    }

  protected:
    void assign_(const IdentificationData& other)
    {
      input_files_ = other.input_files_;
      processing_softwares_ = other.processing_softwares_;
      processing_steps_ = other.processing_steps_;
      db_search_params_ = other.db_search_params_;

      // rebuild db_search_steps
      for (const auto& d : other.db_search_steps_) 
      {
        auto ps = processing_steps_.find(*d.first);
        auto sp = db_search_params_.find(*d.second);
        db_search_steps_.insert({ps, sp});
      }

      score_types_ = other.score_types_;
      data_queries_ = other.data_queries_;

      // copy and modify references (faster)
      parent_molecules_ = other.parent_molecules_;
      updateStepsAndScores_(parent_molecules_);

      // need to rebuild parent_molecule_groupings
      for (const auto grouping : other.parent_molecule_groupings_) 
      {
        ParentMoleculeGrouping new_grouping;
        new_grouping.label = grouping.label; 
        for (const auto& g : grouping.groups)
        {
          ParentMoleculeGroup new_pg;
          // update group member references
          for (ParentMoleculeRef ref : g.parent_molecule_refs)
          {
            auto key = parent_molecules_.find(ref->accession);
            new_pg.parent_molecule_refs.insert(key);
          }

          // update group score references
          for (const auto& pair : g.scores)
          {
            auto key = score_types_.find(*pair.first);
            new_pg.scores.insert({key, pair.second});
          }
          new_grouping.groups.insert(new_pg);
        }
        parent_molecule_groupings_.push_back(new_grouping);
      }

      identified_peptides_ = other.identified_peptides_;
      updateStepsAndScores_(identified_peptides_);
      updateParentMatches_(identified_peptides_);

      identified_compounds_ = other.identified_compounds_;
      updateStepsAndScores_(identified_compounds_);

      identified_oligos_ = other.identified_oligos_;
      updateStepsAndScores_(identified_oligos_);
      updateParentMatches_(identified_oligos_);

      // rebuild molecule query matches
      for (const auto& qm : other.query_matches_)
      {
        const MoleculeType& t = qm.getMoleculeType();
        IdentifiedMoleculeRef new_imr;
        switch(t)
        {
          case MoleculeType::PROTEIN:
          {
            auto value = boost::get<IdentifiedPeptideRef>(&qm.identified_molecule_ref);
            new_imr = identified_peptides_.find((*value)->sequence);
            break;
          }
          case MoleculeType::COMPOUND:
          {
            auto value = boost::get<IdentifiedCompoundRef>(&qm.identified_molecule_ref);
            new_imr = identified_compounds_.find((*value)->identifier);
            break;
          }
          case MoleculeType::RNA:
          {
            auto value = boost::get<IdentifiedOligoRef>(&qm.identified_molecule_ref);
            new_imr = identified_oligos_.find((*value)->sequence);
            break;
          }
          default:
            String msg = "Invalid molecule type.";
            String val(static_cast<Int>(t));
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg, val);
          break;
        }
        DataQueryRef new_dqr = data_queries_.find(*qm.data_query_ref);
        MoleculeQueryMatch new_mqm(new_imr, new_dqr);
        new_mqm.charge = qm.charge;
        for (const auto& pas : qm.peak_annotations)
        {
          if (pas.first) // has peak annotation with processing step information?
          {
            ProcessingStepRef key = processing_steps_.find(**pas.first);
            new_mqm.peak_annotations.insert(make_pair(key, pas.second));
          }
          else
          {
            new_mqm.peak_annotations[boost::none] = pas.second;
          }
        }
        new_mqm.steps_and_scores = qm.steps_and_scores; // references will be updated below        
        query_matches_.insert(new_mqm);
      }
      updateStepsAndScores_(query_matches_);

      for (auto const& qmg : other.query_match_groups_)
      {
        std::set<QueryMatchRef> new_qmr;
        for (auto const& qmr : qmg.query_match_refs)
        {
          DataQueryRef new_dqr = data_queries_.find(*qmr->data_query_ref);
          const MoleculeType& t = qmr->getMoleculeType();
          IdentifiedMoleculeRef new_imr;
          switch(t)
          {
            case MoleculeType::PROTEIN:
            {
              auto value = boost::get<IdentifiedPeptideRef>(&qmr->identified_molecule_ref);
              new_imr = identified_peptides_.find((*value)->sequence);
              break;
            }
            case MoleculeType::COMPOUND:
            {
              auto value = boost::get<IdentifiedCompoundRef>(&qmr->identified_molecule_ref);
              new_imr = identified_compounds_.find((*value)->identifier);
              break;
            }
            case MoleculeType::RNA:
            {
              auto value = boost::get<IdentifiedOligoRef>(&qmr->identified_molecule_ref);
              new_imr = identified_oligos_.find((*value)->sequence);
              break;
            }
            default:
              String msg = "Invalid molecule type.";
              String val(static_cast<Int>(t));
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg, val);
              break;
          }
          
          QueryMatchRef key = query_matches_.find(make_tuple(new_dqr, new_imr));
          new_qmr.insert(key);
        }
        QueryMatchGroup new_qmg;
        new_qmg.query_match_refs = new_qmr;
        new_qmg.steps_and_scores = qmg.steps_and_scores; // updated below
        query_match_groups_.insert(new_qmg);  
      }
      updateStepsAndScores_(query_match_groups_);

      if (other.current_step_ref_ != other.processing_steps_.end())
      {
        current_step_ref_ = processing_steps_.find(*other.current_step_ref_);
      }
      else // set to our end iterator
      {
        current_step_ref_ = processing_steps_.end();
      }

      // look-up tables:
      updateAddressLookup_(data_queries_, data_query_lookup_);
      updateAddressLookup_(parent_molecules_, parent_molecule_lookup_);
      updateAddressLookup_(identified_peptides_, identified_peptide_lookup_);
      updateAddressLookup_(identified_compounds_, identified_compound_lookup_);
      updateAddressLookup_(identified_oligos_, identified_oligo_lookup_);
      updateAddressLookup_(query_matches_, query_match_lookup_);
    }


    // containers:
    InputFiles input_files_;
    DataProcessingSoftwares processing_softwares_;
    DataProcessingSteps processing_steps_;
    DBSearchParams db_search_params_;
    // @TODO: store SearchParamRef inside ProcessingStep? (may not be required
    // for many processing steps)
    DBSearchSteps db_search_steps_;
    ScoreTypes score_types_;
    DataQueries data_queries_;
    ParentMolecules parent_molecules_;
    ParentMoleculeGroupings parent_molecule_groupings_;
    IdentifiedPeptides identified_peptides_;
    IdentifiedCompounds identified_compounds_;
    IdentifiedOligos identified_oligos_;
    MoleculeQueryMatches query_matches_;
    QueryMatchGroups query_match_groups_;

    /// Reference to the current data processing step (see @ref setCurrentProcessingStep())
    ProcessingStepRef current_step_ref_;

    // look-up tables for fast checking of reference validity:
    AddressLookup data_query_lookup_;
    AddressLookup parent_molecule_lookup_;
    // @TODO: just use one "identified_molecule_lookup_" for all molecule types?
    AddressLookup identified_peptide_lookup_;
    AddressLookup identified_compound_lookup_;
    AddressLookup identified_oligo_lookup_;
    AddressLookup query_match_lookup_;

    /// Helper function to check if all score types are valid
    void checkScoreTypes_(const std::map<ScoreTypeRef, double>& scores) const;

    /// Helper function to check if all applied processing steps are valid
    void checkAppliedProcessingSteps_(const AppliedProcessingSteps&
                                      steps_and_scores) const;

    /// Helper function to check if all parent matches are valid
    void checkParentMatches_(const ParentMatches& matches,
                             MoleculeType expected_type) const;

    /*!
      @brief Helper functor for adding processing steps to elements in a @t boost::multi_index_container structure

      The validity of the processing step reference cannot be checked here!
    */
    template <typename ElementType>
    struct ModifyMultiIndexAddProcessingStep
    {
      ModifyMultiIndexAddProcessingStep(ProcessingStepRef step_ref):
        step_ref(step_ref)
      {
      }

      void operator()(ElementType& element)
      {
        element.addProcessingStep(step_ref);
      }

      ProcessingStepRef step_ref;
    };

    /*!
      @brief Helper functor for adding scores to elements in a @t boost::multi_index_container structure

      The validity of the score type reference cannot be checked here!
    */
    template <typename ElementType>
    struct ModifyMultiIndexAddScore
    {
      ModifyMultiIndexAddScore(ScoreTypeRef score_type_ref, double value):
        score_type_ref(score_type_ref), value(value)
      {
      }

      void operator()(ElementType& element)
      {
        if (element.steps_and_scores.empty())
        {
          element.addScore(score_type_ref, value);
        }
        else // add score to most recent step
        {
          element.addScore(score_type_ref, value,
                           element.steps_and_scores.back().processing_step_opt);
        }
      }

      ScoreTypeRef score_type_ref;
      double value;
    };

    /*!
      @brief Helper functor for removing invalid parent matches from elements in a @t boost::multi_index_container structure

      Used during filtering, to update parent matches after parents have been removed.
    */
    template <typename ElementType>
    struct ModifyMultiIndexRemoveParentMatches
    {
      ModifyMultiIndexRemoveParentMatches(const AddressLookup& lookup):
        lookup(lookup)
      {
      }

      void operator()(ElementType& element)
      {
        removeFromSetIf_(element.parent_matches,
                         [&](const ParentMatches::iterator it)
                         {
                           return !lookup.count(it->first);
                         });
      }

      const AddressLookup& lookup;
    };


    /// Helper function for updating steps and score entries and associated references
    template <typename ContainerType>
    void updateStepsAndScores_(ContainerType& target)
    {
      for (auto p = target.begin(); p != target.end(); ++p) 
      {
        AppliedProcessingSteps apss;
        for (const auto& s : p->steps_and_scores)
        {
          AppliedProcessingStep as;

          if (s.processing_step_opt) // boost::optional
          {
            auto ref = processing_steps_.find(**s.processing_step_opt);
            as.processing_step_opt = ref;
          }
          // update group score references
          for (const auto& pair : s.scores)
          {
            auto key = score_types_.find(*pair.first);
            as.scores.insert({key, pair.second});
          }
          apss.push_back(as);
        }
        auto modifier = [&apss](typename ContainerType::value_type& e) { e.steps_and_scores = apss; };
        target.modify(p, modifier);
      }
    }

    /// Helper function for updating references in parent matches
    template <typename ContainerType>
    void updateParentMatches_(ContainerType& target)
    {    
      for (auto p = target.begin(); p != target.end(); ++p) 
      {
        ParentMatches new_pms;

        for (const auto& m : p->parent_matches)
        {
          ParentMoleculeRef key = parent_molecules_.find(m.first->accession);
          new_pms[key] = m.second;
        }
        
        auto modifier = [&new_pms, &target](typename ContainerType::value_type& p) { p.parent_matches = new_pms; };
        target.modify(p, modifier);
      }
    }

    /// Helper function for adding entries (derived from ScoredProcessingResult) to a @t boost::multi_index_container structure
    template <typename ContainerType, typename ElementType>
    typename ContainerType::iterator insertIntoMultiIndex_(
      ContainerType& container, const ElementType& element)
    {
      checkAppliedProcessingSteps_(element.steps_and_scores);

      auto result = container.insert(element);
      if (!result.second) // existing element - merge in new information
      {
        container.modify(result.first, [&element](ElementType& existing)
                         {
                           existing += element;
                         });
      }

      // add current processing step (if necessary):
      if (current_step_ref_ != processing_steps_.end())
      {
        ModifyMultiIndexAddProcessingStep<ElementType>
          modifier(current_step_ref_);
        container.modify(result.first, modifier);
      }

      return result.first;
    }

    /// Variant of insertIntoMultiIndex_() that also updates a look-up table of valid references (addresses)
    template <typename ContainerType, typename ElementType>
    typename ContainerType::iterator insertIntoMultiIndex_(
      ContainerType& container, const ElementType& element,
      AddressLookup& lookup)
    {
      typename ContainerType::iterator ref =
        insertIntoMultiIndex_(container, element);
      lookup.insert(uintptr_t(&(*ref)));
      return ref;
    }

    /// Check whether a reference points to an element in a container
    template <typename RefType, typename ContainerType>
    static bool isValidReference_(RefType ref, ContainerType& container)
    {
      for (auto it = container.begin(); it != container.end(); ++it)
      {
        if (ref == it) return true;
      }
      return false;
    }

    /// Check validity of a reference based on a look-up table of addresses
    template <typename RefType>
    static bool isValidHashedReference_(
      RefType ref, const AddressLookup& lookup)
    {
      return lookup.count(ref);
    }

    /// Remove elements from a set (or ordered multi_index_container) if they fulfill a predicate
    template <typename ContainerType, typename PredicateType>
    // static void removeFromSetIf_(ContainerType& container, std::function<bool(RefType)> predicate)
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

    /// Remove elements from a set (or ordered multi_index_container) if they don't occur in a look-up table
    template <typename ContainerType>
    static void removeFromSetIfNotHashed_(
      ContainerType& container, const AddressLookup& lookup)
    {
      removeFromSetIf_(container, [&lookup](typename ContainerType::iterator it)
                       {
                         return !lookup.count(uintptr_t(&(*it)));
                       });
    }

    /// Recreate the address look-up table for a container
    template <typename ContainerType>
    static void updateAddressLookup_(const ContainerType& container,
                                     AddressLookup& lookup)
    {
      lookup.clear();
      lookup.reserve(container.size());
      for (const auto& element : container)
      {
        lookup.insert(uintptr_t(&element));
      }
    }


    // IDFilter needs access to do its job:
    friend class IDFilter;
  };
}
