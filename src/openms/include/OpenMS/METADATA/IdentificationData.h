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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDENTIFICATIONDATA_H
#define OPENMS_METADATA_IDENTIFICATIONDATA_H

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Software.h>

#include <boost/bimap.hpp>
#include <unordered_map>

namespace OpenMS
{
  class OPENMS_DLLAPI IdentificationData: public MetaInfoInterface
  {
  public:
    typedef UInt64 UniqueKey; // easy to change in case 64 bit isn't enough

    // Input files that were processed:
    typedef UniqueKey InputFileKey;
    typedef boost::bimap<InputFileKey, String> InputFileBimap;
    InputFileBimap input_files;


    /*!
      Information about software used for data processing.

      If the same processing is applied to multiple ID runs, e.g. if multiple files (fractions, replicates) are searched with the same search engine, store the
 software information only once.
    */
    struct DataProcessingSoftware
    {
      Software tool; // also captures CV terms and meta data (MetaInfoInterface)

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      explicit DataProcessingSoftware(
        const Software& tool = Software(),
        std::set<DataProcessing::ProcessingAction> actions =
        std::set<DataProcessing::ProcessingAction>()):
        tool(tool), actions(actions)
      {
      }

      explicit DataProcessingSoftware(
        const String& tool_name, const String& tool_version = "",
        std::set<DataProcessing::ProcessingAction> actions =
        std::set<DataProcessing::ProcessingAction>()):
        tool(), actions(actions)
      {
        tool.setName(tool_name);
        tool.setVersion(tool_version);
      }

      DataProcessingSoftware(const DataProcessingSoftware& other) = default;

      bool operator<(const DataProcessingSoftware& other) const
      {
        return (std::tie(tool.getName(), tool.getVersion(), actions) <
                std::tie(other.tool.getName(), other.tool.getVersion(),
                         other.actions));
      }

      bool operator==(const DataProcessingSoftware& other) const
      {
        return (std::tie(tool.getName(), tool.getVersion(), actions) ==
                std::tie(other.tool.getName(), other.tool.getVersion(),
                         other.actions));
      }
    };

    typedef UniqueKey ProcessingSoftwareKey;
    typedef boost::bimap<ProcessingSoftwareKey,
                         DataProcessingSoftware> SoftwareBimap;
    SoftwareBimap processing_software;


    /*!
      Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct DataProcessingStep: public MetaInfoInterface
    {
      ProcessingSoftwareKey software_key;

      std::vector<InputFileKey> input_files; // reference into "input_files"

      std::vector<String> primary_files; // path(s) to primary MS data

      DateTime date_time;

      explicit DataProcessingStep(
        ProcessingSoftwareKey software_key = 0,
        const std::vector<InputFileKey>& input_files =
        std::vector<InputFileKey>(), const std::vector<String>& primary_files =
        std::vector<String>(), const DateTime& date_time = DateTime::now()):
        software_key(software_key), input_files(input_files),
        primary_files(primary_files), date_time(date_time)
      {
      }

      DataProcessingStep(const DataProcessingStep& other) = default;

      bool operator<(const DataProcessingStep& other) const
      {
        return (std::tie(software_key, input_files, primary_files, date_time) <
                std::tie(other.software_key, other.input_files,
                         other.primary_files, other.date_time));
      }

      bool operator==(const DataProcessingStep& other) const
      {
        return (std::tie(software_key, input_files, primary_files, date_time) ==
                std::tie(other.software_key, other.input_files,
                         other.primary_files, other.date_time));
      }
    };

    typedef UniqueKey ProcessingStepKey;
    typedef boost::bimap<ProcessingStepKey, DataProcessingStep> StepsBimap;
    StepsBimap processing_steps;


    /*!
      Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term;

      String name;

      bool higher_better;

      // reference to the software that assigned the score:
      ProcessingSoftwareKey software_key;
      // @TODO: scores assigned by different software tools/versions are
      // considered as different scores - does that make sense?

      ScoreType():
        higher_better(true), software_key(0)
      {
      }

      explicit ScoreType(const CVTerm& cv_term, bool higher_better,
                         ProcessingSoftwareKey software_key = 0):
        cv_term(cv_term), name(cv_term.getName()), higher_better(higher_better),
        software_key(software_key)
      {
      }

      explicit ScoreType(const String& name, bool higher_better,
                         ProcessingSoftwareKey software_key = 0):
        cv_term(), name(name), higher_better(higher_better),
        software_key(software_key)
      {
      }

      ScoreType(const ScoreType& other) = default;

      // don't include "higher_better" in the comparison:
      bool operator<(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_key) <
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_key));
      }

      // don't include "higher_better" in the comparison:
      bool operator==(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_key) ==
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_key));
      }
    };

    typedef UniqueKey ScoreTypeKey;
    typedef boost::bimap<ScoreTypeKey, ScoreType> ScoreTypeBimap;
    typedef std::vector<std::pair<ScoreTypeKey, double>> ScoreList;
    ScoreTypeBimap score_types;


    /*!
      Search query, e.g. spectrum or feature.
    */
    struct DataQuery: public MetaInfoInterface
    {
      // spectrum or feature ID (from the file reference by "input_file_key"):
      String data_id;

      InputFileKey input_file_key; // reference into "input_files"

      double rt, mz; // position

      explicit DataQuery(const String& data_id = "",
                         InputFileKey input_file_key = 0,
                         double rt = std::numeric_limits<double>::quiet_NaN(),
                         double mz = std::numeric_limits<double>::quiet_NaN()):
        data_id(data_id), input_file_key(input_file_key), rt(rt), mz(mz)
      {
        // @TODO: require "input_file_key"? (see also "DataProcessingStep")
      }

      DataQuery(const DataQuery& other) = default;

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator<(const DataQuery& other) const
      {
        return std::tie(input_file_key, data_id) <
          std::tie(other.input_file_key, other.data_id);
      }

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator==(const DataQuery& other) const
      {
        return std::tie(input_file_key, data_id) ==
          std::tie(other.input_file_key, other.data_id);
      }

      // @TODO: do we need an "experiment label" (used e.g. in pepXML)?
      // if yes, should it be stored here or together with the input file?
    };

    typedef UniqueKey DataQueryKey;
    typedef boost::bimap<DataQueryKey, DataQuery> QueryBimap;
    QueryBimap data_queries;


    /// Base class for data with scores and processing steps (and meta info)
    struct ScoredProcessingResult: public MetaInfoInterface
    {
      ScoreList scores;

      std::vector<ProcessingStepKey> processing_steps;

      ScoredProcessingResult& operator+=(const ScoredProcessingResult& other)
      {
        // merge processing steps:
        for (auto step_key : other.processing_steps)
        {
          if (std::find(processing_steps.begin(), processing_steps.end(),
                        step_key) == processing_steps.end())
          {
            processing_steps.push_back(step_key);
          }
        }
        // merge scores:
        for (auto score_pair : other.scores)
        {
          // @TODO: should we overwrite scores?
          if (std::find(scores.begin(), scores.end(), score_pair) ==
              scores.end())
          {
            scores.push_back(score_pair);
          }
        }
        // merge meta info:
        std::vector<UInt> keys;
        other.getKeys(keys);
        for (const UInt key : keys)
        {
          // @TODO: should we overwrite meta values?
          if (!metaValueExists(key))
          {
            setMetaValue(key, other.getMetaValue(key));
          }
        }

        return *this;
      }

      std::pair<double, bool> getScore(ScoreTypeKey key) const
      {
        // give priority to "later" scores in the list:
        for (ScoreList::const_reverse_iterator it = scores.rbegin();
             it != scores.rend(); ++it)
        {
          if (it->first == key) return std::make_pair(it->second, true);
        }
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);
      }

    protected:
      explicit ScoredProcessingResult(
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepKey>& processing_steps =
        std::vector<ProcessingStepKey>()):
        scores(scores), processing_steps(processing_steps)
      {
      }

      ScoredProcessingResult(const ScoredProcessingResult& other) = default;
    };


    // Identified molecules - at the moment, peptides or small molecules:
    typedef UniqueKey IdentifiedMoleculeKey;
    typedef boost::bimap<IdentifiedMoleculeKey, AASequence> PeptideBimap;
    typedef boost::bimap<IdentifiedMoleculeKey, String> CompoundBimap;
    typedef boost::bimap<IdentifiedMoleculeKey, NASequence> OligoBimap;
    PeptideBimap identified_peptides;
    CompoundBimap identified_compounds;
    OligoBimap identified_oligos;

    enum MoleculeType
    {
      MT_PROTEIN,
      MT_COMPOUND,
      MT_RNA,
      SIZE_OF_MOLECULETYPES
    };

    /*!
      Meta data for an identified molecule.
    */
    struct IdentifiedMetaData: public ScoredProcessingResult
    {
      enum MoleculeType molecule_type; // @TODO: do we need this?

      explicit IdentifiedMetaData(
        enum MoleculeType molecule_type = MT_PROTEIN,
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepKey>& processing_steps =
        std::vector<ProcessingStepKey>()):
        ScoredProcessingResult(scores, processing_steps),
        molecule_type(molecule_type)
      {
      }

      IdentifiedMetaData(const IdentifiedMetaData& other) = default;
    };

    std::unordered_map<IdentifiedMoleculeKey,
                       IdentifiedMetaData> identified_meta_data;


    /*!
      Meta data specific to an identified compound (small molecule).
    */
    struct CompoundMetaData: public MetaInfoInterface
    {
      EmpiricalFormula formula;

      String name;

      String smile;

      String inchi;

      explicit CompoundMetaData(
        const EmpiricalFormula& formula = EmpiricalFormula(),
        const String& name = "", const String& smile = "",
        const String& inchi = ""):
        formula(formula), name(name), smile(smile), inchi(inchi)
      {
      }

      CompoundMetaData(const CompoundMetaData& other) = default;
    };

    std::unordered_map<IdentifiedMoleculeKey,
                       CompoundMetaData> compound_meta_data;


    /*!
      Meta data for a search hit (e.g. peptide-spectrum match).
    */

    // @TODO: move "PeakAnnotation" out of "PeptideHit"
    typedef std::vector<PeptideHit::PeakAnnotation> PeakAnnotations;

    struct MoleculeQueryMatch: public ScoredProcessingResult
    {
      Int charge;

      // peak annotations (fragment ion matches), potentially from different
      // data processing steps:
      std::map<ProcessingStepKey, PeakAnnotations> peak_annotations;

      explicit MoleculeQueryMatch(
        Int charge = 0, const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepKey>& processing_steps =
        std::vector<ProcessingStepKey>(),
        const std::map<ProcessingStepKey, PeakAnnotations>& peak_annotations =
        std::map<ProcessingStepKey, PeakAnnotations>()):
        ScoredProcessingResult(scores, processing_steps), charge(charge),
        peak_annotations(peak_annotations)
      {
      }

      MoleculeQueryMatch(const MoleculeQueryMatch& other) = default;
    };

    // all matches for the same data query should be consecutive!
    typedef std::pair<DataQueryKey, IdentifiedMoleculeKey> QueryMatchKey;
    typedef std::map<QueryMatchKey, MoleculeQueryMatch> QueryMatchMap;
    QueryMatchMap query_matches;


    /*!
      Representation of a parent molecule that is identified only indirectly (e.g. a protein).
    */
    struct ParentMetaData: public ScoredProcessingResult
    {
      enum MoleculeType molecule_type;

      String sequence;

      String description;

      double coverage;

      bool is_decoy;

      explicit ParentMetaData(
        enum MoleculeType molecule_type = MT_PROTEIN,
        const String& sequence = "", const String& description = "",
        double coverage = 0.0, bool is_decoy = false,
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepKey>& processing_steps =
        std::vector<ProcessingStepKey>()):
        ScoredProcessingResult(scores, processing_steps),
        molecule_type(molecule_type), sequence(sequence),
        description(description), coverage(coverage), is_decoy(is_decoy)
      {
      }

      ParentMetaData(const ParentMetaData& other) = default;
    };

    typedef UniqueKey ParentMoleculeKey;
    typedef boost::bimap<ParentMoleculeKey, String> ParentBimap;
    ParentBimap parent_molecules;
    std::unordered_map<ParentMoleculeKey, ParentMetaData> parent_meta_data;


    /*!
      Meta data for the association between an identified molecule (e.g. peptide) and a parent molecule (e.g. protein).
    */
    struct MoleculeParentMatch: public MetaInfoInterface
    {
      // in extraordinary cases (e.g. database searches that allow insertions/
      // deletions), the length of the identified molecule may differ from the
      // length of the subsequence in the parent; therefore, store "end_pos":
      Size start_pos, end_pos;

      char left_neighbor, right_neighbor; // neighboring sequence elements

      static const Size UNKNOWN_POSITION; // = Size(-1)
      static const char UNKNOWN_NEIGHBOR; // = 'X'
      static const char LEFT_TERMINUS; // = '['
      static const char RIGHT_TERMINUS; // = ']'

      explicit MoleculeParentMatch(Size start_pos = UNKNOWN_POSITION,
                                   Size end_pos = UNKNOWN_POSITION,
                                   char left_neighbor = UNKNOWN_NEIGHBOR,
                                   char right_neighbor = UNKNOWN_NEIGHBOR):
        start_pos(start_pos), end_pos(end_pos), left_neighbor(left_neighbor),
        right_neighbor(right_neighbor)
      {
      }

      bool operator<(const MoleculeParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) <
                std::tie(other.start_pos, other.end_pos));
      }

      bool operator==(const MoleculeParentMatch& other) const
      {
        // positions determine neighbors - no need to compare those:
        return (std::tie(start_pos, end_pos) ==
                std::tie(other.start_pos, other.end_pos));
      }

      bool hasValidPositions(Size molecule_length = 0, Size parent_length = 0)
      {
        if ((start_pos == UNKNOWN_POSITION) || (end_pos == UNKNOWN_POSITION))
        {
          return false;
        }
        if (end_pos < start_pos) return false;
        if (molecule_length && (end_pos - start_pos + 1 != molecule_length))
        {
          return false;
        }
        if (parent_length && (end_pos >= parent_length)) return false;
        return true;
      }
    };

    // mapping: identified molecule -> parent molecule -> match information
    typedef std::unordered_map<ParentMoleculeKey,
                               std::set<MoleculeParentMatch>> ParentSubMap;
    typedef std::unordered_map<IdentifiedMoleculeKey,
                               ParentSubMap> ParentMatchMap;
    ParentMatchMap parent_matches;


    /*!
      Parameters specific to a database search step.
    */
    struct DBSearchParameters: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;
      enum ProteinIdentification::PeakMassType peak_mass_type;

      String database;
      String database_version;
      String taxonomy;

      std::set<Int> charges;

      std::set<String> fixed_mods;
      std::set<String> variable_mods;

      double precursor_mass_tolerance;
      double fragment_mass_tolerance;
      bool precursor_tolerance_ppm;
      bool fragment_tolerance_ppm;

      // allow for either "DigestionEnzymeProtein" or "DigestionEnzymeRNA":
      const DigestionEnzyme* digestion_enzyme;
      Size missed_cleavages;
      Size min_length;
      Size max_length;

      DBSearchParameters():
        molecule_type(MT_PROTEIN),
        peak_mass_type(ProteinIdentification::MONOISOTOPIC),
        precursor_mass_tolerance(0.0), fragment_mass_tolerance(0.0),
        precursor_tolerance_ppm(false), fragment_tolerance_ppm(false),
        digestion_enzyme(0), missed_cleavages(0), min_length(0), max_length(0)
      {
      }

      DBSearchParameters(const DBSearchParameters& other) = default;

      bool operator<(const DBSearchParameters& other) const
      {
        return (std::tie(molecule_type, peak_mass_type, database,
                         database_version, taxonomy, charges, fixed_mods,
                         variable_mods, fragment_mass_tolerance,
                         precursor_mass_tolerance, fragment_tolerance_ppm,
                         precursor_tolerance_ppm, digestion_enzyme,
                         missed_cleavages, min_length, max_length) <
                std::tie(other.molecule_type, other.peak_mass_type,
                         other.database, other.database_version, other.taxonomy,
                         other.charges, other.fixed_mods, other.variable_mods,
                         other.fragment_mass_tolerance,
                         other.precursor_mass_tolerance,
                         other.fragment_tolerance_ppm,
                         other.precursor_tolerance_ppm,
                         other.digestion_enzyme, other.missed_cleavages,
                         other.min_length, other.max_length));
      }

      bool operator==(const DBSearchParameters& other) const
      {
        return (std::tie(molecule_type, peak_mass_type, database,
                         database_version, taxonomy, charges, fixed_mods,
                         variable_mods, fragment_mass_tolerance,
                         precursor_mass_tolerance, fragment_tolerance_ppm,
                         precursor_tolerance_ppm, digestion_enzyme,
                         missed_cleavages, min_length, max_length) ==
                std::tie(other.molecule_type, other.peak_mass_type,
                         other.database, other.database_version, other.taxonomy,
                         other.charges, other.fixed_mods, other.variable_mods,
                         other.fragment_mass_tolerance,
                         other.precursor_mass_tolerance,
                         other.fragment_tolerance_ppm,
                         other.precursor_tolerance_ppm,
                         other.digestion_enzyme, other.missed_cleavages,
                         other.min_length, other.max_length));
      }
    };

    typedef UniqueKey SearchParamsKey;
    typedef boost::bimap<SearchParamsKey, DBSearchParameters> SearchParamsBimap;
    SearchParamsBimap db_search_params;
    std::unordered_map<ProcessingStepKey, SearchParamsKey> db_search_steps;


    /// Default constructor
    IdentificationData():
      current_step_key_(0)
    {
    }

    /// Copy constructor
    IdentificationData(const IdentificationData& other) = default;

    /// Import from legacy peptide/protein identifications
    void importIDs(const std::vector<ProteinIdentification>& proteins,
                   const std::vector<PeptideIdentification>& peptides);

    /// Export to legacy peptide/protein identifications
    void exportIDs(std::vector<ProteinIdentification>& proteins,
                   std::vector<PeptideIdentification>& peptides) const;

    /// Export to mzTab format
    MzTab exportMzTab() const;

    /// Helper function to store meta data (derived from ScoredProcessingResult)
    template <typename MetaDataMap>
    void insertMetaData_(const typename MetaDataMap::mapped_type& meta_data,
                         MetaDataMap& meta_data_map,
                         typename MetaDataMap::key_type key, bool new_item)
    {
      checkScoreTypes_(meta_data.scores);
      checkProcessingSteps_(meta_data.processing_steps);

      if (new_item)
      {
        meta_data_map.insert(std::make_pair(key, meta_data));
      }
      else
      {
        meta_data_map.at(key) += meta_data;
      }
      addCurrentProcessingStep_(meta_data_map.at(key).processing_steps);
    }

    std::pair<InputFileKey, bool> registerInputFile(const String& file);

    std::pair<ProcessingSoftwareKey, bool> registerDataProcessingSoftware(
      const DataProcessingSoftware& software);

    std::pair<SearchParamsKey, bool> registerDBSearchParameters(
      const DBSearchParameters& params);

    std::pair<ProcessingStepKey, bool> registerDataProcessingStep(
      const DataProcessingStep& step, SearchParamsKey search_key = 0);

    std::pair<ScoreTypeKey, bool> registerScoreType(const ScoreType& score);

    std::pair<DataQueryKey, bool> registerDataQuery(const DataQuery& query);

    std::pair<IdentifiedMoleculeKey, bool> registerPeptide(
      const AASequence& seq,
      const IdentifiedMetaData& meta_data = IdentifiedMetaData());

    std::pair<IdentifiedMoleculeKey, bool> registerCompound(
      const String& id,
      const CompoundMetaData& compound_meta = CompoundMetaData(),
      const IdentifiedMetaData& id_meta = IdentifiedMetaData(MT_COMPOUND));

    std::pair<IdentifiedMoleculeKey, bool> registerOligo(
      const NASequence& seq,
      const IdentifiedMetaData& meta_data = IdentifiedMetaData(MT_RNA));

    std::pair<ParentMoleculeKey, bool> registerParentMolecule(
      const String& accession,
      const ParentMetaData& meta_data = ParentMetaData());

    // these ones are called "add..." instead of "register..." because they
    // don't return a key:

    bool addMoleculeParentMatch(
      IdentifiedMoleculeKey molecule_key, ParentMoleculeKey parent_key,
      const MoleculeParentMatch& meta_data = MoleculeParentMatch());

    bool addMoleculeQueryMatch(
      IdentifiedMoleculeKey molecule_key, DataQueryKey query_key,
      const MoleculeQueryMatch& meta_data = MoleculeQueryMatch(),
      const PeakAnnotations& peak_annotations = PeakAnnotations());

    /*!
      @brief Set a data processing step that will apply to all subsequent "register..." calls.

      This step will be appended to the list of processing steps for all relevant elements that are registered subsequently (unless it is already the last entry in the list).
      If a score type without a software reference is registered, the software reference of this processing step will be applied.

      Effective until @ref clearCurrentProcessingStep() is called.
     */
    void setCurrentProcessingStep(ProcessingStepKey step_key);

    /*!
      Cancel the effect of @ref setCurrentProcessingStep().
    */
    void clearCurrentProcessingStep();

    std::vector<QueryMatchKey> getBestMatchPerQuery(ScoreTypeKey score_key);

    ScoreTypeKey findScoreType(const String& score_name,
                               ProcessingSoftwareKey software_key = 0);

    bool allParentsAreDecoys(IdentifiedMoleculeKey molecule_key);

  protected:

    /// Key of the current data processing step (see @ref setCurrentProcessingStep())
    ProcessingStepKey current_step_key_;

    /*!
      Helper function for inserting an item into a bidirectional map.

      If the item is not contained yet, a unique key for it will be generated.

      Returns the key of the item and whether an insertion took place.
    */
    template<typename KeyType, typename ItemType>
    std::pair<KeyType, bool> insertIntoBimap_(
      const ItemType& item, boost::bimap<KeyType, ItemType>& container)
    {
      typename boost::bimap<KeyType, ItemType>::right_const_iterator pos =
        container.right.find(item);
      if (pos == container.right.end()) // new entry
      {
        KeyType key = KeyType(UniqueIdGenerator::getUniqueId());
        container.insert(typename boost::bimap<KeyType, ItemType>::
                         value_type(key, item));
        return std::make_pair(key, true);
      }
      return std::make_pair(pos->second, false);
    }

    /// Export a parent molecule (protein or nucleic acid) to mzTab
    template <typename MzTabSectionRow>
    void exportParentMoleculeToMzTab_(const String& accession,
                                      const ParentMetaData& meta_data,
                                      std::vector<MzTabSectionRow>& output,
                                      std::map<ScoreTypeKey, Size>& score_map)
      const
    {
      MzTabSectionRow parent;
      parent.accession.set(accession);
      exportScoresToMzTab_(meta_data.scores, parent.best_search_engine_score,
                           score_map);
      exportProcessingStepsToMzTab_(meta_data.processing_steps,
                                    parent.search_engine);
      parent.description.set(meta_data.description);
      parent.coverage.set(meta_data.coverage);
      if (!meta_data.sequence.empty())
      {
        MzTabOptionalColumnEntry opt_seq;
        opt_seq.first = "opt_sequence";
        opt_seq.second.set(meta_data.sequence);
        parent.opt_.push_back(opt_seq);
      }
      output.push_back(parent);
    }

    /// Export an identified molecule (peptide or oligonucleotide, but not small molecule/compound) to mzTab
    template <typename MzTabSectionRow>
    void exportPeptideOrOligoToMzTab_(const String& sequence,
                                      IdentifiedMoleculeKey key,
                                      std::vector<MzTabSectionRow>& output,
                                      std::map<ScoreTypeKey, Size>& score_map)
      const
    {
      MzTabSectionRow molecule;
      // @TODO: handle modifications properly
      molecule.sequence.set(sequence);
      const IdentifiedMetaData& meta_data = identified_meta_data.at(key);
      exportScoresToMzTab_(meta_data.scores, molecule.best_search_engine_score,
                           score_map);
      exportProcessingStepsToMzTab_(meta_data.processing_steps,
                                    molecule.search_engine);
      // generate one entry (with duplicated data) for every accession:
      ParentMatchMap::const_iterator pos = parent_matches.find(key);
      if (pos != parent_matches.end())
      {
        bool unique = (pos->second.size() == 1);
        for (const std::pair<ParentMoleculeKey, std::set<MoleculeParentMatch>>&
               match_pair : pos->second)
        {
          const String& accession = parent_molecules.left.at(match_pair.first);
          molecule.accession.set(accession);
          molecule.unique.set(unique);
          if (match_pair.second.empty())
          {
            output.push_back(molecule);
          }
          else
          {
            addMzTabMoleculeParentContext_(match_pair.second, molecule, output);
          }
        }
      }
      else // no accession given
      {
        // molecule.unique.set(false); // leave this unset?
        output.push_back(molecule);
      }
    }

    /// Export a molecule-query match (peptide- or oligonucleotide-spectrum match) to mzTab
    template <typename MzTabSectionRow>
    void exportQueryMatchToMzTab_(const String& sequence,
                                  const MoleculeQueryMatch& match,
                                  DataQueryKey query_key, double calc_mass,
                                  std::vector<MzTabSectionRow>& output,
                                  std::map<ScoreTypeKey, Size>& score_map,
                                  std::map<InputFileKey, Size>& file_map)
      const
    {
      MzTabSectionRow xsm; // PSM or OSM
      // @TODO: handle modifications properly
      xsm.sequence.set(sequence);
      exportScoresToMzTab_(match.scores, xsm.search_engine_score, score_map);
      exportProcessingStepsToMzTab_(match.processing_steps, xsm.search_engine);
      const DataQuery& query = data_queries.left.at(query_key);
      std::vector<MzTabDouble> rts(1);
      rts[0].set(query.rt);
      xsm.retention_time.set(rts);
      xsm.charge.set(match.charge);
      xsm.exp_mass_to_charge.set(query.mz);
      xsm.calc_mass_to_charge.set(calc_mass / abs(match.charge));
      xsm.spectra_ref.setMSFile(file_map[query.input_file_key]);
      xsm.spectra_ref.setSpecRef(query.data_id);
      // don't repeat data from the peptide section (e.g. accessions)
      // why are "pre"/"post"/"start"/"end" not in the peptide section?!
      output.push_back(xsm);
    }

    /// Helper function to add search engine scores to MzTab
    void exportScoresToMzTab_(const ScoreList& scores,
                              std::map<Size, MzTabDouble>& output,
                              std::map<ScoreTypeKey, Size>& score_map) const;

    /// Helper function to add processing steps (search engines) to MzTab
    void exportProcessingStepsToMzTab_(
      const std::vector<ProcessingStepKey>& steps, MzTabParameterList& output)
      const;

    /// Helper function to add search engine score entries to MzTab's meta data section
    void addMzTabSEScores_(const std::map<ScoreTypeKey, Size>& scores,
                           std::map<Size, MzTabParameter>& output) const;

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - oligonucleotide variant
    void addMzTabMoleculeParentContext_(
      const std::set<MoleculeParentMatch>& matches,
      const MzTabOligonucleotideSectionRow& molecule,
      std::vector<MzTabOligonucleotideSectionRow>& output) const;

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - peptide variant
    void addMzTabMoleculeParentContext_(
      const std::set<MoleculeParentMatch>& matches,
      const MzTabPeptideSectionRow& molecule,
      std::vector<MzTabPeptideSectionRow>& output) const;

    /// Helper function to compare two scores
    static bool isBetterScore_(double first, double second, bool higher_better)
    {
      if (higher_better) return first > second;
      return first < second;
    }

    /// Helper function to import DB search parameters from legacy format
    SearchParamsKey importDBSearchParameters_(
      const ProteinIdentification::SearchParameters& pisp);

    /// Helper function to export DB search parameters to legacy format
    ProteinIdentification::SearchParameters exportDBSearchParameters_(
      SearchParamsKey key) const;

    /// Helper function to check if all score types are valid and registered
    void checkScoreTypes_(const ScoreList& scores);

    /// Helper function to check if all processing steps are valid and registered
    void checkProcessingSteps_(const std::vector<ProcessingStepKey>& steps);

    /*!
      @brief Helper function to add the current processing step to a list of steps, if applicable.

      @see @ref setCurrentProcessingStep()
    */
    bool addCurrentProcessingStep_(
      std::vector<ProcessingStepKey>& processing_steps);
  };
}

#endif
