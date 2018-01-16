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

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/variant.hpp>

namespace OpenMS
{
  class OPENMS_DLLAPI IdentificationData: public MetaInfoInterface
  {
  public:

    enum MoleculeType
    {
      MT_PROTEIN,
      MT_COMPOUND,
      MT_RNA,
      SIZE_OF_MOLECULETYPES
    };


    // Input files that were processed:
    typedef std::set<String> InputFiles;
    typedef const String* InputFileRef;
    InputFiles input_files;


    /*!
      Information about software used for data processing.

      If the same processing is applied to multiple ID runs, e.g. if multiple files (fractions, replicates) are searched with the same search engine, store the
 software information only once.
    */
    typedef std::set<Software> DataProcessingSoftware;
    typedef const Software* ProcessingSoftwareRef;
    DataProcessingSoftware processing_software;


    /*!
      Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct DataProcessingStep: public MetaInfoInterface
    {
      ProcessingSoftwareRef software_ref;

      std::vector<InputFileRef> input_file_refs;

      std::vector<String> primary_files; // path(s) to primary MS data

      DateTime date_time;

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      explicit DataProcessingStep(
        ProcessingSoftwareRef software_ref,
        const std::vector<InputFileRef>& input_file_refs =
        std::vector<InputFileRef>(), const std::vector<String>& primary_files =
        std::vector<String>(), const DateTime& date_time = DateTime::now(),
        std::set<DataProcessing::ProcessingAction> actions =
        std::set<DataProcessing::ProcessingAction>()):
        software_ref(software_ref), input_file_refs(input_file_refs),
        primary_files(primary_files), date_time(date_time), actions(actions)
      {
      }

      DataProcessingStep(const DataProcessingStep& other) = default;

      // don't compare meta data (?):
      bool operator<(const DataProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, primary_files,
                         date_time, actions) <
                std::tie(other.software_ref, other.input_file_refs,
                         other.primary_files, other.date_time, other.actions));
      }

      // don't compare meta data (?):
      bool operator==(const DataProcessingStep& other) const
      {
        return (std::tie(software_ref, input_file_refs, primary_files,
                         date_time, actions) ==
                std::tie(other.software_ref, other.input_file_refs,
                         other.primary_files, other.date_time, other.actions));
      }
    };

    typedef std::set<DataProcessingStep> DataProcessingSteps;
    typedef const DataProcessingStep* ProcessingStepRef;
    DataProcessingSteps processing_steps;


    /*!
      Parameters specific to a database search step.
    */
    struct DBSearchParam: public MetaInfoInterface
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

      DBSearchParam():
        molecule_type(MT_PROTEIN),
        peak_mass_type(ProteinIdentification::MONOISOTOPIC),
        precursor_mass_tolerance(0.0), fragment_mass_tolerance(0.0),
        precursor_tolerance_ppm(false), fragment_tolerance_ppm(false),
        digestion_enzyme(0), missed_cleavages(0), min_length(0), max_length(0)
      {
      }

      DBSearchParam(const DBSearchParam& other) = default;

      bool operator<(const DBSearchParam& other) const
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

      bool operator==(const DBSearchParam& other) const
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

    typedef std::set<DBSearchParam> DBSearchParams;
    typedef const DBSearchParam* SearchParamRef;
    DBSearchParams db_search_params;
    // @TODO: store SearchParamRef inside ProcessingStep? (may not be required
    // for many processing steps)
    std::map<ProcessingStepRef, SearchParamRef> db_search_steps;


    /*!
      Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term;

      String name;

      bool higher_better;

      // reference to the software that assigned the score:
      ProcessingSoftwareRef software_ref;
      // @TODO: scores assigned by different software tools/versions are
      // considered as different scores (even if they have the same name) -
      // does that make sense?

      ScoreType():
        higher_better(true), software_ref()
      {
      }

      explicit ScoreType(const CVTerm& cv_term, bool higher_better,
                         ProcessingSoftwareRef software_ref = nullptr):
        cv_term(cv_term), name(cv_term.getName()), higher_better(higher_better),
        software_ref(software_ref)
      {
      }

      explicit ScoreType(const String& name, bool higher_better,
                         ProcessingSoftwareRef software_ref = nullptr):
        cv_term(), name(name), higher_better(higher_better),
        software_ref(software_ref)
      {
      }

      ScoreType(const ScoreType& other) = default;

      // don't include "higher_better" in the comparison:
      bool operator<(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_ref) <
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_ref));
      }

      // don't include "higher_better" in the comparison:
      bool operator==(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, software_ref) ==
                std::tie(other.cv_term.getAccession(), other.name,
                         other.software_ref));
      }
    };

    typedef std::set<ScoreType> ScoreTypes;
    typedef const ScoreType* ScoreTypeRef;
    // @TODO: use a "boost::multi_index_container" to allow efficient access in
    // sequence and by key?
    typedef std::vector<std::pair<ScoreTypeRef, double>> ScoreList;
    ScoreTypes score_types;


    /*!
      Search query, e.g. spectrum or feature.
    */
    struct DataQuery: public MetaInfoInterface
    {
      // spectrum or feature ID (from the file referenced by "input_file_key"):
      String data_id;

      InputFileRef input_file_ref;

      double rt, mz; // position

      explicit DataQuery(
        const String& data_id,
        InputFileRef input_file_ref = nullptr,
        double rt = std::numeric_limits<double>::quiet_NaN(),
        double mz = std::numeric_limits<double>::quiet_NaN()):
        data_id(data_id), input_file_ref(input_file_ref), rt(rt), mz(mz)
      {
        // @TODO: require "input_file_ref"? (see also "DataProcessingStep")
      }

      DataQuery(const DataQuery& other) = default;

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator<(const DataQuery& other) const
      {
        return std::tie(input_file_ref, data_id) <
          std::tie(other.input_file_ref, other.data_id);
      }

      // ignore RT and m/z for comparisons to avoid issues with rounding:
      bool operator==(const DataQuery& other) const
      {
        return std::tie(input_file_ref, data_id) ==
          std::tie(other.input_file_ref, other.data_id);
      }

      // @TODO: do we need an "experiment label" (used e.g. in pepXML)?
      // if yes, should it be stored here or together with the input file?
    };

    typedef std::set<DataQuery> DataQueries;
    typedef const DataQuery* DataQueryRef;
    DataQueries data_queries;


    /// Base class for data with scores and processing steps (and meta info)
    struct ScoredProcessingResult: public MetaInfoInterface
    {
      ScoreList scores;

      // @TODO: use a "boost::multi_index_container" here for efficient look-up?
      std::vector<ProcessingStepRef> processing_step_refs;

      ScoredProcessingResult& operator+=(const ScoredProcessingResult& other)
      {
        // merge processing steps:
        for (auto step_ref : other.processing_step_refs)
        {
          if (std::find(processing_step_refs.begin(),
                        processing_step_refs.end(), step_ref) ==
              processing_step_refs.end())
          {
            processing_step_refs.push_back(step_ref);
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

      std::pair<double, bool> getScore(ScoreTypeRef score_ref) const
      {
        // give priority to "later" scores in the list:
        for (ScoreList::const_reverse_iterator it = scores.rbegin();
             it != scores.rend(); ++it)
        {
          if (it->first == score_ref) return std::make_pair(it->second, true);
        }
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), false);
      }

    protected:
      explicit ScoredProcessingResult(
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        scores(scores), processing_step_refs(processing_step_refs)
      {
      }

      ScoredProcessingResult(const ScoredProcessingResult& other) = default;
    };


    /*!
      Representation of a parent molecule that is identified only indirectly (e.g. a protein).
    */
    struct ParentMolecule: public ScoredProcessingResult
    {
      String accession;

      enum MoleculeType molecule_type; // @TODO: do we need this here?

      String sequence;

      String description;

      double coverage;

      bool is_decoy;

      explicit ParentMolecule(
        const String& accession, enum MoleculeType molecule_type = MT_PROTEIN,
        const String& sequence = "", const String& description = "",
        double coverage = 0.0, bool is_decoy = false,
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        ScoredProcessingResult(scores, processing_step_refs),
        accession(accession), molecule_type(molecule_type), sequence(sequence),
        description(description), coverage(coverage), is_decoy(is_decoy)
      {
      }

      ParentMolecule(const ParentMolecule& other) = default;

      ParentMolecule& operator+=(const ParentMolecule& other)
      {
        ScoredProcessingResult::operator+=(other);
        if (sequence.empty()) sequence = other.sequence;
        if (description.empty()) description = other.description;
        if (!is_decoy) is_decoy = other.is_decoy; // believe it when it's set
        // @TODO: what about coverage? (not reliable if we're merging data)

        return *this;
      }
    };

    // parent molecules indexed by their accessions:
    typedef boost::multi_index_container<
      ParentMolecule,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          ParentMolecule, String, &ParentMolecule::accession>>>
      > ParentMolecules;
    typedef const ParentMolecule* ParentMoleculeRef;
    ParentMolecules parent_molecules;


    /*!
      Meta data for the association between an identified molecule (e.g. peptide) and a parent molecule (e.g. protein).
    */
    struct MoleculeParentMatch: public MetaInfoInterface
    {
      // in extraordinary cases (e.g. database searches that allow insertions/
      // deletions), the length of the identified molecule may differ from the
      // length of the subsequence in the parent; therefore, store "end_pos":
      Size start_pos, end_pos;

      // @TODO: does "char" work here - what about modified ribonucleotides?
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

    // mapping: parent molecule -> match information
    typedef std::map<ParentMoleculeRef,
                     std::set<MoleculeParentMatch>> ParentMatches;


    // Identified molecules:
    template <typename SeqType>
    struct IdentifiedSequence: public ScoredProcessingResult
    {
      SeqType sequence;

      ParentMatches parent_matches;

      explicit IdentifiedSequence(
        const SeqType& sequence,
        const ParentMatches& parent_matches = ParentMatches(),
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        ScoredProcessingResult(scores, processing_step_refs),
        sequence(sequence), parent_matches(parent_matches)
      {
      }

      IdentifiedSequence(const IdentifiedSequence& other) = default;

      IdentifiedSequence& operator+=(const IdentifiedSequence& other)
      {
        ScoredProcessingResult::operator+=(other);
        // @TODO: improve merging of parent matches
        parent_matches.insert(other.parent_matches.begin(),
                              other.parent_matches.end());

        return *this;
      }

      bool allParentsAreDecoys() const
      {
        if (parent_matches.empty())
        {
          String msg = "no parent found for identified molecule";
          throw Exception::MissingInformation(__FILE__, __LINE__,
                                              OPENMS_PRETTY_FUNCTION, msg);
        }
        for (const auto& pair : parent_matches)
        {
          if (!pair.first->is_decoy) return false;
        }
        return true;
      }
    };

    typedef IdentifiedSequence<AASequence> IdentifiedPeptide;
    typedef IdentifiedSequence<NASequence> IdentifiedOligo;

    // identified peptides indexed by their sequences:
    typedef boost::multi_index_container<
      IdentifiedPeptide,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedPeptide, AASequence, &IdentifiedPeptide::sequence>>>
      > IdentifiedPeptides;
    typedef const IdentifiedPeptide* IdentifiedPeptideRef;
    IdentifiedPeptides identified_peptides;

    // identified oligos indexed by their sequences:
    typedef boost::multi_index_container<
      IdentifiedOligo,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedOligo, NASequence, &IdentifiedOligo::sequence>>>
      > IdentifiedOligos;
    typedef const IdentifiedOligo* IdentifiedOligoRef;
    IdentifiedOligos identified_oligos;

    struct IdentifiedCompound: public ScoredProcessingResult
    {
      String identifier;

      EmpiricalFormula formula;

      String name;

      String smile;

      String inchi;

      explicit IdentifiedCompound(
        const String& identifier,
        const EmpiricalFormula& formula = EmpiricalFormula(),
        const String& name = "", const String& smile = "",
        const String& inchi = "", const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>()):
        ScoredProcessingResult(scores, processing_step_refs),
        identifier(identifier), formula(formula), name(name), smile(smile),
        inchi(inchi)
      {
      }

      IdentifiedCompound(const IdentifiedCompound& other) = default;
    };

    // identified compounds indexed by their identifiers:
    typedef boost::multi_index_container<
      IdentifiedCompound,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::member<
          IdentifiedCompound, String, &IdentifiedCompound::identifier>>>
      > IdentifiedCompounds;
    typedef const IdentifiedCompound* IdentifiedCompoundRef;
    IdentifiedCompounds identified_compounds;


    /*!
      Meta data for a search hit (e.g. peptide-spectrum match).
    */

    // @TODO: move "PeakAnnotation" out of "PeptideHit"
    typedef std::vector<PeptideHit::PeakAnnotation> PeakAnnotations;

    typedef boost::variant<IdentifiedPeptideRef, IdentifiedCompoundRef,
                           IdentifiedOligoRef> IdentifiedMoleculeRef;

    struct MoleculeQueryMatch: public ScoredProcessingResult
    {
      IdentifiedMoleculeRef identified_molecule_ref;

      DataQueryRef data_query_ref;

      Int charge;

      // peak annotations (fragment ion matches), potentially from different
      // data processing steps:
      std::map<ProcessingStepRef, PeakAnnotations> peak_annotations;

      explicit MoleculeQueryMatch(
        IdentifiedMoleculeRef identified_molecule_ref,
        DataQueryRef data_query_ref, Int charge = 0,
        const ScoreList& scores = ScoreList(),
        const std::vector<ProcessingStepRef>& processing_step_refs =
        std::vector<ProcessingStepRef>(),
        const std::map<ProcessingStepRef, PeakAnnotations>& peak_annotations =
        std::map<ProcessingStepRef, PeakAnnotations>()):
        ScoredProcessingResult(scores, processing_step_refs),
        identified_molecule_ref(identified_molecule_ref),
        data_query_ref(data_query_ref), charge(charge),
        peak_annotations(peak_annotations)
      {
      }

      MoleculeQueryMatch(const MoleculeQueryMatch& other) = default;

      std::pair<DataQueryRef, IdentifiedMoleculeRef> getCombinedKey() const
      {
        return std::make_pair(data_query_ref, identified_molecule_ref);
      }

      enum MoleculeType getMoleculeType() const
      {
        if (boost::get<IdentifiedPeptideRef>(&identified_molecule_ref))
        {
          return MT_PROTEIN;
        }
        if (boost::get<IdentifiedCompoundRef>(&identified_molecule_ref))
        {
          return MT_COMPOUND;
        }
        if (boost::get<IdentifiedOligoRef>(&identified_molecule_ref))
        {
          return MT_RNA;
        }
        return SIZE_OF_MOLECULETYPES; // this shouldn't happen
      }

      IdentifiedPeptideRef getIdentifiedPeptideRef() const
      {
        if (const IdentifiedPeptideRef* ref_ptr =
            boost::get<IdentifiedPeptideRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a peptide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedCompoundRef getIdentifiedCompoundRef() const
      {
        if (const IdentifiedCompoundRef* ref_ptr =
            boost::get<IdentifiedCompoundRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not a compound";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      IdentifiedOligoRef getIdentifiedOligoRef() const
      {
        if (const IdentifiedOligoRef* ref_ptr =
            boost::get<IdentifiedOligoRef>(&identified_molecule_ref))
        {
          return *ref_ptr;
        }
        String msg = "matched molecule is not an oligonucleotide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      MoleculeQueryMatch& operator+=(const MoleculeQueryMatch& other)
      {
        ScoredProcessingResult::operator+=(other);
        if (charge == 0) charge = other.charge;
        peak_annotations.insert(other.peak_annotations.begin(),
                                other.peak_annotations.end());
        return *this;
      }
    };

    // all matches for the same data query should be consecutive!
    // tried using "boost::multi_index::composite_key" here, but that gave weird
    // compiler errors, so we define a combined key ourselves:
    typedef boost::multi_index_container<
      MoleculeQueryMatch,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<boost::multi_index::const_mem_fun<
          MoleculeQueryMatch, std::pair<DataQueryRef, IdentifiedMoleculeRef>,
          &MoleculeQueryMatch::getCombinedKey>>>
      > MoleculeQueryMatches;
    typedef const MoleculeQueryMatch* QueryMatchRef;
    MoleculeQueryMatches query_matches;


    /// Default constructor
    IdentificationData():
      current_step_ref_(nullptr)
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

    InputFileRef registerInputFile(const String& file);

    ProcessingSoftwareRef registerDataProcessingSoftware(
      const Software& software);

    SearchParamRef registerDBSearchParam(const DBSearchParam& param);

    ProcessingStepRef registerDataProcessingStep(
      const DataProcessingStep& step, SearchParamRef search_ref = nullptr);

    ScoreTypeRef registerScoreType(const ScoreType& score);

    DataQueryRef registerDataQuery(const DataQuery& query);

    IdentifiedPeptideRef registerPeptide(const IdentifiedPeptide& peptide);

    IdentifiedCompoundRef registerCompound(const IdentifiedCompound& compound);

    IdentifiedOligoRef registerOligo(const IdentifiedOligo& oligo);

    ParentMoleculeRef registerParentMolecule(const ParentMolecule& parent);

    QueryMatchRef registerMoleculeQueryMatch(const MoleculeQueryMatch& match);

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

    std::vector<QueryMatchRef> getBestMatchPerQuery(ScoreTypeRef
                                                    score_ref) const;

    ScoreTypeRef findScoreType(const String& score_name,
                               ProcessingSoftwareRef software_ref = nullptr)
      const;

    /// Helper function to compare two scores
    static bool isBetterScore(double first, double second, bool higher_better)
    {
      if (higher_better) return first > second;
      return first < second;
    }

  protected:

    /// Reference to the current data processing step (see @ref setCurrentProcessingStep())
    ProcessingStepRef current_step_ref_;

    /// Export a parent molecule (protein or nucleic acid) to mzTab
    template <typename MzTabSectionRow>
    void exportParentMoleculeToMzTab_(const ParentMolecule& parent,
                                      std::vector<MzTabSectionRow>& output,
                                      std::map<ScoreTypeRef, Size>& score_map)
      const
    {
      MzTabSectionRow row;
      row.accession.set(parent.accession);
      exportScoresToMzTab_(parent.scores, row.best_search_engine_score,
                           score_map);
      exportProcessingStepsToMzTab_(parent.processing_step_refs,
                                    row.search_engine);
      row.description.set(parent.description);
      row.coverage.set(parent.coverage);
      if (!parent.sequence.empty())
      {
        MzTabOptionalColumnEntry opt_seq;
        opt_seq.first = "opt_sequence";
        opt_seq.second.set(parent.sequence);
        row.opt_.push_back(opt_seq);
      }
      output.push_back(row);
    }

    /// Export an identified sequence (peptide or oligonucleotide, but not small molecule/compound) to mzTab
    template <typename MzTabSectionRow, typename IdentSeq>
    void exportPeptideOrOligoToMzTab_(const IdentSeq& identified,
                                      std::vector<MzTabSectionRow>& output,
                                      std::map<ScoreTypeRef, Size>& score_map)
      const
    {
      MzTabSectionRow row;
      // @TODO: handle modifications properly
      row.sequence.set(identified.sequence.toString());
      exportScoresToMzTab_(identified.scores, row.best_search_engine_score,
                           score_map);
      exportProcessingStepsToMzTab_(identified.processing_step_refs,
                                    row.search_engine);
      // generate one entry (with duplicated data) for every accession:
      bool unique = (identified.parent_matches.size() == 1);
      for (const std::pair<ParentMoleculeRef, std::set<MoleculeParentMatch>>&
             match_pair : identified.parent_matches)
      {
        const String& accession = match_pair.first->accession;
        row.accession.set(accession);
        row.unique.set(unique);
        if (match_pair.second.empty())
        {
          output.push_back(row);
        }
        else
        {
          addMzTabMoleculeParentContext_(match_pair.second, row, output);
        }
      }
      if (identified.parent_matches.empty())
      {
        // row.unique.set(false); // leave this unset?
        output.push_back(row);
      }
    }

    /// Export a molecule-query match (peptide- or oligonucleotide-spectrum match) to mzTab
    template <typename MzTabSectionRow>
    void exportQueryMatchToMzTab_(const String& sequence,
                                  const MoleculeQueryMatch& match,
                                  double calc_mass,
                                  std::vector<MzTabSectionRow>& output,
                                  std::map<ScoreTypeRef, Size>& score_map,
                                  std::map<InputFileRef, Size>& file_map)
      const
    {
      MzTabSectionRow xsm; // PSM or OSM
      // @TODO: handle modifications properly
      xsm.sequence.set(sequence);
      exportScoresToMzTab_(match.scores, xsm.search_engine_score, score_map);
      exportProcessingStepsToMzTab_(match.processing_step_refs,
                                    xsm.search_engine);
      const DataQuery& query = *match.data_query_ref;
      std::vector<MzTabDouble> rts(1);
      rts[0].set(query.rt);
      xsm.retention_time.set(rts);
      xsm.charge.set(match.charge);
      xsm.exp_mass_to_charge.set(query.mz);
      xsm.calc_mass_to_charge.set(calc_mass / abs(match.charge));
      if (query.input_file_ref)
      {
        xsm.spectra_ref.setMSFile(file_map[query.input_file_ref]);
      }
      xsm.spectra_ref.setSpecRef(query.data_id);
      // don't repeat data from the peptide section (e.g. accessions)
      // why are "pre"/"post"/"start"/"end" not in the peptide section?!
      output.push_back(xsm);
    }

    /// Helper function to add search engine scores to MzTab
    void exportScoresToMzTab_(const ScoreList& scores,
                              std::map<Size, MzTabDouble>& output,
                              std::map<ScoreTypeRef, Size>& score_map) const;

    /// Helper function to add processing steps (search engines) to MzTab
    void exportProcessingStepsToMzTab_(
      const std::vector<ProcessingStepRef>& steps, MzTabParameterList& output)
      const;

    /// Helper function to add search engine score entries to MzTab's meta data section
    void addMzTabSEScores_(const std::map<ScoreTypeRef, Size>& scores,
                           std::map<Size, MzTabParameter>& output) const;

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - oligonucleotide variant
    void addMzTabMoleculeParentContext_(
      const std::set<MoleculeParentMatch>& matches,
      const MzTabOligonucleotideSectionRow& row,
      std::vector<MzTabOligonucleotideSectionRow>& output) const;

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - peptide variant
    void addMzTabMoleculeParentContext_(
      const std::set<MoleculeParentMatch>& matches,
      const MzTabPeptideSectionRow& row,
      std::vector<MzTabPeptideSectionRow>& output) const;

    /// Helper function to import DB search parameters from legacy format
    SearchParamRef importDBSearchParameters_(
      const ProteinIdentification::SearchParameters& pisp);

    /// Helper function to export DB search parameters to legacy format
    ProteinIdentification::SearchParameters exportDBSearchParameters_(
      SearchParamRef ref) const;

    /// Helper function to check if all score types are valid
    void checkScoreTypes_(const ScoreList& scores);

    /// Helper function to check if all processing steps are valid
    void checkProcessingSteps_(const std::vector<ProcessingStepRef>& step_refs);

    /// Helper function to check if all parent matches are valid
    void checkParentMatches_(const ParentMatches& matches,
                             enum MoleculeType expected_type);

    /*!
      @brief Helper function to add the current processing step to a list of steps, if applicable.

      @see @ref setCurrentProcessingStep()
    */
    bool addCurrentProcessingStep_(
      std::vector<ProcessingStepRef>& processing_step_refs);

    /// Helper functor for augmenting entries (derived from ScoredProcessingResult) in a @t boost::multi_index_container structure
    template <typename ElementType>
    struct ModifyMultiIndexMergeElements
    {
      ModifyMultiIndexMergeElements(const ElementType& update):
        update(update)
      {
      }

      void operator()(ElementType& element)
      {
        element += update;
      }

      const ElementType& update;
    };

    /// Helper functor for adding the current processing step to elements in a @t boost::multi_index_container structure
    template <typename ElementType>
    struct ModifyMultiIndexAddProcessingStep
    {
      ModifyMultiIndexAddProcessingStep(ProcessingStepRef step_ref):
        step_ref(step_ref)
      {
      }

      void operator()(ElementType& element)
      {
        if ((step_ref != nullptr) &&
            (element.processing_step_refs.empty() ||
             (element.processing_step_refs.back() != step_ref)))
        {
          element.processing_step_refs.push_back(step_ref);
        }
      }

      ProcessingStepRef step_ref;
    };

    /// Helper function for adding entries (derived from ScoredProcessingResult) to a @t boost::multi_index_container structure
    template <typename ContainerType, typename ElementType>
    const ElementType* insertIntoMultiIndex_(ContainerType& container,
                                             const ElementType& element)
    {
      checkScoreTypes_(element.scores);
      checkProcessingSteps_(element.processing_step_refs);

      auto result = container.insert(element);
      if (!result.second) // existing element - merge in new information
      {
        ModifyMultiIndexMergeElements<ElementType> modifier(element);
        container.modify(result.first, modifier);
      }

      // add current processing step (if necessary):
      ModifyMultiIndexAddProcessingStep<ElementType>
        modifier(current_step_ref_);
      container.modify(result.first, modifier);

      return &(*result.first);
    }


    /// Check whether a pointer references an element in a container
    template <typename RefType, typename ContainerType>
    bool isValidReference_(RefType ref, const ContainerType& container)
    {
      for (const auto& element : container)
      {
        if (ref == &element) return true;
      }
      return false;
    }

  };
}

#endif
