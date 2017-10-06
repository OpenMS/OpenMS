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

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Software.h>

#include <boost/bimap.hpp>
#include <boost/functional/hash.hpp> // for "unordered_map<pair<...>, ...>"
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{

  class OPENMS_DLLAPI IdentificationData: public MetaInfoInterface
  {
  protected:
    typedef UInt64 UniqueKey; // in case 64 bit isn't enough

    // Input files that were processed:
    typedef UniqueKey InputFileKey;
    typedef boost::bimap<InputFileKey, String> InputFileBimap;
    InputFileBimap input_files;


    /*!
      Data processing parameters.

      If the same processing is applied to multiple ID runs, e.g. if multiple files (fractions, replicates) are searched with the same search engine, store the
 parameters only once.
    */
    struct DataProcessingParameters: public MetaInfoInterface
    {
      Software tool;

      // @TODO: add processing actions that are relevant for ID data
      std::set<DataProcessing::ProcessingAction> actions;

      bool operator<(const DataProcessingParameters& other) const
      {
        return (std::tie(tool.getName(), tool.getVersion(), actions) <
                std::tie(other.tool.getName(), other.tool.getVersion(),
                         other.actions));
      }

      bool operator==(const DataProcessingParameters& other) const
      {
        return (std::tie(tool.getName(), tool.getVersion(), actions) ==
                std::tie(other.tool.getName(), other.tool.getVersion(),
                         other.actions));
      }
    };

    typedef UniqueKey ProcessingParamsKey;
    typedef boost::bimap<ProcessingParamsKey,
                         DataProcessingParameters> ParamsBimap;
    ParamsBimap processing_params;


    /*!
      Data processing step that is applied to the data (e.g. database search, PEP calculation, filtering, ConsensusID).
    */
    struct DataProcessingStep: public MetaInfoInterface
    {
      ProcessingParamsKey params_key;

      std::vector<String> ms_data_path; // path(s) to primary data

      std::vector<InputFileKey> input_files; // reference into "input_files"

      DateTime date_time;
    };

    typedef UniqueKey ProcessingStepKey;
    std::unordered_map<ProcessingStepKey, DataProcessingStep> processing_steps;


    /*!
      Information about a score type.
    */
    struct ScoreType: public MetaInfoInterface
    {
      CVTerm cv_term;

      String name;

      // Reference to the software that assigned the score:
      ProcessingParamsKey params_key;

      bool higher_better;

      bool operator<(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, params_key) <
                std::tie(other.cv_term.getAccession(), other.name,
                         other.params_key));
      }

      bool operator==(const ScoreType& other) const
      {
        return (std::tie(cv_term.getAccession(), name, params_key) ==
                std::tie(other.cv_term.getAccession(), other.name,
                         other.params_key));
      }
    };

    typedef UniqueKey ScoreTypeKey;
    typedef boost::bimap<ScoreTypeKey, ScoreType> ScoreTypeBimap;
    ScoreTypeBimap score_types;


    /*!
      Search query, e.g. spectrum or feature.
    */
    struct DataQuery: public MetaInfoInterface
    {
      InputFileKey input_file_key; // reference into "input_files"

      // spectrum or feature ID (from the file reference by "input_file_key"):
      String data_id;

      double rt, mz; // position

      // @TODO: ignore RT and m/z for comparisons?
      bool operator<(const DataQuery& other) const
        {
          return std::tie(input_file_key, data_id, rt, mz) <
            std::tie(other.input_file_key, other.data_id, other.rt, other.mz);
        }

      // @TODO: ignore RT and m/z for comparisons?
      bool operator==(const DataQuery& other) const
        {
          return std::tie(input_file_key, data_id, rt, mz) ==
            std::tie(other.input_file_key, other.data_id, other.rt, other.mz);
        }
    };

    typedef UniqueKey DataQueryKey;
    typedef boost::bimap<DataQueryKey, DataQuery> QueryBimap;
    QueryBimap data_queries;


    // Identified molecules - at the moment, peptides or small molecules:
    typedef UniqueKey IdentifiedMoleculeKey;
    typedef boost::bimap<IdentifiedMoleculeKey, AASequence> PeptideBimap;
    typedef boost::bimap<IdentifiedMoleculeKey, String> CompoundBimap;
    PeptideBimap identified_peptides;
    CompoundBimap identified_compounds;

    enum MoleculeType
    {
      MT_PROTEIN,
      MT_COMPOUND,
      // MT_RNA,
      SIZE_OF_MOLECULETYPES
    };

    /*!
      Meta data for an identified molecule.
    */
    struct IdentifiedMetaData: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;

      std::unordered_map<ScoreTypeKey, double> scores;

      std::vector<ProcessingStepKey> processing_steps;
    };

    std::unordered_map<IdentifiedMoleculeKey,
                       IdentifiedMetaData> identified_meta_data;


    // Evidence linking identified molecules to parent molecules (e.g. peptides
    // to proteins):
    // @TODO: rename "PeptideEvidence" to "MoleculeEvidence" (incl. members)
    typedef std::unordered_map<IdentifiedMoleculeKey,
                               std::vector<PeptideEvidence>> EvidenceMap;
    EvidenceMap parent_evidence;


    /*!
      Meta data specific to an identified compound (small molecule).
    */
    struct CompoundMetaData: public MetaInfoInterface
    {
      EmpiricalFormula formula;

      String smile;

      String inchi;
    };

    std::unordered_map<IdentifiedMoleculeKey,
                       CompoundMetaData> compound_meta_data;


    /*!
      Meta data for a search hit (e.g. peptide-spectrum match).
    */
    struct MatchMetaData: public MetaInfoInterface
    {
      Int charge;

      std::unordered_map<ScoreTypeKey, double> scores;

      // is it useful to store this, as different processing steps/score types
      // may give different rankings?
      Size rank; // rank among matches for the same spectrum/feature

      // ordered list of references to data processing steps:
      std::vector<ProcessingStepKey> processing_steps;

      // @TODO: move "PeakAnnotation" out of "PeptideHit"
      std::vector<PeptideHit::PeakAnnotation> peak_annotations;
    };

    // standard lib. doesn't include a hash function for pairs, but Boost does:
    typedef boost::hash<std::pair<DataQueryKey,
                                  IdentifiedMoleculeKey>> PairHash;
    typedef std::unordered_map<std::pair<DataQueryKey, IdentifiedMoleculeKey>,
                               MatchMetaData, PairHash> MatchMap;
    MatchMap matches;


    /*!
      Representation of a parent molecule that is identified only indirectly (e.g. a protein).
    */
    struct ParentMetaData: public MetaInfoInterface
    {
      enum MoleculeType molecule_type;

      String sequence;

      String description;

      double coverage;

      std::unordered_map<ScoreTypeKey, double> scores;

      // ordered list of references to data processing steps:
      std::vector<ProcessingStepKey> processing_steps;
    };

    typedef UniqueKey ParentMoleculeKey;
    typedef boost::bimap<ParentMoleculeKey, String> ParentBimap;
    ParentBimap parent_molecules;
    std::unordered_map<ParentMoleculeKey, ParentMetaData> parent_meta_data;


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


  public:

    void importIDs(const std::vector<ProteinIdentification>& proteins,
                   const std::vector<PeptideIdentification>& peptides);

    // "export" is a reserved keyword!
    void exportIDs(std::vector<ProteinIdentification>& proteins,
                   std::vector<PeptideIdentification>& peptides) const;
  };
}

#endif
