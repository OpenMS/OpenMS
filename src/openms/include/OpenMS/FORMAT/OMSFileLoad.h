// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/FORMAT/OMSFileStore.h>

#include <QtCore/QJsonArray> // for JSON export

namespace SQLite
{
  class Database;
} // namespace SQLite

namespace OpenMS
{
  class FeatureMap;
  class ConsensusMap;

  namespace Internal
  {
    /*!
      @brief Helper class for loading .oms files (SQLite format)

      This class encapsulates the SQLite database stored in a .oms file and allows to load data from it.
    */
    class OMSFileLoad: public ProgressLogger
    {
    public:
      using Key = OMSFileStore::Key; ///< Type used for database keys

      /*!
        @brief Constructor

        Opens the connection to the database file (in read-only mode).

        @param filename Path to the .oms input file (SQLite database)
        @param log_type Type of logging to use

        @throw Exception::FailedAPICall Database cannot be opened
      */
      OMSFileLoad(const String& filename, LogType log_type);

      /*!
        @brief Destructor

        Closes the connection to the database file.
      */
      ~OMSFileLoad();

      /// Load data from database and populate an IdentificationData object
      void load(IdentificationData& id_data);

      /// Load data from database and populate a FeatureMap object
      void load(FeatureMap& features);

      /// Load data from database and populate a ConsensusMap object
      void load(ConsensusMap& consensus);

      /// Export database contents in JSON format, write to stream
      void exportToJSON(std::ostream& output);

    private:
      /// Does the @p query contain an empty SQL statement (signifying that it shouldn't be executed)?
      static bool isEmpty_(const SQLite::Statement& query);

      /// Generate a DataValue with information returned by an SQL query
      static DataValue makeDataValue_(const SQLite::Statement& query);

      // static CVTerm loadCVTerm_(int id);

      /// Load information on score type from the database into IdentificationData
      void loadScoreTypes_(IdentificationData& id_data);

      /// Load information on input files from the database into IdentificationData
      void loadInputFiles_(IdentificationData& id_data);

      /// Load information on data processing software from the database into IdentificationData
      void loadProcessingSoftwares_(IdentificationData& id_data);

      /// Load information on sequence database search parameters from the database into IdentificationData
      void loadDBSearchParams_(IdentificationData& id_data);

      /// Load information on data processing steps from the database into IdentificationData
      void loadProcessingSteps_(IdentificationData& id_data);

      /// Load information on observations (e.g. spectra) from the database into IdentificationData
      void loadObservations_(IdentificationData& id_data);

      /// Load information on parent sequences (e.g. proteins) from the database into IdentificationData
      void loadParentSequences_(IdentificationData& id_data);

      /// Load information on parent group sets (e.g. protein groups) from the database into IdentificationData
      void loadParentGroupSets_(IdentificationData& id_data);

      /// Load information on identified compounds from the database into IdentificationData
      void loadIdentifiedCompounds_(IdentificationData& id_data);

      /// Load information on identified sequences (peptides or oligonucleotides) from the database into IdentificationData
      void loadIdentifiedSequences_(IdentificationData& id_data);

      /// Load information on adducts from the database into IdentificationData
      void loadAdducts_(IdentificationData& id_data);

      /// Load information on observation matches (e.g. PSMs) from the database into IdentificationData
      void loadObservationMatches_(IdentificationData& id_data);

      /// Helper function for loading meta data on feature/consensus maps from the database
      template <class MapType> String loadMapMetaDataTemplate_(MapType& features);

      /// Load feature map meta data from the database
      void loadMapMetaData_(FeatureMap& features);

      /// Load consensus map meta data from the database
      void loadMapMetaData_(ConsensusMap& consensus);

      /// Load information on data processing for feature/consensus maps from the database
      void loadDataProcessing_(std::vector<DataProcessing>& data_processing);

      /// Load information on features from the database into a feature map
      void loadFeatures_(FeatureMap& features);

      /// Generate a feature (incl. subordinate features) from data returned by SQL queries
      Feature loadFeatureAndSubordinates_(SQLite::Statement& query_feat,
                                          SQLite::Statement& query_meta,
                                          SQLite::Statement& query_match,
                                          SQLite::Statement& query_hull);

      /// Load consensus map column headers from the database
      void loadConsensusColumnHeaders_(ConsensusMap& consensus);

      /// Load information on consensus features from the database into a consensus map
      void loadConsensusFeatures_(ConsensusMap& consensus);

      /// Generate a BaseFeature (parent class) from data returned by SQL queries
      BaseFeature makeBaseFeature_(int id, SQLite::Statement& query_feat,
                                   SQLite::Statement& query_meta,
                                   SQLite::Statement& query_match);

      /// Prepare SQL queries for loading (meta) data on BaseFeatures from the database
      void prepareQueriesBaseFeature_(SQLite::Statement& query_meta,
                                      SQLite::Statement& query_match);

      /// Prepare SQL query for loading meta values associated with a particular class (stored in @p parent_table)
      bool prepareQueryMetaInfo_(SQLite::Statement& query, const String& parent_table);

      /// Store results from an SQL query on meta values in a MetaInfoInterface(-derived) object
      void handleQueryMetaInfo_(SQLite::Statement& query, MetaInfoInterface& info,
                                Key parent_id);

      /// Prepare SQL query for loading processing metadata associated with a particular class (stored in @p parent_table)
      bool prepareQueryAppliedProcessingStep_(SQLite::Statement& query,
                                              const String& parent_table);

      /// Store results from an SQL query on processing metadata in a ScoredProcessingResult(-derived) object
      void handleQueryAppliedProcessingStep_(
        SQLite::Statement& query,
        IdentificationDataInternal::ScoredProcessingResult& result,
        Key parent_id);

      /// Store results from an SQL query on parent matches
      void handleQueryParentMatch_(
        SQLite::Statement& query, IdentificationData::ParentMatches& parent_matches,
        Key molecule_id);

      /// Store results from an SQL query on peak annotations in an observation match
      void handleQueryPeakAnnotation_(
        SQLite::Statement& query, IdentificationData::ObservationMatch& match,
        Key parent_id);

      /// Export the contents of a database table to JSON
      QJsonArray exportTableToJSON_(const QString& table, const QString& order_by);

      /// The database connection (read)
      std::unique_ptr<SQLite::Database> db_;

      int version_number_; ///< schema version number

      QString subquery_score_; ///< query for score types used in JSON export

      // mappings between database keys and loaded data:
      std::unordered_map<Key, IdentificationData::ScoreTypeRef> score_type_refs_;
      std::unordered_map<Key, IdentificationData::InputFileRef> input_file_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingSoftwareRef> processing_software_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingStepRef> processing_step_refs_;
      std::unordered_map<Key, IdentificationData::SearchParamRef> search_param_refs_;
      std::unordered_map<Key, IdentificationData::ObservationRef> observation_refs_;
      std::unordered_map<Key, IdentificationData::ParentSequenceRef> parent_sequence_refs_;
      std::unordered_map<Key, IdentificationData::IdentifiedMolecule> identified_molecule_vars_;
      std::unordered_map<Key, IdentificationData::ObservationMatchRef> observation_match_refs_;
      std::unordered_map<Key, IdentificationData::AdductRef> adduct_refs_;

      // mapping: table name -> ordering critera (for JSON export)
      // @TODO: could use 'unordered_map' here, but would need to specify hash function for 'QString'
      static std::map<QString, QString> export_order_by_;
    };
  }
}
