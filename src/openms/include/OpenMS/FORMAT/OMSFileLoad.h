// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OMSFileStore.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

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

        @param[in] filename Path to the .oms input file (SQLite database)
        @param[in] log_type Type of logging to use

        @throw Exception::FailedAPICall Database cannot be opened
      */
      OMSFileLoad(const std::string& filename, LogType log_type);

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

      /// Helper function for loading meta data on feature/consensus maps from the database
      template <class MapType> std::string loadMapMetaDataTemplate_(MapType& features);

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
      bool prepareQueryMetaInfo_(SQLite::Statement& query, const std::string& parent_table);

      /// Store results from an SQL query on meta values in a MetaInfoInterface(-derived) object
      void handleQueryMetaInfo_(SQLite::Statement& query, MetaInfoInterface& info,
                                Key parent_id);

      /// The database connection (read)
      std::unique_ptr<SQLite::Database> db_;

      int version_number_; ///< schema version number

      std::string subquery_score_; ///< query for score types used in JSON export

      void loadLegacyIdentifications_(IdentificationData& data);
      std::map<Key, IdentificationData::MoleculeIdentity> compatibility_molecules_;
      std::map<Key, IdentificationData::MatchReference> compatibility_matches_;
      const IdentificationData* identification_data_ = nullptr;
      // mapping: table name -> ordering critera (for JSON export)
      static std::map<std::string, std::string> export_order_by_;
    };
  }
}
