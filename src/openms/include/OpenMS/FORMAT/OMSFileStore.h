// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace SQLite
{
  class Database;
  class Exception;
  class Statement;
}

namespace OpenMS
{
  namespace Internal
  {
    /*!
      @brief Raise a more informative database error

      Add context to an SQL error encountered by Qt and throw it as a FailedAPICall exception.

      @param[in] error The error that occurred
      @param[in] line Line in the code where error occurred
      @param[in] function Name of the function where error occurred
      @param[in] context Context for the error
      @param[in] query Text of the query that was executed (optional)

      @throw Exception::FailedAPICall Throw this exception
    */
    void raiseDBError_(const std::string& error, int line, const char* function, const std::string& context, const std::string& query = "");

    /*!
      @brief Execute and reset an SQL query

      @return Whether the number of modifications made by the query matches the expected number
    */
    bool execAndReset(SQLite::Statement& query, int expected_modifications);

    /// If execAndReset() returns false, call raiseDBError_()
    void execWithExceptionAndReset(SQLite::Statement& query, int expected_modifications, int line, const char* function, const char* context);


    /*!
      @brief Helper class for storing .oms files (SQLite format)

      This class encapsulates the SQLite database in a .oms file and allows to write data to it.
    */
    class OMSFileStore: public ProgressLogger
    {
    public:
       ///< Type used for database keys
       using Key = int64_t; //std::decltype(((SQLite::Database*)nullptr)->getLastInsertRowid());

       /*!
        @brief Constructor

        Deletes the output file if it exists, then creates an SQLite database in its place.
        Opens the database and configures it for fast writing.

        @param[out] filename Path to the .oms output file (SQLite database)
        @param[in] log_type Type of logging to use

        @throw Exception::FailedAPICall Database cannot be opened
      */
      OMSFileStore(const std::string& filename, LogType log_type);

      /*!
        @brief Destructor

        Closes the connection to the database file.
      */
      ~OMSFileStore();

      /// Write data from an IdentificationData object to database
      void store(const IdentificationData& id_data);

      /// Write data from a FeatureMap object to database
      void store(const FeatureMap& features);

      /// Write data from a ConsensusMap object to database
      void store(const ConsensusMap& consensus);

    private:
      /*!
        @brief Helper function to create a database table

        @param[in] name Name of the new table
        @param[in] definition Table definition in SQL
        @param[in] may_exist If true, the table may already exist (otherwise this is an error)
      */
      void createTable_(const std::string& name, const std::string& definition, bool may_exist = false);

      /// Create a database table for the data types used in DataValue
      void createTableDataValue_DataType_();

      /// Create a database table (and prepare a query) for storing meta values
      void createTableMetaInfo_(const std::string& parent_table,
                                const std::string& key_column = "id");

      /// Store version information and current date/time in the database
      void storeVersionAndDate_();

      /// Store meta values (associated with one object) in the database
      void storeMetaInfo_(const MetaInfoInterface& info, const std::string& parent_table,
                          Key parent_id);

      /// Store meta values (for all objects in a container) in the database
      template<class MetaInfoInterfaceContainer, class DBKeyTable>
      void storeMetaInfos_(const MetaInfoInterfaceContainer& container,
                           const std::string& parent_table, const DBKeyTable& db_keys)
      {
        bool table_created = false;
        for (const auto& element : container)
        {
          if (!element.isMetaEmpty())
          {
            if (!table_created)
            {
              createTableMetaInfo_(parent_table);
              table_created = true;
            }
            storeMetaInfo_(element, parent_table, db_keys.at(&element));
          }
        }
      }

      /// @name Helper functions for storing (consensus) feature data
      ///@{
      /// Create a table for storing feature information
      void createTableBaseFeature_(bool with_metainfo, bool with_idmatches);

      /// Store information on a feature in the database
      void storeBaseFeature_(const BaseFeature& feature, int feature_id, int parent_id);

      /// Store information on features from a feature map in the database
      void storeFeatures_(const FeatureMap& features);

      /// Store a feature (incl. its subordinate features) in the database
      void storeFeatureAndSubordinates_(
        const Feature& feature, int& feature_id, int parent_id);

      /// check whether a predicate is true for any feature (or subordinate thereof) in a container
      template <class FeatureContainer, class Predicate>
      bool anyFeaturePredicate_(const FeatureContainer& features, const Predicate& pred)
      {
        if (features.empty()) return false;
        for (const Feature& feature : features)
        {
          if (pred(feature)) return true;
          if (anyFeaturePredicate_(feature.getSubordinates(), pred)) return true;
        }
        return false;
      }

      /// Store feature/consensus map meta data in the database
      template <class MapType>
      void storeMapMetaData_(const MapType& features, const std::string& experiment_type = "");

      /// Store information on data processing from a feature/consensus map in the database
      void storeDataProcessing_(const std::vector<DataProcessing>& data_processing);

      /// Store information on consensus features from a consensus map in the database
      void storeConsensusFeatures_(const ConsensusMap& consensus);

      /// Store information on column headers from a consensus map in the database
      void storeConsensusColumnHeaders_(const ConsensusMap& consensus);
      ///@}

      /// The database connection (read/write)
      std::unique_ptr<SQLite::Database> db_;

      /// Prepared queries for inserting data into different tables
      std::map<std::string, std::unique_ptr<SQLite::Statement>> prepared_queries_;

      const IdentificationData* identification_data_ = nullptr;
      // for feature/consensus maps:
      std::map<const DataProcessing*, Key> feat_processing_keys_;
    };
  }
}
