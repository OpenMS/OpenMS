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

      @param error The error that occurred
      @param line Line in the code where error occurred
      @param function Name of the function where error occurred
      @param context Context for the error
      @param query Text of the query that was executed (optional)

      @throw Exception::FailedAPICall Throw this exception
    */
    void raiseDBError_(const String& error, int line, const char* function, const String& context, const String& query = "");

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

        @param filename Path to the .oms output file (SQLite database)
        @param log_type Type of logging to use

        @throw Exception::FailedAPICall Database cannot be opened
      */
      OMSFileStore(const String& filename, LogType log_type);

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

        @param name Name of the new table
        @param definition Table definition in SQL
        @param may_exist If true, the table may already exist (otherwise this is an error)
      */
      void createTable_(const String& name, const String& definition, bool may_exist = false);

      /// Create a database table for the data types used in DataValue
      void createTableDataValue_DataType_();

      /// Create a database table (and prepare a query) for storing CV terms
      void createTableCVTerm_();

      /// Create a database table (and prepare a query) for storing meta values
      void createTableMetaInfo_(const String& parent_table,
                                const String& key_column = "id");

      /// Store version information and current date/time in the database
      void storeVersionAndDate_();

      /// Store a CV term in the database
      Key storeCVTerm_(const CVTerm& cv_term);

      /// Store meta values (associated with one object) in the database
      void storeMetaInfo_(const MetaInfoInterface& info, const String& parent_table,
                          Key parent_id);

      /// Store meta values (for all objects in a container) in the database
      template<class MetaInfoInterfaceContainer, class DBKeyTable>
      void storeMetaInfos_(const MetaInfoInterfaceContainer& container,
                           const String& parent_table, const DBKeyTable& db_keys)
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

      /// @name Helper functions for storing identification data
      ///@{
      /// Store score type information from IdentificationData in the database
      void storeScoreTypes_(const IdentificationData& id_data);

      /// Store input file information from IdentificationData in the database
      void storeInputFiles_(const IdentificationData& id_data);

      /// Store information on data processing software from IdentificationData in the database
      void storeProcessingSoftwares_(const IdentificationData& id_data);

      /// Store sequence database search parameters from IdentificationData in the database
      void storeDBSearchParams_(const IdentificationData& id_data);

      /// Store information on data processing steps from IdentificationData in the database
      void storeProcessingSteps_(const IdentificationData& id_data);

      /// Store information on observations (e.g. spectra) from IdentificationData in the database
      void storeObservations_(const IdentificationData& id_data);

      /// Store information on parent sequences (e.g. proteins) from IdentificationData in the database
      void storeParentSequences_(const IdentificationData& id_data);

      /// Store information on parent group sets (e.g. protein groups) from IdentificationData in the database
      void storeParentGroupSets_(const IdentificationData& id_data);

      /// Store information on identified compounds from IdentificationData in the database
      void storeIdentifiedCompounds_(const IdentificationData& id_data);

      /// Store information on identified sequences (peptides or oligonucleotides) from IdentificationData in the database
      void storeIdentifiedSequences_(const IdentificationData& id_data);

      /// Store information on adducts from IdentificationData in the database
      void storeAdducts_(const IdentificationData& id_data);

      /// Store information on observation matches (e.g. PSMs) from IdentificationData in the database
      void storeObservationMatches_(const IdentificationData& id_data);

      /// Create a database table for molecule types (proteins, compounds, RNA)
      void createTableMoleculeType_();

      /// Create a database table for storing processing metadata
      void createTableAppliedProcessingStep_(const String& parent_table);

      /// Store processing metadata for a particular class (stored in @p parent_table) in the database
      void storeAppliedProcessingStep_(
        const IdentificationData::AppliedProcessingStep& step, Size step_order,
        const String& parent_table, Key parent_id);

      /// Create a database table for storing identified molecules (peptides, compounds, oligonucleotides)
      void createTableIdentifiedMolecule_();

      /// Return the database key used for an identified molecule (peptide, compound or oligonucleotide)
      Key getDatabaseKey_(const IdentificationData::IdentifiedMolecule& molecule_var);

      /// Create a database table for storing parent matches (e.g. proteins for a peptide)
      void createTableParentMatches_();

      /// Store information on parent matches in the database
      void storeParentMatches_(
        const IdentificationData::ParentMatches& matches, Key molecule_id);

      /// Store metadata on scores/processing steps (for all objects in a container) in the database
      template<class ScoredProcessingResultContainer, class DBKeyTable>
      void storeScoredProcessingResults_(const ScoredProcessingResultContainer& container,
                                         const String& parent_table, const DBKeyTable& db_keys)
      {
        bool table_created = false;
        for (const auto& element : container)
        {
          if (!element.steps_and_scores.empty())
          {
            if (!table_created)
            {
              createTableAppliedProcessingStep_(parent_table);
              table_created = true;
            }
            Size counter = 0;
            for (const IdentificationData::AppliedProcessingStep& step : element.steps_and_scores)
            {
              storeAppliedProcessingStep_(step, ++counter, parent_table, db_keys.at(&element));
            }
          }
        }
        storeMetaInfos_(container, parent_table, db_keys);
      }
      ///@}

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
      void storeMapMetaData_(const MapType& features, const String& experiment_type = "");

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

      // mapping between loaded data and database keys:
      // @NOTE: in principle we could use `unordered_map` here for efficiency,
      // but that gives compiler errors when pointers or iterators (`...Ref`)
      // are used as keys (because they aren't hashable?)
      std::map<const IdentificationData::ScoreType*, Key> score_type_keys_;
      std::map<const IdentificationData::InputFile*, Key> input_file_keys_;
      std::map<const IdentificationData::ProcessingSoftware*, Key> processing_software_keys_;
      std::map<const IdentificationData::ProcessingStep*, Key> processing_step_keys_;
      std::map<const IdentificationData::DBSearchParam*, Key> search_param_keys_;
      std::map<const IdentificationData::Observation*, Key> observation_keys_;
      std::map<const IdentificationData::ParentSequence*, Key> parent_sequence_keys_;
      std::map<const IdentificationData::ParentGroupSet*, Key> parent_grouping_keys_;
      std::map<const IdentificationData::IdentifiedCompound*, Key> identified_compound_keys_;
      std::map<const IdentificationData::IdentifiedPeptide*, Key> identified_peptide_keys_;
      std::map<const IdentificationData::IdentifiedOligo*, Key> identified_oligo_keys_;
      std::map<const AdductInfo*, Key> adduct_keys_;
      std::map<const IdentificationData::ObservationMatch*, Key> observation_match_keys_;
      // for feature/consensus maps:
      std::map<const DataProcessing*, Key> feat_processing_keys_;
    };
  }
}
