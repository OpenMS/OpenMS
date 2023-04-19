// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
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

    bool execAndReset(SQLite::Statement& query, int expected_modifications);

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

    private:
      void storeVersionAndDate_();

      void storeScoreTypes_(const IdentificationData& id_data);

      void storeInputFiles_(const IdentificationData& id_data);

      void storeProcessingSoftwares_(const IdentificationData& id_data);

      void storeDBSearchParams_(const IdentificationData& id_data);

      void storeProcessingSteps_(const IdentificationData& id_data);

      void storeObservations_(const IdentificationData& id_data);

      void storeParentSequences_(const IdentificationData& id_data);

      void storeParentGroupSets_(const IdentificationData& id_data);

      void storeIdentifiedCompounds_(const IdentificationData& id_data);

      void storeIdentifiedSequences_(const IdentificationData& id_data);

      void storeAdducts_(const IdentificationData& id_data);

      void storeObservationMatches_(const IdentificationData& id_data);

      void storeFeatures_(const FeatureMap& features);

      void createTable_(const String& name, const String& definition,
                        bool may_exist = false);

      void createTableMoleculeType_();

      void createTableDataValue_DataType_();

      void createTableCVTerm_();

      Key storeCVTerm_(const CVTerm& cv_term);

      void createTableMetaInfo_(const String& parent_table,
                                const String& key_column = "id");

      void storeMetaInfo_(const MetaInfoInterface& info, const String& parent_table,
                          Key parent_id);

      void createTableAppliedProcessingStep_(const String& parent_table);

      void storeAppliedProcessingStep_(
        const IdentificationData::AppliedProcessingStep& step, Size step_order,
        const String& parent_table, Key parent_id);

      void createTableIdentifiedMolecule_();

      Key getDatabaseKey_(const IdentificationData::IdentifiedMolecule& molecule_var);

      void createTableParentMatches_();

      void storeParentMatches_(
        const IdentificationData::ParentMatches& matches, Key molecule_id);

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

      void storeFeature_(const FeatureMap& features);

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

      void storeMapMetaData_(const FeatureMap& features);

      void storeDataProcessing_(const FeatureMap& features);

      /// The database connection (read/write)
      std::unique_ptr<SQLite::Database> db_;

      /// prepared queries for inserting data into different tables
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
      // for feature maps:
      std::map<const DataProcessing*, Key> feat_processing_keys_;
    };
  }
}
