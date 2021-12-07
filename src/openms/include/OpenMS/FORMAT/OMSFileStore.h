// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <QtSql/QSqlQuery>

class QSqlError;

namespace OpenMS
{
  namespace Internal
  {
    /*!
      @brief Check if a specified database table exists

      @param db_name Name of the database (as used by Qt/QSqlDatabase)
      @param table_name Name of the table to check

      @return True if table exists, false if not
    */
    bool tableExists_(const String& db_name, const String& table_name);

    /*!
      @brief Raise a more informative database error

      Add context to an SQL error encountered by Qt and throw it as a FailedAPICall exception.

      @param error The error that occurred
      @param line Line in the code where error occurred
      @param function Name of the function where error occurred
      @param context Context for the error

      @throw Exception::FailedAPICall Throw this exception
    */
    void raiseDBError_(const QSqlError& error, int line,
                       const char* function, const String& context);


    /*!
      @brief Helper class for storing .oms files (SQLite format)

      This class encapsulates the SQLite database in a .oms file and allows to write data to it.
    */
    class OMSFileStore: public ProgressLogger
    {
    public:
      using Key = qint64; ///< Type used for database keys

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

      void createTableDataValue_();

      Key storeDataValue_(const DataValue& value);

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

      Key getAddress_(const IdentificationData::IdentifiedMolecule& molecule_var);

      void createTableParentMatches_();

      void storeParentMatches_(
        const IdentificationData::ParentMatches& matches, Key molecule_id);

      template<class MetaInfoInterfaceContainer>
      void storeMetaInfos_(const MetaInfoInterfaceContainer& container,
                           const String& parent_table)
      {
        bool table_created = false;
        QSqlQuery query; // prepare query only once and only if needed
        for (const auto& element : container)
        {
          if (!element.isMetaEmpty())
          {
            if (!table_created)
            {
              createTableMetaInfo_(parent_table);
              table_created = true;
            }
            storeMetaInfo_(element, parent_table, Key(&element));
          }
        }
      }

      template<class ScoredProcessingResultContainer>
      void storeScoredProcessingResults_(
        const ScoredProcessingResultContainer& container,
        const String& parent_table)
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
            for (const IdentificationData::AppliedProcessingStep& step :
                   element.steps_and_scores)
            {
              storeAppliedProcessingStep_(step, ++counter, parent_table,
                                          Key(&element));
            }
          }
        }
        storeMetaInfos_(container, parent_table);
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

      // store name, not database connection itself (see https://stackoverflow.com/a/55200682):
      QString db_name_;

      /// prepared queries for inserting data into different tables
      std::map<std::string, QSqlQuery> prepared_queries_;
    };
  }
}
