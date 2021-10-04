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

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlError>
#include <QtSql/QSqlQuery>

namespace OpenMS
{
  /**
      @brief This class supports reading and writing of OMS files.

      OMS files are SQLite databases consisting of several tables.
  */
  class OPENMS_DLLAPI OMSFile: public ProgressLogger
  {
  public:

    /// Constructor (with option to set log type)
    explicit OMSFile(LogType log_type = LogType::NONE):
      log_type_(log_type)
    {
      setLogType(log_type); // @TODO: move logging to OMSFileLoad
    }

    /** @brief Write out an IdentificationData object to SQL-based OMS file
     *
     * @param filename The output file
     * @param id_data The IdentificationData object
     */
    void store(const String& filename, const IdentificationData& id_data);

    /** @brief Write out a feature map to SQL-based OMS file
     *
     * @param filename The output file
     * @param features The feature map
     */
    void store(const String& filename, const FeatureMap& features);

    /** @brief Read in a OMS file and construct an IdentificationData object
     *
     * @param filename The input file
     * @param id_data The IdentificationData object
     */
    void load(const String& filename, IdentificationData& id_data);

    /** @brief Read in a OMS file and construct a feature map
     *
     * @param filename The input file
     * @param features The feature map
     */
    void load(const String& filename, FeatureMap& features);

  protected:

    using Key = qint64;

    LogType log_type_;

    // convenience functions:

    static bool tableExists_(const String& db_name, const String& table_name);

    static void raiseDBError_(const QSqlError& error, int line,
                              const char* function, const String& context);

    // store helper class:
    class OMSFileStore: public ProgressLogger
    {
    public:
      OMSFileStore(const String& filename, LogType log_type);

      ~OMSFileStore();

      void store(const IdentificationData& id_data);

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

      QSqlQuery getQueryMetaInfo_(const String& parent_table);

      void storeMetaInfo_(const MetaInfoInterface& info, Key parent_id,
                          QSqlQuery& query);

      void createTableAppliedProcessingStep_(const String& parent_table);

      void storeAppliedProcessingStep_(
        const IdentificationData::AppliedProcessingStep& step, Size step_order,
        const String& parent_table, Key parent_id);

      void createTableIdentifiedMolecule_();

      void createTableIdentifiedCompound_();

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
              query = getQueryMetaInfo_(parent_table);
            }
            storeMetaInfo_(element, Key(&element), query);
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
        const Feature& feature, int& feature_id, int parent_id,
        QSqlQuery& query_feat, QSqlQuery& query_meta, QSqlQuery& query_hull,
        QSqlQuery& query_match);

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
    };

    // load helper class:
    class OMSFileLoad: public ProgressLogger
    {
    public:
      OMSFileLoad(const String& filename, LogType log_type);

      ~OMSFileLoad();

      void load(IdentificationData& id_data);

      void load(FeatureMap& features);

    private:
      // static CVTerm loadCVTerm_(int id);

      void loadScoreTypes_(IdentificationData& id_data);

      void loadInputFiles_(IdentificationData& id_data);

      void loadProcessingSoftwares_(IdentificationData& id_data);

      void loadDBSearchParams_(IdentificationData& id_data);

      void loadProcessingSteps_(IdentificationData& id_data);

      void loadObservations_(IdentificationData& id_data);

      void loadParentSequences_(IdentificationData& id_data);

      void loadParentGroupSets_(IdentificationData& id_data);

      void loadIdentifiedCompounds_(IdentificationData& id_data);

      void loadIdentifiedSequences_(IdentificationData& id_data);

      void loadAdducts_(IdentificationData& id_data);

      void loadObservationMatches_(IdentificationData& id_data);

      void loadMapMetaData_(FeatureMap& features);

      void loadDataProcessing_(FeatureMap& features);

      void loadFeatures_(FeatureMap& features);

      Feature loadFeatureAndSubordinates_(QSqlQuery& query_feat,
                                          boost::optional<QSqlQuery>& query_meta,
                                          boost::optional<QSqlQuery>& query_hull,
                                          boost::optional<QSqlQuery>& query_match);

      static DataValue makeDataValue_(const QSqlQuery& query);

      bool prepareQueryMetaInfo_(QSqlQuery& query, const String& parent_table);

      void handleQueryMetaInfo_(QSqlQuery& query, MetaInfoInterface& info,
                                Key parent_id);

      bool prepareQueryAppliedProcessingStep_(QSqlQuery& query,
                                              const String& parent_table);

      void handleQueryAppliedProcessingStep_(
        QSqlQuery& query,
        IdentificationDataInternal::ScoredProcessingResult& result,
        Key parent_id);

      void handleQueryParentMatch_(
        QSqlQuery& query, IdentificationData::ParentMatches& parent_matches,
        Key molecule_id);

      void handleQueryPeakAnnotation_(
        QSqlQuery& query, IdentificationData::ObservationMatch& match,
        Key parent_id);

      // store name, not database connection itself (see https://stackoverflow.com/a/55200682):
      QString db_name_;

      std::unordered_map<Key, IdentificationData::ScoreTypeRef> score_type_refs_;
      std::unordered_map<Key, IdentificationData::InputFileRef> input_file_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingSoftwareRef> processing_software_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingStepRef> processing_step_refs_;
      std::unordered_map<Key, IdentificationData::SearchParamRef> search_param_refs_;
      std::unordered_map<Key, IdentificationData::ObservationRef> observation_refs_;
      std::unordered_map<Key, IdentificationData::ParentSequenceRef> parent_refs_;
      std::unordered_map<Key, IdentificationData::IdentifiedMolecule> identified_molecule_vars_;
      std::unordered_map<Key, IdentificationData::ObservationMatchRef> observation_match_refs_;
      std::unordered_map<Key, IdentificationData::AdductRef> adduct_refs_;
    };

  };
} // namespace OpenMS
