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

#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlError>
#include <QtSql/QSqlQuery>

namespace OpenMS
{
  class IdentificationData;
  /**
      @brief This class supports reading and writing of OMS files.

      OMS files are SQLite databases consisting of several tables.
  */
  class OPENMS_DLLAPI OMSFile
  {
  public:

    /** @brief Write out an IdentificationData object to SQL-based OMS file
     *
     * @param filename The output file
     * @param id_data The IdentificationData object
     *
    */
    void store(const String& filename, const IdentificationData& id_data);

    /** @brief Read in a OMS file and construct an IdentificationData
     *
     * @param filename The input file
     * @param id_data The IdentificationData object
     *
    */
    void load(const String& filename, IdentificationData& id_data);

  protected:

    using Key = qint64;

    // convenience functions:

    static bool tableExists_(const String& db_name, const String& table_name);

    static void raiseDBError_(const QSqlError& error, int line,
                              const char* function, const String& context);

    // store helper class:
    class OMSFileStore
    {
    public:
      OMSFileStore(const String& filename, const IdentificationData& id_data);

      ~OMSFileStore();

      void storeVersionAndDate();

      void storeScoreTypes();

      void storeInputFiles();

      void storeDataProcessingSoftwares();

      void storeDBSearchParams();

      void storeDataProcessingSteps();

      void storeDataQueries();

      void storeParentMolecules();

      void storeParentMoleculeGroupings();

      void storeIdentifiedCompounds();

      void storeIdentifiedSequences();

      void storeMoleculeQueryMatches();

    private:
      void createTable_(const String& name, const String& definition,
                        bool may_exist = false);

      void createTableMoleculeType_();

      void createTableDataValue_();

      Key storeDataValue_(const DataValue& value);

      void createTableCVTerm_();

      Key storeCVTerm_(const CVTerm& cv_term);

      void createTableMetaInfo_(const String& parent_table);

      void storeMetaInfo_(const MetaInfoInterface& info,
                          const String& parent_table, Key parent_id);

      void createTableAppliedProcessingStep_(const String& parent_table);

      void storeAppliedProcessingStep_(
        const IdentificationData::AppliedProcessingStep& step, Size step_order,
        const String& parent_table, Key parent_id);

      void createTableIdentifiedMolecule_();

      void createTableMoleculeParentMatches_();

      void storeMoleculeParentMatches_(
        const IdentificationData::ParentMatches& matches, Key molecule_id);

      template<class MetaInfoInterfaceContainer>
      void storeMetaInfos_(const MetaInfoInterfaceContainer& container,
                           const String& parent_table)
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

      // store name, not database connection itself (see https://stackoverflow.com/a/55200682):
      QString db_name_;

      const IdentificationData& id_data_;
    };

    // load helper class:
    class OMSFileLoad
    {
    public:
      OMSFileLoad(const String& filename, IdentificationData& id_data);

      ~OMSFileLoad();

      // static CVTerm loadCVTerm_(int id);

      void loadScoreTypes();

      void loadInputFiles();

      void loadDataProcessingSoftwares();

      void loadDBSearchParams();

      void loadDataProcessingSteps();

      void loadDataQueries();

      void loadParentMolecules();

      void loadParentMoleculeGroupings();

      void loadIdentifiedCompounds();

      void loadIdentifiedSequences();

      void loadMoleculeQueryMatches();

    private:
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
        QSqlQuery& query, IdentificationData::MoleculeQueryMatch& match,
        Key parent_id);

      // store name, not database connection itself (see https://stackoverflow.com/a/55200682):
      QString db_name_;

      IdentificationData& id_data_;

      std::unordered_map<Key, IdentificationData::ScoreTypeRef> score_type_refs_;
      std::unordered_map<Key, IdentificationData::InputFileRef> input_file_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingSoftwareRef> processing_software_refs_;
      std::unordered_map<Key, IdentificationData::ProcessingStepRef> processing_step_refs_;
      std::unordered_map<Key, IdentificationData::SearchParamRef> search_param_refs_;
      std::unordered_map<Key, IdentificationData::DataQueryRef> data_query_refs_;
      std::unordered_map<Key, IdentificationData::ParentMoleculeRef> parent_molecule_refs_;
      std::unordered_map<Key, IdentificationData::IdentifiedMoleculeRef> identified_molecule_refs_;
      std::unordered_map<Key, IdentificationData::QueryMatchRef> query_match_refs_;
    };

  };
} // namespace OpenMS


