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

class QSqlQuery;

namespace OpenMS
{
  namespace Internal
  {
    /*!
      @brief Helper class for loading .oms files (SQLite format)

      This class encapsulates the SQLite database stored in a .oms file and allows to load data from it.
    */
    class OMSFileLoad: public ProgressLogger
    {
    public:
      using Key = qint64; ///< Type used for database keys

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
                                          std::optional<QSqlQuery>& query_meta,
                                          std::optional<QSqlQuery>& query_hull,
                                          std::optional<QSqlQuery>& query_match);

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

      int version_number_; ///< schema version number

      // mappings between database keys and loaded data:
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
  }
}
