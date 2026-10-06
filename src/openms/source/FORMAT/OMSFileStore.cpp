// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include "OMSIdentificationData.h"

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/OMSFileStore.h>
#include <OpenMS/SYSTEM/File.h>
#include <SQLiteCpp/Database.h>
#include <SQLiteCpp/Transaction.h>
#include <sqlite3.h>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
constexpr int version_number = 6; // increase this whenever the DB schema changes!

void raiseDBError_(const std::string& error, int line, const char* function, const std::string& context, const std::string& query)
{
  std::string msg = context + ": " + error;
  if (! query.empty()) { msg += std::string("\nQuery was: ") + query; }
  throw Exception::FailedAPICall(__FILE__, line, function, msg);
}

  bool execAndReset(SQLite::Statement& query, int expected_modifications)
  {
    auto ret = query.exec();
    query.reset();
    return ret == expected_modifications;
  }

  void execWithExceptionAndReset(SQLite::Statement& query, int expected_modifications, int line, const char* function, const char* context)
  {
    if (!execAndReset(query, expected_modifications))
    {
      raiseDBError_(query.getErrorMsg(), line, function, context);
    }
  }

  OMSFileStore::OMSFileStore(const std::string& filename, LogType log_type)
  {
    setLogType(log_type);
    File::remove(filename); // nuke the file (SQLite cannot overwrite it)
    db_ = make_unique<SQLite::Database>(filename, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE); // throws on error
    // foreign key constraints are disabled by default - turn them on:
    // @TODO: performance impact? (seems negligible, but should be tested more)
    db_->exec("PRAGMA foreign_keys = ON");
    // disable synchronous filesystem access and the rollback journal to greatly
    // increase write performance - since we write a new output file every time,
    // we don't have to worry about database consistency:
    db_->exec("PRAGMA synchronous = OFF");
    db_->exec("PRAGMA journal_mode = OFF");
  }

  OMSFileStore::~OMSFileStore() = default;

  void OMSFileStore::createTable_(const std::string& name, const std::string& definition, bool may_exist)
  {
    std::string sql_create = "CREATE TABLE ";
    if (may_exist) sql_create += "IF NOT EXISTS ";
    sql_create += name + " (" + definition + ")";
    db_->exec(sql_create);
  }


  void OMSFileStore::storeVersionAndDate_()
  {
    createTable_("version",
                 "OMSFile INT NOT NULL, "       \
                 "date TEXT NOT NULL, "         \
                 "OpenMS TEXT, "                \
                 "build_date TEXT");

    SQLite::Statement query(*db_, "INSERT INTO version VALUES ("  \
                                 ":format_version, "             \
                                 "datetime('now'), "             \
                                 ":openms_version, "             \
                                 ":build_date)");
    query.bind(":format_version", version_number);
    query.bind(":openms_version", VersionInfo::getVersion());
    query.bind(":build_date", VersionInfo::getTime());
    query.exec();
  }


  void OMSFileStore::createTableDataValue_DataType_()
  {
    createTable_("DataValue_DataType",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "data_type TEXT UNIQUE NOT NULL");
    auto sql_insert =
      "INSERT INTO DataValue_DataType VALUES " \
      "(1, 'STRING_VALUE'), "                  \
      "(2, 'INT_VALUE'), "                     \
      "(3, 'DOUBLE_VALUE'), "                  \
      "(4, 'STRING_LIST'), "                   \
      "(5, 'INT_LIST'), "                      \
      "(6, 'DOUBLE_LIST')";
    db_->exec(sql_insert);
  }

  void OMSFileStore::createTableMetaInfo_(const std::string& parent_table, const std::string& key_column)
  {
    if (!db_->tableExists("DataValue_DataType")) createTableDataValue_DataType_();

    std::string parent_ref = parent_table + " (" + key_column + ")";
    std::string table = parent_table + "_MetaInfo";
    // for the data_type_id, empty values are represented using NULL
    createTable_(table, "parent_id INTEGER NOT NULL, "
                        "name TEXT NOT NULL, "
                        "data_type_id INTEGER, "
                        "value TEXT, unit_type INTEGER NOT NULL, unit INTEGER NOT NULL, "
                        "FOREIGN KEY (parent_id) REFERENCES "
                          + parent_ref
                          + ", "
                            "FOREIGN KEY (data_type_id) REFERENCES DataValue_DataType (id), "
                            "PRIMARY KEY (parent_id, name)");

    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO " + table
                                                        + " VALUES ("
                                                          ":parent_id, "
                                                          ":name, "
                                                          ":data_type_id, "
                                                          ":value, :unit_type, :unit)");
    prepared_queries_.emplace(table, std::move(query));
  }


  void OMSFileStore::storeMetaInfo_(const MetaInfoInterface& info, const std::string& parent_table, Key parent_id)
  {
    if (info.isMetaEmpty()) return;

    // this assumes the "..._MetaInfo" and "DataValue_DataType" tables exist already!
    auto& query = *prepared_queries_[parent_table + "_MetaInfo"];
    query.bind(":parent_id", parent_id);
    // this is inefficient, but MetaInfoInterface doesn't support iteration:
    vector<std::string> info_keys;
    info.getKeys(info_keys);
    for (const std::string& info_key : info_keys)
    {
      query.bind(":name", info_key);

      const DataValue& value = info.getMetaValue(info_key);
      if (value.isEmpty()) // use NULL as the type for empty values
      {
        query.bind(":data_type_id");
      }
      else
      {
        query.bind(":data_type_id", int(value.valueType()) + 1);
      }
      if (value.valueType() == DataValue::STRING_LIST)
      {
        std::string encoded;
        const auto append_integer = [&](UInt64 number) {
          for (Size i = 0; i < 8; ++i)
            encoded.push_back(static_cast<char>(number >> (i * 8)));
        };
        const auto values = value.toStringList();
        append_integer(values.size());
        for (const auto& item : values)
        {
          append_integer(item.size());
          encoded.append(item);
        }
        query.bind(":value", encoded.data(), static_cast<int>(encoded.size()));
      }
      else
        query.bind(":value", value.toString());
      query.bind(":unit_type", static_cast<int>(value.getUnitType()));
      query.bind(":unit", value.getUnit());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
    }
  }


  void OMSFileStore::store(const IdentificationData& data)
  {
    data.validate();
    identification_data_ = &data;
    const auto body = [&]() {
      storeVersionAndDate_();
      storeOMSIdentifications(*db_, data);
    };
    if (sqlite3_get_autocommit(db_->getHandle()))
    {
      SQLite::Transaction transaction(*db_);
      body();
      transaction.commit();
    }
    else
      body();
  }

  void OMSFileStore::createTableBaseFeature_(bool with_metainfo, bool with_idmatches)
  {
    createTable_("FEAT_BaseFeature",
                 "id INTEGER PRIMARY KEY NOT NULL, "
                 "rt REAL, "
                 "mz REAL, "
                 "intensity REAL, "
                 "charge INTEGER, "
                 "width REAL, "
                 "quality REAL, "
                 "unique_id INTEGER, "
                 "primary_encoding INTEGER, primary_representation TEXT, "
                 "subordinate_of INTEGER, "
                 "FOREIGN KEY (subordinate_of) REFERENCES FEAT_BaseFeature (id), "
                 "CHECK (id > subordinate_of)"); // check to prevent cycles

    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_BaseFeature VALUES ("
                                                      ":id, "
                                                      ":rt, "
                                                      ":mz, "
                                                      ":intensity, "
                                                      ":charge, "
                                                      ":width, "
                                                      ":quality, "
                                                      ":unique_id, "
                                                      ":primary_encoding, :primary_representation, "
                                                      ":subordinate_of)");
    prepared_queries_.emplace("FEAT_BaseFeature", std::move(query));

    if (with_metainfo)
    {
      createTableMetaInfo_("FEAT_BaseFeature");
    }
    if (with_idmatches)
    {
      createTable_("FEAT_ObservationMatch", "feature_id INTEGER NOT NULL, "
                                            "run_uuid TEXT NOT NULL, match_id TEXT NOT NULL, "
                                            "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
      query = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_ObservationMatch VALUES ("
                                                   ":feature_id, "
                                                   ":run_uuid, :match_id)");
      prepared_queries_.emplace("FEAT_ObservationMatch", std::move(query));
    }
    createTable_("FEAT_Query", "feature_id INTEGER NOT NULL, run_uuid TEXT NOT NULL, query_id TEXT NOT NULL, PRIMARY KEY(feature_id, run_uuid, "
                               "query_id), FOREIGN KEY(feature_id) REFERENCES FEAT_BaseFeature(id)");
    prepared_queries_["FEAT_Query"] = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_Query VALUES (:feature_id, :run_uuid, :query_id)");
  }


  void OMSFileStore::storeBaseFeature_(const BaseFeature& feature, int feature_id, int parent_id)
  {
    auto& query_feat = *prepared_queries_["FEAT_BaseFeature"];
    query_feat.bind(":id", feature_id);
    query_feat.bind(":rt", feature.getRT());
    query_feat.bind(":mz", feature.getMZ());
    query_feat.bind(":intensity", feature.getIntensity());
    query_feat.bind(":charge", feature.getCharge());
    query_feat.bind(":width", feature.getWidth());
    query_feat.bind(":quality", feature.getQuality());
    query_feat.bind(":unique_id", int64_t(feature.getUniqueId()));
    if (feature.hasPrimaryID())
    {
      query_feat.bind(":primary_encoding", static_cast<int>(feature.getPrimaryID().encoding));
      query_feat.bind(":primary_representation", feature.getPrimaryID().representation);
    }
    else // use NULL value
    {
      query_feat.bind(":primary_encoding");
      query_feat.bind(":primary_representation");
    }
    if (parent_id >= 0) // feature is a subordinate
    {
      query_feat.bind(":subordinate_of", parent_id);
    }
    else // use NULL value
    {
      query_feat.bind(":subordinate_of");
    }
    execWithExceptionAndReset(query_feat, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    // store ID observation matches:
    if (!feature.getIDMatches().empty())
    {
      auto& query_match = *prepared_queries_["FEAT_ObservationMatch"];
      query_match.bind(":feature_id", feature_id);
      for (const auto& ref : feature.getIDMatches())
      {
        const auto* run = identification_data_ ? identification_data_->findRunByUuid(ref.run_uuid) : nullptr;
        if (! run || ! run->findMatch(ref.match))
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature link refers to a missing match");
        query_match.bind(":run_uuid", ref.run_uuid);
        query_match.bind(":match_id", std::to_string(ref.match.value));
        execWithExceptionAndReset(query_match, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
    }

    auto& query_link = *prepared_queries_["FEAT_Query"];
    for (const auto& ref : feature.getIDQueries())
    {
      const auto* run = identification_data_ ? identification_data_->findRunByUuid(ref.run_uuid) : nullptr;
      if (! run || ! run->findIdentification(ref.query))
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature link refers to a missing query");
      query_link.bind(":feature_id", feature_id);
      query_link.bind(":run_uuid", ref.run_uuid);
      query_link.bind(":query_id", std::to_string(ref.query.value));
      execWithExceptionAndReset(query_link, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "Error storing feature query link");
    }
    storeMetaInfo_(feature, "FEAT_BaseFeature", feature_id);
  }


  void OMSFileStore::storeFeatureAndSubordinates_(
    const Feature& feature, int& feature_id, int parent_id)
  {
    storeBaseFeature_(feature, feature_id, parent_id);

    auto& query_feat = *prepared_queries_["FEAT_Feature"];
    query_feat.bind(":feature_id", feature_id);
    query_feat.bind(":rt_quality", feature.getQuality(0));
    query_feat.bind(":mz_quality", feature.getQuality(1));
    execWithExceptionAndReset(query_feat, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    // store convex hulls:
    const vector<ConvexHull2D>& hulls = feature.getConvexHulls();
    if (!hulls.empty())
    {
      auto& query_hull = *prepared_queries_["FEAT_ConvexHull"];
      query_hull.bind(":feature_id", feature_id);
      for (size_t i = 0; i < hulls.size(); ++i)
      {
        query_hull.bind(":hull_index", (int64_t)i);
        for (size_t j = 0; j < hulls[i].getHullPoints().size(); ++j)
        {
          const ConvexHull2D::PointType& point = hulls[i].getHullPoints()[j];
          query_hull.bind(":point_index", (int64_t)j);
          query_hull.bind(":point_x", point.getX());
          query_hull.bind(":point_y", point.getY());
          execWithExceptionAndReset(query_hull, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }

      }
    }
    // recurse into subordinates:
    parent_id = feature_id;
    ++feature_id; // variable is passed by reference, so effect is global
    for (const Feature& sub : feature.getSubordinates())
    {
      storeFeatureAndSubordinates_(sub, feature_id, parent_id);
    }
  }


  void OMSFileStore::storeFeatures_(const FeatureMap& features)
  {
    if (features.empty()) return;

    // create table(s) for BaseFeature parent class:
    // any meta infos on features?
    bool any_metainfo = anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.isMetaEmpty();
    });
    // any ID observations on features?
    bool any_idmatches = anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.getIDMatches().empty();
    });
    createTableBaseFeature_(any_metainfo, any_idmatches);

    createTable_("FEAT_Feature",
                 "feature_id INTEGER NOT NULL, "                        \
                 "rt_quality REAL, "                                    \
                 "mz_quality REAL, "                                    \
                 "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
    auto query = make_unique<SQLite::Statement>(*db_,
                                                "INSERT INTO FEAT_Feature VALUES (" \
                                                ":feature_id, "         \
                                                ":rt_quality, "         \
                                                ":mz_quality)");
    prepared_queries_.emplace("FEAT_Feature", std::move(query));

    // any convex hulls on features?
    if (anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.getConvexHulls().empty();
    }))
    {
      createTable_("FEAT_ConvexHull",
                   "feature_id INTEGER NOT NULL, "                      \
                   "hull_index INTEGER NOT NULL CHECK (hull_index >= 0), " \
                   "point_index INTEGER NOT NULL CHECK (point_index >= 0), " \
                   "point_x REAL, "                                     \
                   "point_y REAL, "                                     \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
      auto query2 = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_ConvexHull VALUES (" \
                                                   ":feature_id, "      \
                                                   ":hull_index, "      \
                                                   ":point_index, "     \
                                                   ":point_x, "         \
                                                   ":point_y)");
      prepared_queries_.emplace("FEAT_ConvexHull", std::move(query2));
    }

    // features and their subordinates are stored in DFS-like order:
    int feature_id = 0;
    for (const Feature& feat : features)
    {
      storeFeatureAndSubordinates_(feat, feature_id, -1);
      nextProgress();
    }
  }


  template <class MapType>
  void OMSFileStore::storeMapMetaData_(const MapType& features,
                                       const std::string& experiment_type)
  {
    createTable_("FEAT_MapMetaData",
                 "unique_id INTEGER PRIMARY KEY, "  \
                 "identifier TEXT, "                \
                 "file_path TEXT, "                 \
                 "file_type TEXT, "
                 "experiment_type TEXT"); // ConsensusMap only
    // @TODO: worth using a prepared query for just one insert?
    SQLite::Statement query(*db_,
                            "INSERT INTO FEAT_MapMetaData VALUES (" \
                            ":unique_id, "                          \
                            ":identifier, "                         \
                            ":file_path, "                          \
                            ":file_type, "                          \
                            ":experiment_type)");
    query.bind(":unique_id", int64_t(features.getUniqueId()));
    query.bind(":identifier", features.getIdentifier());
    query.bind(":file_path", features.getLoadedFilePath());
    std::string file_type = FileTypes::typeToName(features.getLoadedFileType());
    query.bind(":file_type", file_type);
    if (!experiment_type.empty())
    {
      query.bind(":experiment_type", experiment_type);
    }

    execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    if (!features.isMetaEmpty())
    {
      createTableMetaInfo_("FEAT_MapMetaData", "unique_id");
      storeMetaInfo_(features, "FEAT_MapMetaData", int64_t(features.getUniqueId()));
    }
  }

  // template specializations:
  template void OMSFileStore::storeMapMetaData_<FeatureMap>(const FeatureMap&, const std::string&);
  template void OMSFileStore::storeMapMetaData_<ConsensusMap>(const ConsensusMap&, const std::string&);


  void OMSFileStore::storeDataProcessing_(const vector<DataProcessing>& data_processing)
  {
    if (data_processing.empty()) return;

    createTable_("FEAT_DataProcessing",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "software_name TEXT, "                 \
                 "software_version TEXT, "              \
                 "processing_actions TEXT, "            \
                 "completion_time TEXT");
    // "id" is needed to connect to meta info table (see "storeMetaInfos_");
    // "position" is position in the vector ("index" is a reserved word in SQL)
    SQLite::Statement query(*db_, "INSERT INTO FEAT_DataProcessing VALUES (" \
                            ":id, "                                     \
                            ":software_name, "                          \
                            ":software_version, "                       \
                            ":processing_actions, "                     \
                            ":completion_time)");

    Key id = 1;
    for (const DataProcessing& proc : data_processing)
    {
      query.bind(":id", id);
      query.bind(":software_name", proc.getSoftware().getName());
      query.bind(":software_version", proc.getSoftware().getVersion());
      std::string actions;
      for (DataProcessing::ProcessingAction action : proc.getProcessingActions())
      {
        if (!actions.empty()) actions += ","; // @TODO: use different separator?
        actions += DataProcessing::NamesOfProcessingAction[action];
      }
      query.bind(":processing_actions", actions);
      query.bind(":completion_time", proc.getCompletionTime().get());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      feat_processing_keys_[&proc] = id;
      ++id;
    }
    storeMetaInfos_(data_processing, "FEAT_DataProcessing", feat_processing_keys_);
  }


  void OMSFileStore::store(const FeatureMap& features)
  {
    SQLite::Transaction transaction(*db_); // avoid SQLite's "implicit transactions", improve runtime
    store(features.getIdentificationData());
    startProgress(0, features.size() + 2, "Writing feature data to file");
    storeMapMetaData_(features);
    nextProgress();
    storeDataProcessing_(features.getDataProcessing());
    nextProgress();
    storeFeatures_(features);
    transaction.commit();
    endProgress();
  }


  void OMSFileStore::storeConsensusColumnHeaders_(const ConsensusMap& consensus)
  {
    if (consensus.getColumnHeaders().empty()) return; // shouldn't be empty in practice

    createTable_("FEAT_ConsensusColumnHeader",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "filename TEXT, "                      \
                 "label TEXT, "                         \
                 "size INTEGER, "                       \
                 "unique_id INTEGER");
    if (any_of(consensus.getColumnHeaders().begin(), consensus.getColumnHeaders().end(),
               [](const auto& pair){
                 return !pair.second.isMetaEmpty();
               }))
    {
      createTableMetaInfo_("FEAT_ConsensusColumnHeader");
    }

    SQLite::Statement query(*db_,
                            "INSERT INTO FEAT_ConsensusColumnHeader VALUES (" \
                            ":id, "                                     \
                            ":filename, "                               \
                            ":label, "                                  \
                            ":size, "                                   \
                            ":unique_id)");
    for (const auto& pair : consensus.getColumnHeaders())
    {
      Key id = int64_t(pair.first);
      query.bind(":id", id);
      query.bind(":filename", pair.second.filename);
      query.bind(":label", pair.second.label);
      query.bind(":size", int64_t(pair.second.size));
      query.bind(":unique_id", int64_t(pair.second.unique_id));

      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      storeMetaInfo_(pair.second, "FEAT_ConsensusColumnHeader", id);
    }
  }


  void OMSFileStore::storeConsensusFeatures_(const ConsensusMap& consensus)
  {
    if (consensus.empty()) return;

    // create table(s) for BaseFeature parent class:
    // any meta infos on features?
    bool any_metainfo = any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.isMetaEmpty();
    });
    // any ID observations on features?
    bool any_idmatches = any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.getIDMatches().empty();
    });
    createTableBaseFeature_(any_metainfo, any_idmatches);

    createTable_("FEAT_FeatureHandle",
                 "feature_id INTEGER NOT NULL, "                        \
                 "map_index INTEGER NOT NULL, "                                    \
                 "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
    SQLite::Statement query_handle(*db_,
                                   "INSERT INTO FEAT_FeatureHandle VALUES (" \
                                   ":feature_id, "                      \
                                   ":map_index)");

    // any ratios on consensus features?
    unique_ptr<SQLite::Statement> query_ratio; // only assign if needed below
    if (any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.getRatios().empty();
    }))
    {
      createTable_("FEAT_ConsensusRatio",
                   "feature_id INTEGER NOT NULL, "                      \
                   "ratio_index INTEGER NOT NULL CHECK (ratio_index >= 0), " \
                   "ratio_value REAL, "                                 \
                   "denominator_ref TEXT, "                             \
                   "numerator_ref TEXT, "                               \
                   "description TEXT, "                                 \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
      query_ratio = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_ConsensusRatio VALUES (" \
                                                   ":feature_id, "      \
                                                   ":ratio_index, "     \
                                                   ":ratio_value, "     \
                                                   ":denominator_ref, " \
                                                   ":numerator_ref, "   \
                                                   ":description)");
    }

    // consensus features and their subfeatures are stored in DFS-like order:
    int feature_id = 0;
    for (const ConsensusFeature& feat : consensus)
    {
      storeBaseFeature_(feat, feature_id, -1);
      int parent_id = feature_id;
      for (const FeatureHandle& handle : feat.getFeatures())
      {
        storeBaseFeature_(BaseFeature(handle), ++feature_id, parent_id);
        query_handle.bind(":feature_id", feature_id);
        query_handle.bind(":map_index", int64_t(handle.getMapIndex()));
        execWithExceptionAndReset(query_handle, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
      for (uint32_t i = 0; i < feat.getRatios().size(); ++i)
      {
        const ConsensusFeature::Ratio& ratio = feat.getRatios()[i];
        query_ratio->bind(":feature_id", feature_id);
        query_ratio->bind(":ratio_index", i);
        query_ratio->bind(":ratio_value", ratio.ratio_value_);
        query_ratio->bind(":denominator_ref", ratio.denominator_ref_);
        query_ratio->bind(":numerator_ref", ratio.numerator_ref_);
        query_ratio->bind(":description", ListUtils::concatenate(ratio.description_, ","));
        execWithExceptionAndReset(*query_ratio, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
      nextProgress();
      ++feature_id;
    }
  }


  void OMSFileStore::store(const ConsensusMap& consensus)
  {
    SQLite::Transaction transaction(*db_); // avoid SQLite's "implicit transactions", improve runtime
    store(consensus.getIdentificationData());
    startProgress(0, consensus.size() + 3, "Writing consensus feature data to file");
    storeMapMetaData_(consensus, consensus.getExperimentType());
    nextProgress();
    storeConsensusColumnHeaders_(consensus);
    nextProgress();
    storeDataProcessing_(consensus.getDataProcessing());
    nextProgress();
    storeConsensusFeatures_(consensus);
    transaction.commit();
    endProgress();
  }
}
