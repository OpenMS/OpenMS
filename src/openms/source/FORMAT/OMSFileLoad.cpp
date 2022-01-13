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

#include <OpenMS/FORMAT/OMSFileLoad.h>
#include <OpenMS/FORMAT/OMSFileStore.h> // for "tableExists_" and "raiseDBError_"
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlError>
#include <QtSql/QSqlQuery>
// strangely, this is needed for type conversions in "QSqlQuery::bindValue":
#include <QtSql/QSqlQueryModel>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
  OMSFileLoad::OMSFileLoad(const String& filename, LogType log_type):
    db_name_("load_" + filename.toQString() + "_" + QString::number(UniqueIdGenerator::getUniqueId()))
  {
    setLogType(log_type);

    // open database:
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", db_name_);
    db.setDatabaseName(filename.toQString());
    db.setConnectOptions("QSQLITE_OPEN_READONLY");
    if (!db.open())
    {
      raiseDBError_(db.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error opening SQLite database");
      // if d'tor doesn't get called, DB connection (db_name_) doesn't get
      // removed, but that shouldn't be a big problem
    }
  }


  OMSFileLoad::~OMSFileLoad()
  {
    QSqlDatabase::database(db_name_).close();
    QSqlDatabase::removeDatabase(db_name_);
  }


  // currently not needed:
  // CVTerm OMSFileLoad::loadCVTerm_(int id)
  // {
  //   // this assumes that the "CVTerm" table exists!
  //   QSqlQuery query(db_);
  //   query.setForwardOnly(true);
  //   QString sql_select = "SELECT * FROM CVTerm WHERE id = " + QString(id);
  //   if (!query.exec(sql_select) || !query.next())
  //   {
  //     raiseDBError_(model.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
  //                   "error reading from database");
  //   }
  //   return CVTerm(query.value("accession").toString(),
  //                 query.value("name").toString(),
  //                 query.value("cv_identifier_ref").toString());
  // }


  void OMSFileLoad::loadScoreTypes_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ScoreType")) return;
    if (!tableExists_(db_name_, "CVTerm")) // every score type is a CV term
    {
      String msg = "required database table 'CVTerm' not found";
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    // careful - both joined tables have an "id" field, need to exclude one:
    if (!query.exec("SELECT S.*, C.accession, C.name, C.cv_identifier_ref " \
                    "FROM ID_ScoreType AS S JOIN CVTerm AS C "          \
                    "ON S.cv_term_id = C.id"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      CVTerm cv_term(query.value("accession").toString(),
                     query.value("name").toString(),
                     query.value("cv_identifier_ref").toString());
      bool higher_better = query.value("higher_better").toInt();
      ID::ScoreType score_type(cv_term, higher_better);
      ID::ScoreTypeRef ref = id_data.registerScoreType(score_type);
      score_type_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFileLoad::loadInputFiles_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_InputFile")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_InputFile"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::InputFile input(query.value("name").toString(),
                          query.value("experimental_design_id").toString());
      String primary_files = query.value("primary_files").toString();
      vector<String> pf_list = ListUtils::create<String>(primary_files);
      input.primary_files.insert(pf_list.begin(), pf_list.end());
      ID::InputFileRef ref = id_data.registerInputFile(input);
      input_file_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFileLoad::loadProcessingSoftwares_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ProcessingSoftware")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ProcessingSoftware"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    bool have_scores = tableExists_(db_name_,
                                    "ID_ProcessingSoftware_AssignedScore");
    QSqlQuery subquery(db);
    if (have_scores)
    {
      subquery.setForwardOnly(true);
      subquery.prepare("SELECT score_type_id "                         \
                       "FROM ID_ProcessingSoftware_AssignedScore " \
                       "WHERE software_id = :id ORDER BY score_type_order ASC");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      ID::ProcessingSoftware software(query.value("name").toString(),
                                          query.value("version").toString());
      if (have_scores)
      {
        subquery.bindValue(":id", id);
        if (!subquery.exec())
        {
          raiseDBError_(subquery.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                        "error reading from database");
        }
        while (subquery.next())
        {
          Key score_type_id = subquery.value(0).toLongLong();
          software.assigned_scores.push_back(score_type_refs_[score_type_id]);
        }
      }
      ID::ProcessingSoftwareRef ref = id_data.registerProcessingSoftware(software);
      processing_software_refs_[id] = ref;
    }
  }


  DataValue OMSFileLoad::makeDataValue_(const QSqlQuery& query)
  {
    DataValue::DataType type = DataValue::EMPTY_VALUE;
    int type_index = query.value("data_type_id").toInt();
    if (type_index > 0) type = DataValue::DataType(type_index - 1);
    String value = query.value("value").toString();
    switch (type)
    {
    case DataValue::STRING_VALUE:
      return DataValue(value);
    case DataValue::INT_VALUE:
      return DataValue(value.toInt());
    case DataValue::DOUBLE_VALUE:
      return DataValue(value.toDouble());
    // converting lists to String adds square brackets - remove them:
    case DataValue::STRING_LIST:
      value = value.substr(1, value.size() - 2);
      return DataValue(ListUtils::create<String>(value));
    case DataValue::INT_LIST:
      value = value.substr(1, value.size() - 2);
      return DataValue(ListUtils::create<int>(value));
    case DataValue::DOUBLE_LIST:
      value = value.substr(1, value.size() - 2);
      return DataValue(ListUtils::create<double>(value));
    default: // DataValue::EMPTY_VALUE (avoid warning about missing return)
      return DataValue();
    }
  }


  bool OMSFileLoad::prepareQueryMetaInfo_(QSqlQuery& query,
                                          const String& parent_table)
  {
    String table_name = parent_table + "_MetaInfo";
    if (!tableExists_(db_name_, table_name)) return false;

    query.setForwardOnly(true);
    QString sql_select =
      "SELECT * FROM " + table_name.toQString() + " AS MI " \
      "JOIN DataValue AS DV ON MI.data_value_id = DV.id "   \
      "WHERE MI.parent_id = :id";
    query.prepare(sql_select);
    return true;
  }


  bool OMSFileLoad::prepareQueryAppliedProcessingStep_(QSqlQuery& query,
                                                       const String& parent_table)
  {
    String table_name = parent_table + "_AppliedProcessingStep";
    if (!tableExists_(db_name_, table_name)) return false;

    // query.setForwardOnly(true);
    QString sql_select = "SELECT * FROM " + table_name.toQString() +
      " WHERE parent_id = :id ORDER BY processing_step_order ASC";
    query.prepare(sql_select);
    return true;
  }


  void OMSFileLoad::handleQueryMetaInfo_(QSqlQuery& query,
                                         MetaInfoInterface& info,
                                         Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      DataValue value = makeDataValue_(query);
      info.setMetaValue(query.value("name").toString(), value);
    }
  }


  void OMSFileLoad::handleQueryAppliedProcessingStep_(
    QSqlQuery& query,
    IdentificationDataInternal::ScoredProcessingResult& result,
    Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::AppliedProcessingStep step;
      QVariant step_id_opt = query.value("processing_step_id");
      if (!step_id_opt.isNull())
      {
        step.processing_step_opt =
          processing_step_refs_[step_id_opt.toLongLong()];
      }
      QVariant score_type_opt = query.value("score_type_id");
      if (!score_type_opt.isNull())
      {
        step.scores[score_type_refs_[score_type_opt.toLongLong()]] =
          query.value("score").toDouble();
      }
      result.addProcessingStep(step); // this takes care of merging the steps
    }
  }


  void OMSFileLoad::loadDBSearchParams_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_DBSearchParam")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_DBSearchParam"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      ID::DBSearchParam param;
      int molecule_type_index = query.value("molecule_type_id").toInt() - 1;
      param.molecule_type = ID::MoleculeType(molecule_type_index);
      int mass_type_index = query.value("mass_type_average").toInt();
      param.mass_type = ID::MassType(mass_type_index);
      param.database = query.value("database").toString();
      param.database_version = query.value("database_version").toString();
      param.taxonomy = query.value("taxonomy").toString();
      vector<Int> charges =
        ListUtils::create<Int>(query.value("charges").toString());
      param.charges.insert(charges.begin(), charges.end());
      vector<String> fixed_mods =
        ListUtils::create<String>(query.value("fixed_mods").toString());
      param.fixed_mods.insert(fixed_mods.begin(), fixed_mods.end());
      vector<String> variable_mods =
        ListUtils::create<String>(query.value("variable_mods").toString());
      param.variable_mods.insert(variable_mods.begin(), variable_mods.end());
      param.precursor_mass_tolerance =
        query.value("precursor_mass_tolerance").toDouble();
      param.fragment_mass_tolerance =
        query.value("fragment_mass_tolerance").toDouble();
      param.precursor_tolerance_ppm =
        query.value("precursor_tolerance_ppm").toInt();
      param.fragment_tolerance_ppm =
        query.value("fragment_tolerance_ppm").toInt();
      String enzyme = query.value("digestion_enzyme").toString();
      if (!enzyme.empty())
      {
        if (param.molecule_type == ID::MoleculeType::PROTEIN)
        {
          param.digestion_enzyme = ProteaseDB::getInstance()->getEnzyme(enzyme);
        }
        else if (param.molecule_type == ID::MoleculeType::RNA)
        {
          param.digestion_enzyme = RNaseDB::getInstance()->getEnzyme(enzyme);
        }
      }
      param.missed_cleavages = query.value("missed_cleavages").toUInt();
      param.min_length = query.value("min_length").toUInt();
      param.max_length = query.value("max_length").toUInt();
      ID::SearchParamRef ref = id_data.registerDBSearchParam(param);
      search_param_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadProcessingSteps_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ProcessingStep")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ProcessingStep"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_file(db);
    bool have_input_files = tableExists_(db_name_,
                                         "ID_ProcessingStep_InputFile");
    if (have_input_files)
    {
      subquery_file.setForwardOnly(true);
      subquery_file.prepare("SELECT input_file_id "                 \
                            "FROM ID_ProcessingStep_InputFile " \
                            "WHERE processing_step_id = :id");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "ID_ProcessingStep");
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      Key software_id = query.value("software_id").toLongLong();
      ID::ProcessingStep step(processing_software_refs_[software_id]);
      String date_time = query.value("date_time").toString();
      if (!date_time.empty()) step.date_time.set(date_time);
      if (have_input_files)
      {
        subquery_file.bindValue(":id", id);
        if (!subquery_file.exec())
        {
          raiseDBError_(subquery_file.lastError(), __LINE__,
                        OPENMS_PRETTY_FUNCTION, "error reading from database");
        }
        while (subquery_file.next())
        {
          Key input_file_id = subquery_file.value(0).toLongLong();
          // the foreign key constraint should ensure that look-up succeeds:
          step.input_file_refs.push_back(input_file_refs_[input_file_id]);
        }
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, step, id);
      }
      ID::ProcessingStepRef ref;
      QVariant opt_search_param_id = query.value("search_param_id");
      if (opt_search_param_id.isNull()) // no DB search params available
      {
        ref = id_data.registerProcessingStep(step);
      }
      else
      {
        ID::SearchParamRef search_param_ref =
          search_param_refs_[opt_search_param_id.toLongLong()];
        ref = id_data.registerProcessingStep(step, search_param_ref);
      }
      processing_step_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadObservations_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_Observation")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_Observation"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_Observation");

    while (query.next())
    {
      QVariant input_file_id = query.value("input_file_id");
      ID::Observation obs(query.value("data_id").toString(),
                          input_file_refs_[input_file_id.toLongLong()]);
      QVariant rt = query.value("rt");
      if (!rt.isNull()) obs.rt = rt.toDouble();
      QVariant mz = query.value("mz");
      if (!mz.isNull()) obs.mz = mz.toDouble();
      Key id = query.value("id").toLongLong();
      if (have_meta_info) handleQueryMetaInfo_(subquery_info, obs, id);
      ID::ObservationRef ref = id_data.registerObservation(obs);
      observation_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadParentSequences_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ParentSequence")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentSequence"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentSequence");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step, "ID_ParentSequence");

    while (query.next())
    {
      String accession = query.value("accession").toString();
      ID::ParentSequence parent(accession);
      int molecule_type_index = query.value("molecule_type_id").toInt() - 1;
      parent.molecule_type = ID::MoleculeType(molecule_type_index);
      parent.sequence = query.value("sequence").toString();
      parent.description = query.value("description").toString();
      parent.coverage = query.value("coverage").toDouble();
      parent.is_decoy = query.value("is_decoy").toInt();
      Key id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, parent, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, parent, id);
      }
      ID::ParentSequenceRef ref = id_data.registerParentSequence(parent);
      parent_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadParentGroupSets_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ParentGroupSet")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentGroupSet ORDER BY grouping_order ASC"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentGroupSet");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_ParentGroupSet");

    QSqlQuery subquery_group(db);
    subquery_group.setForwardOnly(true);
    subquery_group.prepare("SELECT * FROM ID_ParentGroup WHERE grouping_id = :id");

    QSqlQuery subquery_parent(db);
    subquery_parent.setForwardOnly(true);
    subquery_parent.prepare(
      "SELECT parent_id FROM ID_ParentGroup_ParentSequence WHERE group_id = :id");

    while (query.next())
    {
      ID::ParentGroupSet grouping(query.value("label").toString());
      Key grouping_id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, grouping, grouping_id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, grouping, grouping_id);
      }

      subquery_group.bindValue(":id", grouping_id);
      if (!subquery_group.exec())
      {
        raiseDBError_(subquery_group.lastError(), __LINE__,
                      OPENMS_PRETTY_FUNCTION, "error reading from database");
      }

      // get all groups in this grouping:
      map<Key, ID::ParentGroup> groups_map;
      while (subquery_group.next())
      {
        Key group_id = subquery_group.value("id").toLongLong();
        QVariant score_type_id = subquery_group.value("score_type_id");
        if (score_type_id.isNull()) // no scores
        {
          groups_map[group_id]; // insert empty group
        }
        else
        {
          ID::ScoreTypeRef ref = score_type_refs_[score_type_id.toLongLong()];
          groups_map[group_id].scores[ref] =
            subquery_group.value("score").toDouble();
        }
      }
      // get parent sequences in each group:
      for (auto& pair : groups_map)
      {
        subquery_parent.bindValue(":id", pair.first);
        if (!subquery_parent.exec())
        {
          raiseDBError_(subquery_parent.lastError(), __LINE__,
                        OPENMS_PRETTY_FUNCTION, "error reading from database");
        }
        while (subquery_parent.next())
        {
          Key parent_id = subquery_parent.value(0).toLongLong();
          pair.second.parent_refs.insert(
            parent_refs_[parent_id]);
        }
        grouping.groups.insert(pair.second);
      }

      id_data.registerParentGroupSet(grouping);
    }
  }


  void OMSFileLoad::loadIdentifiedCompounds_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_IdentifiedCompound")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    QString sql_select =
      "SELECT * FROM ID_IdentifiedMolecule JOIN ID_IdentifiedCompound " \
      "ON ID_IdentifiedMolecule.id = ID_IdentifiedCompound.molecule_id";
    if (!query.exec(sql_select))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "ID_IdentifiedMolecule");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step, "ID_IdentifiedMolecule");

    while (query.next())
    {
      ID::IdentifiedCompound compound(
        query.value("identifier").toString(),
        EmpiricalFormula(query.value("formula").toString()),
        query.value("name").toString(),
        query.value("smile").toString(),
        query.value("inchi").toString());
      Key id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, compound, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, compound, id);
      }
      ID::IdentifiedCompoundRef ref = id_data.registerIdentifiedCompound(compound);
      identified_molecule_vars_[id] = ref;
    }
  }


  void OMSFileLoad::handleQueryParentMatch_(QSqlQuery& query,
                                            IdentificationData::ParentMatches& parent_matches,
                                            Key molecule_id)
  {
    query.bindValue(":id", molecule_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::ParentSequenceRef ref =
        parent_refs_[query.value("parent_id").toLongLong()];
      ID::ParentMatch match;
      QVariant start_pos = query.value("start_pos");
      QVariant end_pos = query.value("end_pos");
      if (!start_pos.isNull()) match.start_pos = start_pos.toInt();
      if (!end_pos.isNull()) match.end_pos = end_pos.toInt();
      match.left_neighbor = query.value("left_neighbor").toString();
      match.right_neighbor = query.value("right_neighbor").toString();
      parent_matches[ref].insert(match);
    }
  }


  void OMSFileLoad::loadIdentifiedSequences_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_IdentifiedMolecule")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    query.prepare("SELECT * FROM ID_IdentifiedMolecule "          \
                  "WHERE molecule_type_id = :molecule_type_id");
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_IdentifiedMolecule");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_IdentifiedMolecule");
    QSqlQuery subquery_parent(db);
    bool have_parent_matches = tableExists_(db_name_,
                                            "ID_ParentMatch");
    if (have_parent_matches)
    {
      subquery_parent.setForwardOnly(true);
      subquery_parent.prepare("SELECT * FROM ID_ParentMatch " \
                              "WHERE molecule_id = :id");
    }

    // load peptides:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      String sequence = query.value("identifier").toString();
      ID::IdentifiedPeptide peptide(AASequence::fromString(sequence));
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, peptide, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, peptide, id);
      }
      if (have_parent_matches)
      {
        handleQueryParentMatch_(subquery_parent, peptide.parent_matches, id);
      }
      ID::IdentifiedPeptideRef ref = id_data.registerIdentifiedPeptide(peptide);
      identified_molecule_vars_[id] = ref;
    }

    // load RNA oligos:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      String sequence = query.value("identifier").toString();
      ID::IdentifiedOligo oligo(NASequence::fromString(sequence));
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, oligo, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, oligo, id);
      }
      if (have_parent_matches)
      {
        handleQueryParentMatch_(subquery_parent, oligo.parent_matches, id);
      }
      ID::IdentifiedOligoRef ref = id_data.registerIdentifiedOligo(oligo);
      identified_molecule_vars_[id] = ref;
    }
  }


  void OMSFileLoad::handleQueryPeakAnnotation_(QSqlQuery& query,
                                               ID::ObservationMatch& match,
                                               Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      QVariant processing_step_id = query.value("processing_step_id");
      std::optional<ID::ProcessingStepRef> processing_step_opt = std::nullopt;
      if (!processing_step_id.isNull())
      {
        processing_step_opt =
          processing_step_refs_[processing_step_id.toLongLong()];
      }
      PeptideHit::PeakAnnotation ann;
      ann.annotation = query.value("peak_annotation").toString();
      ann.charge = query.value("peak_charge").toInt();
      ann.mz = query.value("peak_mz").toDouble();
      ann.intensity = query.value("peak_intensity").toDouble();
      match.peak_annotations[processing_step_opt].push_back(ann);
    }
  }


  void OMSFileLoad::loadAdducts_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "AdductInfo")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM AdductInfo"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }

    while (query.next())
    {
      EmpiricalFormula formula(query.value("formula").toString());
      AdductInfo adduct(query.value("name").toString(), formula,
                        query.value("charge").toInt(),
                        query.value("mol_multiplier").toInt());
      ID::AdductRef ref = id_data.registerAdduct(adduct);
      adduct_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFileLoad::loadObservationMatches_(IdentificationData& id_data)
  {
    if (!tableExists_(db_name_, "ID_ObservationMatch")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ObservationMatch"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ObservationMatch");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_ObservationMatch");
    QSqlQuery subquery_ann(db);
    bool have_peak_annotations =
      tableExists_(db_name_, "ID_ObservationMatch_PeakAnnotation");
    if (have_peak_annotations)
    {
      subquery_ann.setForwardOnly(true);
      subquery_ann.prepare(
        "SELECT * FROM ID_ObservationMatch_PeakAnnotation " \
        "WHERE parent_id = :id");
    }

    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      Key molecule_id = query.value("identified_molecule_id").toLongLong();
      Key query_id = query.value("observation_id").toLongLong();
      ID::ObservationMatch match(identified_molecule_vars_[molecule_id],
                                   observation_refs_[query_id],
                                   query.value("charge").toInt());
      QVariant adduct_id = query.value("adduct_id"); // adduct is optional
      if (!adduct_id.isNull())
      {
        match.adduct_opt = adduct_refs_[adduct_id.toLongLong()];
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, match, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, match, id);
      }
      if (have_peak_annotations)
      {
        handleQueryPeakAnnotation_(subquery_ann, match, id);
      }
      ID::ObservationMatchRef ref = id_data.registerObservationMatch(match);
      observation_match_refs_[id] = ref;
    }
  }


  void OMSFileLoad::load(IdentificationData& id_data)
  {
    startProgress(0, 12, "Reading identification data from file");
    loadInputFiles_(id_data);
    nextProgress();
    loadScoreTypes_(id_data);
    nextProgress();
    loadProcessingSoftwares_(id_data);
    nextProgress();
    loadDBSearchParams_(id_data);
    nextProgress();
    loadProcessingSteps_(id_data);
    nextProgress();
    loadObservations_(id_data);
    nextProgress();
    loadParentSequences_(id_data);
    nextProgress();
    loadParentGroupSets_(id_data);
    nextProgress();
    loadIdentifiedCompounds_(id_data);
    nextProgress();
    loadIdentifiedSequences_(id_data);
    nextProgress();
    loadAdducts_(id_data);
    nextProgress();
    loadObservationMatches_(id_data);
    endProgress();
    // @TODO: load input match groups
  }


  void OMSFileLoad::loadMapMetaData_(FeatureMap& features)
  {
    if (!tableExists_(db_name_, "FEAT_MapMetaData")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM FEAT_MapMetaData"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }

    query.next(); // there should be only one row
    Key id = query.value("unique_id").toLongLong();
    features.setUniqueId(id);
    features.setIdentifier(query.value("identifier").toString());
    features.setLoadedFilePath(query.value("file_path").toString());
    String file_type = query.value("file_type").toString();
    features.setLoadedFilePath(FileTypes::nameToType(file_type));
    QSqlQuery query_meta(QSqlDatabase::database(db_name_));
    if (prepareQueryMetaInfo_(query_meta, "FEAT_MapMetaData"))
    {
      handleQueryMetaInfo_(query_meta, features, id);
    }
  }


  void OMSFileLoad::loadDataProcessing_(FeatureMap& features)
  {
    if (!tableExists_(db_name_, "FEAT_DataProcessing")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM FEAT_DataProcessing ORDER BY position ASC"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }

    QSqlQuery subquery_info(QSqlDatabase::database(db_name_));
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "FEAT_DataProcessing");

    while (query.next())
    {
      DataProcessing proc;
      Software sw(query.value("software_name").toString(),
                  query.value("software_version").toString());
      proc.setSoftware(sw);
      vector<String> actions =
        ListUtils::create<String>(query.value("processing_actions").toString());
      for (const String& action : actions)
      {
        auto pos = find(begin(DataProcessing::NamesOfProcessingAction),
                          end(DataProcessing::NamesOfProcessingAction), action);
        if (pos != end(DataProcessing::NamesOfProcessingAction))
        {
          Size index = pos - begin(DataProcessing::NamesOfProcessingAction);
          proc.getProcessingActions().insert(DataProcessing::ProcessingAction(index));
        }
        else // @TODO: throw an exception here?
        {
          OPENMS_LOG_ERROR << "Error: unknown data processing action '" << action << "' - skipping";
        }
      }
      DateTime time;
      time.set(query.value("completion_time").toString());
      proc.setCompletionTime(time);
      if (have_meta_info)
      {
        Key id = query.value("id").toLongLong();
        handleQueryMetaInfo_(subquery_info, proc, id);
      }
      features.getDataProcessing().push_back(proc);
    }
  }


  Feature OMSFileLoad::loadFeatureAndSubordinates_(
    QSqlQuery& query_feat, std::optional<QSqlQuery>& query_meta,
    std::optional<QSqlQuery>& query_hull, std::optional<QSqlQuery>& query_match)
  {
    Feature feature;
    int id = query_feat.value("id").toInt();
    feature.setRT(query_feat.value("rt").toDouble());
    feature.setMZ(query_feat.value("mz").toDouble());
    feature.setIntensity(query_feat.value("intensity").toDouble());
    feature.setCharge(query_feat.value("charge").toInt());
    feature.setWidth(query_feat.value("width").toDouble());
    feature.setOverallQuality(query_feat.value("overall_quality").toDouble());
    feature.setQuality(0, query_feat.value("rt_quality").toDouble());
    feature.setQuality(1, query_feat.value("mz_quality").toDouble());
    feature.setUniqueId(query_feat.value("unique_id").toLongLong());
    QVariant primary_id = query_feat.value("primary_molecule_id"); // optional
    if (!primary_id.isNull())
    {
      feature.setPrimaryID(identified_molecule_vars_[primary_id.toLongLong()]);
    }
    // meta data:
    if (query_meta)
    {
      handleQueryMetaInfo_(*query_meta, feature, id);
    }
    // convex hulls:
    if (query_hull)
    {
      query_hull->bindValue(":id", id);
      if (!query_hull->exec())
      {
        raiseDBError_(query_hull->lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error reading from database");
      }
      while (query_hull->next())
      {
        Size hull_index = query_hull->value("hull_index").toUInt();
        // first row should have max. hull index (sorted descending):
        if (feature.getConvexHulls().size() <= hull_index)
        {
          feature.getConvexHulls().resize(hull_index + 1);
        }
        ConvexHull2D::PointType point(query_hull->value("point_x").toDouble(),
                                      query_hull->value("point_y").toDouble());
        // @TODO: this may be inefficient (see implementation of "addPoint"):
        feature.getConvexHulls()[hull_index].addPoint(point);
      }
    }
    // ID matches:
    if (query_match)
    {
      query_match->bindValue(":id", id);
      if (!query_match->exec())
      {
        raiseDBError_(query_match->lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error reading from database");
      }
      while (query_match->next())
      {
        Key match_id = query_match->value("observation_match_id").toLongLong();
        feature.addIDMatch(observation_match_refs_[match_id]);
      }
    }
    // subordinates:
    QSqlQuery query_sub(QSqlDatabase::database(db_name_));
    query_sub.setForwardOnly(true);
    QString sql = "SELECT * FROM FEAT_Feature WHERE subordinate_of = " +
      QString::number(id) + " ORDER BY id ASC";
    if (!query_sub.exec(sql))
    {
      raiseDBError_(query_sub.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query_sub.next())
    {
      Feature sub = loadFeatureAndSubordinates_(query_sub, query_meta,
                                                query_hull, query_match);
      feature.getSubordinates().push_back(sub);
    }
    return feature;
  }


  void OMSFileLoad::loadFeatures_(FeatureMap& features)
  {
    if (!tableExists_(db_name_, "FEAT_Feature")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);

    // start with top-level features only:
    QSqlQuery query_feat(db);
    query_feat.setForwardOnly(true);
    if (!query_feat.exec("SELECT * FROM FEAT_Feature WHERE subordinate_of IS NULL ORDER BY id ASC"))
    {
      raiseDBError_(query_feat.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // prepare sub-queries (optional - corresponding tables may not be present):
    std::optional<QSqlQuery> query_meta(db);
    if (!prepareQueryMetaInfo_(*query_meta, "FEAT_Feature"))
    {
      query_meta = std::nullopt;
    }
    std::optional<QSqlQuery> query_hull;
    if (tableExists_(db_name_, "FEAT_ConvexHull"))
    {
      query_hull = QSqlQuery(db);
      query_hull->prepare("SELECT * FROM FEAT_ConvexHull WHERE feature_id = :id " \
                         "ORDER BY hull_index DESC, point_index ASC");
    }
    std::optional<QSqlQuery> query_match;
    if (tableExists_(db_name_, "FEAT_ObservationMatch"))
    {
      query_match = QSqlQuery(db);
      query_match->prepare("SELECT * FROM FEAT_ObservationMatch WHERE feature_id = :id");
    }

    while (query_feat.next())
    {
      Feature feature = loadFeatureAndSubordinates_(query_feat, query_meta,
                                                    query_hull, query_match);
      features.push_back(feature);
    }
  }


  void OMSFileLoad::load(FeatureMap& features)
  {
    load(features.getIdentificationData()); // load IDs, if any
    startProgress(0, 3, "Reading feature data from file");
    loadMapMetaData_(features);
    nextProgress();
    loadDataProcessing_(features);
    nextProgress();
    loadFeatures_(features);
    endProgress();
  }
}
