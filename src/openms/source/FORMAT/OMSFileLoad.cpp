// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFileLoad.h>
#include <OpenMS/FORMAT/OMSFileStore.h> // for "raiseDBError_"
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <QtCore/QString>
// JSON export:
#include <QtCore/QJsonDocument>
#include <QtCore/QJsonObject>

#include <SQLiteCpp/Database.h>

#include <sqlite3.h>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
  // initialize lookup table:
  map<QString, QString> OMSFileLoad::export_order_by_ = {
    {"version", ""},
    {"ID_IdentifiedCompound", "molecule_id"},
    {"ID_ParentMatch", "molecule_id, parent_id, start_pos, end_pos"},
    {"ID_ParentGroup_ParentSequence", "group_id, parent_id"},
    {"ID_ProcessingStep_InputFile", "processing_step_id, input_file_id"},
    {"ID_ProcessingSoftware_AssignedScore", "software_id, score_type_order"},
    {"ID_ObservationMatch_PeakAnnotation", "parent_id, processing_step_id, peak_mz, peak_annotation"},
    {"FEAT_ConvexHull", "feature_id, hull_index, point_index"},
    {"FEAT_ObservationMatch", "feature_id, observation_match_id"},
    {"FEAT_MapMetaData", "unique_id"}
  };


  OMSFileLoad::OMSFileLoad(const String& filename, LogType log_type):
    db_(make_unique<SQLite::Database>(filename))
  {
    setLogType(log_type);

    // read version number:
    try
    {
      auto version = db_->execAndGet("SELECT OMSFile FROM version");
      version_number_ = version.getInt();
    }
    catch (...)
    {
      raiseDBError_(db_->getErrorMsg(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading file format version number");
    }
  }


  OMSFileLoad::~OMSFileLoad()
  {
  }


  bool OMSFileLoad::isEmpty_(const SQLite::Statement& query)
  {
    return query.getQuery().empty();
  }


  // currently not needed:
  // CVTerm OMSFileLoad::loadCVTerm_(int id)
  // {
  //   // this assumes that the "CVTerm" table exists!
  //   SQLite::Statement query(db_);
  //
  //   QString sql_select = "SELECT * FROM CVTerm WHERE id = " + QString(id);
  //   if (!query.exec(sql_select) || !query.executeStep())
  //   {
  //     raiseDBError_(model.getErrorMsg(), __LINE__, OPENMS_PRETTY_FUNCTION,
  //                   "error reading from database");
  //   }
  //   return CVTerm(query.getColumn("accession").getString(),
  //                 query.getColumn("name").getString(),
  //                 query.getColumn("cv_identifier_ref").getString());
  // }


  void OMSFileLoad::loadScoreTypes_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ScoreType")) return;
    if (!db_->tableExists("CVTerm")) // every score type is a CV term
    {
      String msg = "required database table 'CVTerm' not found";
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    // careful - both joined tables have an "id" field, need to exclude one:
    SQLite::Statement query(*db_, "SELECT S.*, C.accession, C.name, C.cv_identifier_ref " \
                    "FROM ID_ScoreType AS S JOIN CVTerm AS C "          \
                    "ON S.cv_term_id = C.id");
    while (query.executeStep())
    {
      CVTerm cv_term(query.getColumn("accession").getString(),
                     query.getColumn("name").getString(),
                     query.getColumn("cv_identifier_ref").getString());
      bool higher_better = query.getColumn("higher_better").getInt();
      ID::ScoreType score_type(cv_term, higher_better);
      ID::ScoreTypeRef ref = id_data.registerScoreType(score_type);
      score_type_refs_[query.getColumn("id").getInt64()] = ref;
    }
  }


  void OMSFileLoad::loadInputFiles_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_InputFile")) return;

    SQLite::Statement query(*db_, "SELECT * FROM ID_InputFile");
    while (query.executeStep())
    {
      ID::InputFile input(query.getColumn("name").getString(),
                          query.getColumn("experimental_design_id").getString());
      String primary_files = query.getColumn("primary_files").getString();
      vector<String> pf_list = ListUtils::create<String>(primary_files);
      input.primary_files.insert(pf_list.begin(), pf_list.end());
      ID::InputFileRef ref = id_data.registerInputFile(input);
      input_file_refs_[query.getColumn("id").getInt64()] = ref;
    }
  }


  void OMSFileLoad::loadProcessingSoftwares_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ProcessingSoftware")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_ProcessingSoftware");
    bool have_scores = db_->tableExists("ID_ProcessingSoftware_AssignedScore");
    SQLite::Statement subquery(*db_, "");
    if (have_scores)
    {
      subquery = SQLite::Statement(*db_, "SELECT score_type_id "                         \
                       "FROM ID_ProcessingSoftware_AssignedScore " \
                       "WHERE software_id = :id ORDER BY score_type_order ASC");
    }
    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      ID::ProcessingSoftware software(query.getColumn("name").getString(),
                                          query.getColumn("version").getString());
      if (have_scores)
      {
        subquery.bind(":id", id);
        while (subquery.executeStep())
        {
          Key score_type_id = subquery.getColumn(0).getInt64();
          software.assigned_scores.push_back(score_type_refs_[score_type_id]);
        }
        subquery.reset(); // get ready for new executeStep()
      }
      ID::ProcessingSoftwareRef ref = id_data.registerProcessingSoftware(software);
      processing_software_refs_[id] = ref;
    }
  }


  DataValue OMSFileLoad::makeDataValue_(const SQLite::Statement& query)
  {
    DataValue::DataType type = DataValue::EMPTY_VALUE;
    int type_index = query.getColumn("data_type_id").getInt();
    if (type_index > 0) type = DataValue::DataType(type_index - 1);
    String value = query.getColumn("value").getString();
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


  bool OMSFileLoad::prepareQueryMetaInfo_(SQLite::Statement& query,
                                          const String& parent_table)
  {
    String table_name = parent_table + "_MetaInfo";
    if (!db_->tableExists(table_name)) return false;

    String sql_select =
    "SELECT * FROM " + table_name.toQString() + " AS MI " \
    "WHERE MI.parent_id = :id";

    if (version_number_ < 4)
    {
      sql_select =
      "SELECT * FROM " + table_name.toQString() + " AS MI " \
      "JOIN DataValue AS DV ON MI.data_value_id = DV.id "   \
      "WHERE MI.parent_id = :id";
    }
    query = SQLite::Statement(*db_, sql_select);
    return true;
  }


  bool OMSFileLoad::prepareQueryAppliedProcessingStep_(SQLite::Statement& query,
                                                       const String& parent_table)
  {
    String table_name = parent_table + "_AppliedProcessingStep";
    if (!db_->tableExists(table_name)) return false;

    //
    String sql_select = "SELECT * FROM " + table_name.toQString() +
      " WHERE parent_id = :id ORDER BY processing_step_order ASC";
    query = SQLite::Statement(*db_, sql_select);
    return true;
  }


  void OMSFileLoad::handleQueryMetaInfo_(SQLite::Statement& query,
                                         MetaInfoInterface& info,
                                         Key parent_id)
  {
    query.bind(":id", parent_id);
    while (query.executeStep())
    {
      DataValue value = makeDataValue_(query);
      info.setMetaValue(query.getColumn("name").getString(), value);
    }
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::handleQueryAppliedProcessingStep_(
    SQLite::Statement& query,
    IdentificationDataInternal::ScoredProcessingResult& result,
    Key parent_id)
  {
    query.bind(":id", parent_id);
    while (query.executeStep())
    {
      ID::AppliedProcessingStep step;
      auto step_id_opt = query.getColumn("processing_step_id");
      if (!step_id_opt.isNull())
      {
        step.processing_step_opt =
          processing_step_refs_[step_id_opt.getInt64()];
      }
      auto score_type_opt = query.getColumn("score_type_id");
      if (!score_type_opt.isNull())
      {
        step.scores[score_type_refs_[score_type_opt.getInt64()]] =
          query.getColumn("score").getDouble();
      }
      result.addProcessingStep(step); // this takes care of merging the steps
    }
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::loadDBSearchParams_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_DBSearchParam")) return;

    SQLite::Statement query(*db_, "SELECT * FROM ID_DBSearchParam");
    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      ID::DBSearchParam param;
      int molecule_type_index = query.getColumn("molecule_type_id").getInt() - 1;
      param.molecule_type = ID::MoleculeType(molecule_type_index);
      int mass_type_index = query.getColumn("mass_type_average").getInt();
      param.mass_type = ID::MassType(mass_type_index);
      param.database = query.getColumn("database").getString();
      param.database_version = query.getColumn("database_version").getString();
      param.taxonomy = query.getColumn("taxonomy").getString();
      vector<Int> charges =
        ListUtils::create<Int>(query.getColumn("charges").getString());
      param.charges.insert(charges.begin(), charges.end());
      vector<String> fixed_mods =
        ListUtils::create<String>(query.getColumn("fixed_mods").getString());
      param.fixed_mods.insert(fixed_mods.begin(), fixed_mods.end());
      vector<String> variable_mods =
        ListUtils::create<String>(query.getColumn("variable_mods").getString());
      param.variable_mods.insert(variable_mods.begin(), variable_mods.end());
      param.precursor_mass_tolerance =
        query.getColumn("precursor_mass_tolerance").getDouble();
      param.fragment_mass_tolerance =
        query.getColumn("fragment_mass_tolerance").getDouble();
      param.precursor_tolerance_ppm =
        query.getColumn("precursor_tolerance_ppm").getInt();
      param.fragment_tolerance_ppm =
        query.getColumn("fragment_tolerance_ppm").getInt();
      String enzyme = query.getColumn("digestion_enzyme").getString();
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
      if (version_number_ > 1)
      {
        String spec = query.getColumn("enzyme_term_specificity").getString();
        param.enzyme_term_specificity = EnzymaticDigestion::getSpecificityByName(spec);
      }
      param.missed_cleavages = query.getColumn("missed_cleavages").getUInt();
      param.min_length = query.getColumn("min_length").getUInt();
      param.max_length = query.getColumn("max_length").getUInt();
      ID::SearchParamRef ref = id_data.registerDBSearchParam(param);
      search_param_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadProcessingSteps_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ProcessingStep")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_ProcessingStep");
    SQLite::Statement subquery_file(*db_, "");
    bool have_input_files = db_->tableExists(
                                         "ID_ProcessingStep_InputFile");
    if (have_input_files)
    {
      subquery_file = SQLite::Statement(*db_, "SELECT input_file_id "                 \
                            "FROM ID_ProcessingStep_InputFile " \
                            "WHERE processing_step_id = :id");
    }
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "ID_ProcessingStep");
    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      Key software_id = query.getColumn("software_id").getInt64();
      ID::ProcessingStep step(processing_software_refs_[software_id]);
      String date_time = query.getColumn("date_time").getString();
      if (!date_time.empty()) step.date_time.set(date_time);
      if (have_input_files)
      {
        subquery_file.bind(":id", id);
        while (subquery_file.executeStep())
        {
          Key input_file_id = subquery_file.getColumn(0).getInt64();
          // the foreign key constraint should ensure that look-up succeeds:
          step.input_file_refs.push_back(input_file_refs_[input_file_id]);
        }
        subquery_file.reset(); // get ready for new executeStep()
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, step, id);
      }
      ID::ProcessingStepRef ref;
      auto opt_search_param_id = query.getColumn("search_param_id");
      if (opt_search_param_id.isNull()) // no DB search params available
      {
        ref = id_data.registerProcessingStep(step);
      }
      else
      {
        ID::SearchParamRef search_param_ref =
          search_param_refs_[opt_search_param_id.getInt64()];
        ref = id_data.registerProcessingStep(step, search_param_ref);
      }
      processing_step_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadObservations_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_Observation")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_Observation");
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_Observation");

    while (query.executeStep())
    {
      auto input_file_id = query.getColumn("input_file_id");
      ID::Observation obs(query.getColumn("data_id").getString(),
                          input_file_refs_[input_file_id.getInt64()]);
      auto rt = query.getColumn("rt");
      if (!rt.isNull()) obs.rt = rt.getDouble();
      auto mz = query.getColumn("mz");
      if (!mz.isNull()) obs.mz = mz.getDouble();
      Key id = query.getColumn("id").getInt64();
      if (have_meta_info) handleQueryMetaInfo_(subquery_info, obs, id);
      ID::ObservationRef ref = id_data.registerObservation(obs);
      observation_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadParentSequences_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ParentSequence")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_ParentSequence");
    // @TODO: can we combine handling of meta info and applied processing steps?
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentSequence");
    SQLite::Statement subquery_step(*db_, "");
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step, "ID_ParentSequence");

    while (query.executeStep())
    {
      String accession = query.getColumn("accession").getString();
      ID::ParentSequence parent(accession);
      int molecule_type_index = query.getColumn("molecule_type_id").getInt() - 1;
      parent.molecule_type = ID::MoleculeType(molecule_type_index);
      parent.sequence = query.getColumn("sequence").getString();
      parent.description = query.getColumn("description").getString();
      parent.coverage = query.getColumn("coverage").getDouble();
      parent.is_decoy = query.getColumn("is_decoy").getInt();
      Key id = query.getColumn("id").getInt64();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, parent, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, parent, id);
      }
      ID::ParentSequenceRef ref = id_data.registerParentSequence(parent);
      parent_sequence_refs_[id] = ref;
    }
  }


  void OMSFileLoad::loadParentGroupSets_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ParentGroupSet")) return;

    // "grouping_order" column was removed in schema version 3:
    String order_by = version_number_ > 2 ? "id" : "grouping_order";

    SQLite::Statement query(*db_, "SELECT * FROM ID_ParentGroupSet ORDER BY " + order_by + " ASC");
    // @TODO: can we combine handling of meta info and applied processing steps?
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentGroupSet");
    SQLite::Statement subquery_step(*db_, "");
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_ParentGroupSet");

    SQLite::Statement subquery_group(*db_, "SELECT * FROM ID_ParentGroup WHERE grouping_id = :id");

    SQLite::Statement subquery_parent(*db_, "SELECT parent_id FROM ID_ParentGroup_ParentSequence WHERE group_id = :id");

    while (query.executeStep())
    {
      ID::ParentGroupSet grouping(query.getColumn("label").getString());
      Key grouping_id = query.getColumn("id").getInt64();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, grouping, grouping_id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, grouping, grouping_id);
      }

      subquery_group.bind(":id", grouping_id);
      // get all groups in this grouping:
      map<Key, ID::ParentGroup> groups_map;
      while (subquery_group.executeStep())
      {
        Key group_id = subquery_group.getColumn("id").getInt64();
        auto score_type_id = subquery_group.getColumn("score_type_id");
        if (score_type_id.isNull()) // no scores
        {
          groups_map[group_id]; // insert empty group
        }
        else
        {
          ID::ScoreTypeRef ref = score_type_refs_[score_type_id.getInt64()];
          groups_map[group_id].scores[ref] =
            subquery_group.getColumn("score").getDouble();
        }
      }
      subquery_group.reset(); // get ready for new executeStep()
      // get parent sequences in each group:
      for (auto& pair : groups_map)
      {
        subquery_parent.bind(":id", pair.first);
        while (subquery_parent.executeStep())
        {
          Key parent_id = subquery_parent.getColumn(0).getInt64();
          pair.second.parent_refs.insert(
            parent_sequence_refs_[parent_id]);
        }
        subquery_parent.reset(); // get ready for new executeStep()
        grouping.groups.insert(pair.second);
      }

      id_data.registerParentGroupSet(grouping);
    }
  }


  void OMSFileLoad::loadIdentifiedCompounds_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_IdentifiedCompound")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_IdentifiedMolecule JOIN ID_IdentifiedCompound " \
      "ON ID_IdentifiedMolecule.id = ID_IdentifiedCompound.molecule_id");
    // @TODO: can we combine handling of meta info and applied processing steps?
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "ID_IdentifiedMolecule");
    SQLite::Statement subquery_step(*db_, "");
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step, "ID_IdentifiedMolecule");

    while (query.executeStep())
    {
      ID::IdentifiedCompound compound(
        query.getColumn("identifier").getString(),
        EmpiricalFormula(query.getColumn("formula").getString()),
        query.getColumn("name").getString(),
        query.getColumn("smile").getString(),
        query.getColumn("inchi").getString());
      Key id = query.getColumn("id").getInt64();
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


  void OMSFileLoad::handleQueryParentMatch_(SQLite::Statement& query,
                                            IdentificationData::ParentMatches& parent_matches,
                                            Key molecule_id)
  {
    query.bind(":id", molecule_id);
    while (query.executeStep())
    {
      ID::ParentSequenceRef ref = parent_sequence_refs_[query.getColumn("parent_id").getInt64()];
      ID::ParentMatch match;
      auto start_pos = query.getColumn("start_pos");
      auto end_pos = query.getColumn("end_pos");
      if (!start_pos.isNull()) match.start_pos = start_pos.getInt();
      if (!end_pos.isNull()) match.end_pos = end_pos.getInt();
      match.left_neighbor = query.getColumn("left_neighbor").getString();
      match.right_neighbor = query.getColumn("right_neighbor").getString();
      parent_matches[ref].insert(match);
    }
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::loadIdentifiedSequences_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_IdentifiedMolecule")) return;

    SQLite::Statement query(*db_, "SELECT * FROM ID_IdentifiedMolecule "          \
                            "WHERE molecule_type_id = :molecule_type_id");
    // @TODO: can we combine handling of meta info and applied processing steps?
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_IdentifiedMolecule");
    SQLite::Statement subquery_step(*db_, "");
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_IdentifiedMolecule");
    SQLite::Statement subquery_parent(*db_, "");
    bool have_parent_matches = db_->tableExists(
                                            "ID_ParentMatch");
    if (have_parent_matches)
    {
      subquery_parent = SQLite::Statement(*db_, "SELECT * FROM ID_ParentMatch " \
                              "WHERE molecule_id = :id");
    }

    // load peptides:
    query.bind(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      String sequence = query.getColumn("identifier").getString();
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
    query.reset(); // get ready for new executeStep()

    // load RNA oligos:
    query.bind(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      String sequence = query.getColumn("identifier").getString();
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
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::handleQueryPeakAnnotation_(SQLite::Statement& query,
                                               ID::ObservationMatch& match,
                                               Key parent_id)
  {
    query.bind(":id", parent_id);
    while (query.executeStep())
    {
      auto processing_step_id = query.getColumn("processing_step_id");
      std::optional<ID::ProcessingStepRef> processing_step_opt = std::nullopt;
      if (!processing_step_id.isNull())
      {
        processing_step_opt =
          processing_step_refs_[processing_step_id.getInt64()];
      }
      PeptideHit::PeakAnnotation ann;
      ann.annotation = query.getColumn("peak_annotation").getString();
      ann.charge = query.getColumn("peak_charge").getInt();
      ann.mz = query.getColumn("peak_mz").getDouble();
      ann.intensity = query.getColumn("peak_intensity").getDouble();
      match.peak_annotations[processing_step_opt].push_back(ann);
    }
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::loadAdducts_(IdentificationData& id_data)
  {
    if (!db_->tableExists("AdductInfo")) return;

    SQLite::Statement query(*db_, "SELECT * FROM AdductInfo");
    while (query.executeStep())
    {
      EmpiricalFormula formula(query.getColumn("formula").getString());
      AdductInfo adduct(query.getColumn("name").getString(), formula,
                        query.getColumn("charge").getInt(),
                        query.getColumn("mol_multiplier").getInt());
      ID::AdductRef ref = id_data.registerAdduct(adduct);
      adduct_refs_[query.getColumn("id").getInt64()] = ref;
    }
  }


  void OMSFileLoad::loadObservationMatches_(IdentificationData& id_data)
  {
    if (!db_->tableExists("ID_ObservationMatch")) return;


    SQLite::Statement query(*db_, "SELECT * FROM ID_ObservationMatch");
    // @TODO: can we combine handling of meta info and applied processing steps?
    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ObservationMatch");
    SQLite::Statement subquery_step(*db_, "");
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_ObservationMatch");
    SQLite::Statement subquery_ann(*db_, "");
    bool have_peak_annotations =
      db_->tableExists("ID_ObservationMatch_PeakAnnotation");
    if (have_peak_annotations)
    {
      subquery_ann = SQLite::Statement(*db_,
        "SELECT * FROM ID_ObservationMatch_PeakAnnotation " \
        "WHERE parent_id = :id");
    }

    while (query.executeStep())
    {
      Key id = query.getColumn("id").getInt64();
      Key molecule_id = query.getColumn("identified_molecule_id").getInt64();
      Key query_id = query.getColumn("observation_id").getInt64();
      ID::ObservationMatch match(identified_molecule_vars_[molecule_id],
                                   observation_refs_[query_id],
                                   query.getColumn("charge").getInt());
      auto adduct_id = query.getColumn("adduct_id"); // adduct is optional
      if (!adduct_id.isNull())
      {
        match.adduct_opt = adduct_refs_[adduct_id.getInt64()];
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

  template <class MapType>
  String OMSFileLoad::loadMapMetaDataTemplate_(MapType& features)
  {
    if (!db_->tableExists("FEAT_MapMetaData")) return "";

    SQLite::Statement query(*db_, "SELECT * FROM FEAT_MapMetaData");
    query.executeStep(); // there should be only one row
    Key id = query.getColumn("unique_id").getInt64();
    features.setUniqueId(id);
    features.setIdentifier(query.getColumn("identifier").getString());
    features.setLoadedFilePath(query.getColumn("file_path").getString());
    String file_type = query.getColumn("file_type").getString();
    features.setLoadedFilePath(FileTypes::nameToType(file_type));
    SQLite::Statement query_meta(*db_, "");
    if (prepareQueryMetaInfo_(query_meta, "FEAT_MapMetaData"))
    {
      handleQueryMetaInfo_(query_meta, features, id);
    }
    if (version_number_ < 5) return ""; // "experiment_type" column doesn't exist yet
    return query.getColumn("experiment_type").getString(); // for consensus map only
  }

  // template specializations:
  template String OMSFileLoad::loadMapMetaDataTemplate_<FeatureMap>(FeatureMap&);
  template String OMSFileLoad::loadMapMetaDataTemplate_<ConsensusMap>(ConsensusMap&);


  void OMSFileLoad::loadMapMetaData_(FeatureMap& features)
  {
    loadMapMetaDataTemplate_(features);
  }

  void OMSFileLoad::loadMapMetaData_(ConsensusMap& consensus)
  {
    String experiment_type = loadMapMetaDataTemplate_(consensus);
    consensus.setExperimentType(experiment_type);
  }


  void OMSFileLoad::loadDataProcessing_(vector<DataProcessing>& data_processing)
  {
    if (!db_->tableExists("FEAT_DataProcessing")) return;

    // "position" column was removed in schema version 3:
    String order_by = version_number_ > 2 ? "id" : "position";
    SQLite::Statement query(*db_, "SELECT * FROM FEAT_DataProcessing ORDER BY " + order_by + " ASC");

    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "FEAT_DataProcessing");

    while (query.executeStep())
    {
      DataProcessing proc;
      Software sw(query.getColumn("software_name").getString(),
                  query.getColumn("software_version").getString());
      proc.setSoftware(sw);
      vector<String> actions =
        ListUtils::create<String>(query.getColumn("processing_actions").getString());
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
      time.set(query.getColumn("completion_time").getString());
      proc.setCompletionTime(time);
      if (have_meta_info)
      {
        Key id = query.getColumn("id").getInt64();
        handleQueryMetaInfo_(subquery_info, proc, id);
      }
      data_processing.push_back(proc);
    }
  }


  BaseFeature OMSFileLoad::makeBaseFeature_(int id, SQLite::Statement& query_feat,
                                            SQLite::Statement& query_meta,
                                            SQLite::Statement& query_match)
  {
    BaseFeature feature;
    feature.setRT(query_feat.getColumn("rt").getDouble());
    feature.setMZ(query_feat.getColumn("mz").getDouble());
    feature.setIntensity(query_feat.getColumn("intensity").getDouble());
    feature.setCharge(query_feat.getColumn("charge").getInt());
    feature.setWidth(query_feat.getColumn("width").getDouble());
    string quality_column = (version_number_ < 5) ? "overall_quality" : "quality";
    feature.setQuality(query_feat.getColumn(quality_column.c_str()).getDouble());
    feature.setUniqueId(query_feat.getColumn("unique_id").getInt64());
    if (id == -1) return feature; // stop here for feature handles (in consensus maps)

    auto primary_id = query_feat.getColumn("primary_molecule_id"); // optional
    if (!primary_id.isNull())
    {
      feature.setPrimaryID(identified_molecule_vars_[primary_id.getInt64()]);
    }
    // meta data:
    if (!isEmpty_(query_meta))
    {
      handleQueryMetaInfo_(query_meta, feature, id);
    }
    // ID matches:
    if (!isEmpty_(query_match))
    {
      query_match.bind(":id", id);
      while (query_match.executeStep())
      {
        Key match_id = query_match.getColumn("observation_match_id").getInt64();
        feature.addIDMatch(observation_match_refs_[match_id]);
      }
      query_match.reset(); // get ready for new executeStep()
    }
    return feature;
  }


  void OMSFileLoad::prepareQueriesBaseFeature_(SQLite::Statement& query_meta,
                                               SQLite::Statement& query_match)
  {
    // the "main" query is different for Feature/ConsensusFeature, so don't include it here
    string main_table = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature";
    prepareQueryMetaInfo_(query_meta, main_table);
    if (db_->tableExists("FEAT_ObservationMatch"))
    {
      query_match = SQLite::Statement(*db_, "SELECT * FROM FEAT_ObservationMatch WHERE feature_id = :id");
    }
  }


  Feature OMSFileLoad::loadFeatureAndSubordinates_(
    SQLite::Statement& query_feat, SQLite::Statement& query_meta,
    SQLite::Statement& query_match, SQLite::Statement& query_hull)
  {
    int id = query_feat.getColumn("id").getInt();
    Feature feature(makeBaseFeature_(id, query_feat, query_meta, query_match));
    // Feature-specific attributes:
    feature.setQuality(0, query_feat.getColumn("rt_quality").getDouble());
    feature.setQuality(1, query_feat.getColumn("mz_quality").getDouble());
    // convex hulls:
    if (!isEmpty_(query_hull))
    {
      query_hull.bind(":id", id);
      while (query_hull.executeStep())
      {
        Size hull_index = query_hull.getColumn("hull_index").getUInt();
        // first row should have max. hull index (sorted descending):
        if (feature.getConvexHulls().size() <= hull_index)
        {
          feature.getConvexHulls().resize(hull_index + 1);
        }
        ConvexHull2D::PointType point(query_hull.getColumn("point_x").getDouble(),
                                      query_hull.getColumn("point_y").getDouble());
        // @TODO: this may be inefficient (see implementation of "addPoint"):
        feature.getConvexHulls()[hull_index].addPoint(point);
      }
      query_hull.reset(); // get ready for new executeStep()
    }
    // subordinates:
    string from = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature JOIN FEAT_Feature ON id = feature_id";
    SQLite::Statement query_sub(*db_, "SELECT * FROM " + from + " WHERE subordinate_of = " + String(id) + " ORDER BY id ASC");
    while (query_sub.executeStep())
    {
      Feature sub = loadFeatureAndSubordinates_(query_sub, query_meta,
                                                query_match, query_hull);
      feature.getSubordinates().push_back(sub);
    }
    return feature;
  }


  void OMSFileLoad::loadFeatures_(FeatureMap& features)
  {
    if (!db_->tableExists("FEAT_Feature")) return;

    // start with top-level features only:
    string from = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature JOIN FEAT_Feature ON id = feature_id";
    SQLite::Statement query_feat(*db_, "SELECT * FROM " + from + " WHERE subordinate_of IS NULL ORDER BY id ASC");
    // prepare sub-queries (optional - corresponding tables may not be present):
    SQLite::Statement query_meta(*db_, "");
    SQLite::Statement query_match(*db_, "");
    prepareQueriesBaseFeature_(query_meta, query_match);
    SQLite::Statement query_hull(*db_, "");
    if (db_->tableExists("FEAT_ConvexHull"))
    {
      query_hull = SQLite::Statement(*db_, "SELECT * FROM FEAT_ConvexHull WHERE feature_id = :id " \
                                     "ORDER BY hull_index DESC, point_index ASC");
    }

    while (query_feat.executeStep())
    {
      Feature feature = loadFeatureAndSubordinates_(query_feat, query_meta,
                                                    query_match, query_hull);
      features.push_back(feature);
    }
  }


  void OMSFileLoad::load(FeatureMap& features)
  {
    load(features.getIdentificationData()); // load IDs, if any
    startProgress(0, 3, "Reading feature data from file");
    loadMapMetaData_(features);
    nextProgress();
    loadDataProcessing_(features.getDataProcessing());
    nextProgress();
    loadFeatures_(features);
    endProgress();
  }


  void OMSFileLoad::loadConsensusFeatures_(ConsensusMap& consensus)
  {
    if (!db_->tableExists("FEAT_FeatureHandle")) return;

    // start with top-level features only:
    SQLite::Statement query_feat(*db_, "SELECT * FROM FEAT_BaseFeature LEFT JOIN FEAT_FeatureHandle ON id = feature_id ORDER BY id ASC");
    // prepare sub-queries (optional - corresponding tables may not be present):
    SQLite::Statement query_meta(*db_, "");
    SQLite::Statement query_match(*db_, "");
    prepareQueriesBaseFeature_(query_meta, query_match);
    SQLite::Statement query_ratio(*db_, "");
    if (db_->tableExists("FEAT_ConsensusRatio"))
    {
      query_ratio = SQLite::Statement(*db_, "SELECT * FROM FEAT_ConsensusRatio WHERE feature_id = :id " \
                                      "ORDER BY ratio_index DESC");
    }

    while (query_feat.executeStep())
    {
      if (query_feat.getColumn("subordinate_of").isNull()) // ConsensusFeature
      {
        int id = query_feat.getColumn("id").getInt();
        ConsensusFeature feature(makeBaseFeature_(id, query_feat, query_meta, query_match));
        consensus.push_back(feature);
        if (!isEmpty_(query_ratio))
        {
          query_ratio.bind(":id", id);
          while (query_ratio.executeStep())
          {
            Size ratio_index = query_ratio.getColumn("ratio_index").getUInt();
            // first row should have max. hull index (sorted descending):
            if (feature.getRatios().size() <= ratio_index)
            {
              feature.getRatios().resize(ratio_index + 1);
            }
            ConsensusFeature::Ratio& ratio = feature.getRatios()[ratio_index];
            ratio.ratio_value_ = query_ratio.getColumn("ratio_value").getDouble();
            ratio.denominator_ref_ = query_ratio.getColumn("denominator_ref").getString();
            ratio.numerator_ref_ = query_ratio.getColumn("numerator_ref").getString();
            ratio.description_ = ListUtils::create<String>(query_ratio.getColumn("description").getString());
          }
          query_ratio.reset(); // get ready for new executeStep()
        }
      }
      else // FeatureHandle
      {
        BaseFeature feature(makeBaseFeature_(-1, query_feat, query_meta, query_match));
        UInt64 map_index = query_feat.getColumn("map_index").getInt64();
        FeatureHandle handle(map_index, feature);
        consensus.back().insert(handle);
      }
    }
  }


  void OMSFileLoad::loadConsensusColumnHeaders_(ConsensusMap& consensus)
  {
    consensus.getColumnHeaders().clear();
    if (!db_->tableExists("FEAT_ConsensusColumnHeader")) return;

    SQLite::Statement query(*db_, "SELECT * FROM FEAT_ConsensusColumnHeader");
    SQLite::Statement query_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(query_info, "FEAT_ConsensusColumnHeader");
    while (query.executeStep())
    {
      UInt64 id = query.getColumn("id").getInt64();
      ConsensusMap::ColumnHeader header;
      header.filename = query.getColumn("filename").getString();
      header.label = query.getColumn("label").getString();
      header.size = query.getColumn("size").getInt64();
      header.unique_id = query.getColumn("unique_id").getInt64();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(query_info, header, id);
      }
      consensus.getColumnHeaders()[id] = header;
    }
  }


  void OMSFileLoad::load(ConsensusMap& consensus)
  {
    load(consensus.getIdentificationData()); // load IDs, if any
    startProgress(0, 4, "Reading feature data from file");
    loadMapMetaData_(consensus);
    nextProgress();
    loadConsensusColumnHeaders_(consensus);
    nextProgress();
    loadDataProcessing_(consensus.getDataProcessing());
    nextProgress();
    loadConsensusFeatures_(consensus);
    endProgress();
  }


  QJsonArray OMSFileLoad::exportTableToJSON_(const QString& table, const QString& order_by)
  {
    // code based on: https://stackoverflow.com/a/18067555
    String sql = "SELECT * FROM " + table;
    if (!order_by.isEmpty())
    {
      sql += " ORDER BY " + order_by;
    }

    SQLite::Statement query(*db_, sql);

    QJsonArray array;
    while (query.executeStep())
    {
      QJsonObject record;
      for (int i = 0; i < query.getColumnCount(); ++i)
      {
        // @TODO: this will repeat field names for every row -
        // avoid this with separate "header" and "rows" (array)?

        // sqlite stores each cell based on the actual value, not the declared column type;
        // thus, we could use query.getColumnDeclaredType(i), but it would incur conversion
        switch (query.getColumn(i).getType())
        {
          case SQLITE_INTEGER: record.insert(query.getColumnName(i), qint64(query.getColumn(i).getInt64())); break;
          case SQLITE_FLOAT: record.insert(query.getColumnName(i), query.getColumn(i).getDouble()); break;
          case SQLITE_BLOB: record.insert(query.getColumnName(i), query.getColumn(i).getText()); break;
          case SQLITE_NULL: record.insert(query.getColumnName(i), ""); break;
          case SQLITE3_TEXT: record.insert(query.getColumnName(i), query.getColumn(i).getText()); break;
          default:
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
      }
      array.push_back(record);
    }
    return array;
  }


  void OMSFileLoad::exportToJSON(ostream& output)
  {
    // @TODO: this constructs the whole JSON file in memory - write directly to stream instead?
    // (more code, but would use less memory)
    QJsonObject json_data;
    // get names of all tables (except SQLite-internal ones) in the database:
    SQLite::Statement query(*db_, "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name");
    while (query.executeStep())
    {
      String table = query.getColumn("name").getString();
      QString order_by = "id"; // row order for most tables
      // special cases regarding ordering, e.g. tables without "id" column:
      if (table.hasSuffix("_MetaInfo"))
      {
        order_by = "parent_id, name";
      }
      else if (table.hasSuffix("_AppliedProcessingStep"))
      {
        order_by = "parent_id, processing_step_order, score_type_id";
      }
      else if (auto pos = export_order_by_.find(table.toQString()); pos != export_order_by_.end())
      {
        order_by = pos->second;
      }
      json_data.insert(table.toQString(), exportTableToJSON_(table.toQString(), order_by));
    }

    QJsonDocument json_doc;
    json_doc.setObject(json_data);
    output << json_doc.toJson().toStdString();
  }
}
