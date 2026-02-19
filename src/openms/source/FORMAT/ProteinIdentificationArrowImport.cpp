// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProteinIdentificationArrowImport.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{

namespace // anonymous
{
  /// Read a single Parquet file into an Arrow table.
  std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to open file: " << filename << std::endl;
      return nullptr;
    }
    auto infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to create Parquet reader for: " << filename << std::endl;
      return nullptr;
    }
    auto reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to read table: " << read_status.ToString() << std::endl;
      return nullptr;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to combine chunks" << std::endl;
      return nullptr;
    }

    return *combined;
  }

  /// Fetch a named column from a table, combining chunks if needed.
  /// Returns nullptr if column not found (logs error if required).
  std::shared_ptr<arrow::Array> getColumn_(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    bool required = true)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      if (required)
      {
        OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Missing required column '" << name << "'" << std::endl;
      }
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Column '" << name << "' has no chunks" << std::endl;
      return nullptr;
    }
    if (column->num_chunks() == 1)
    {
      return column->chunk(0);
    }
    auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
    if (!combined.ok())
    {
      OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to combine chunks for column '" << name << "'" << std::endl;
      return nullptr;
    }
    return *combined;
  }

  /// Get string value at a row, returning empty string if null.
  String getStringValue_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array || array->IsNull(row)) return "";
    return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
  }

  /// Get double value at a row, returning default_val if null.
  double getDoubleValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val = 0.0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
  }

  /// Get int32 value at a row, returning default_val if null.
  int32_t getInt32Value_(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val = 0)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
  }

  /// Get boolean value at a row, returning default_val if null.
  bool getBoolValue_(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val = false)
  {
    if (!array || array->IsNull(row)) return default_val;
    return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
  }

  /// Check if value at row is null.
  bool isNull_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    return !array || array->IsNull(row);
  }

  /// Read a list<utf8> column at a given row into a vector of Strings.
  std::vector<String> readStringList_(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    std::vector<String> result;
    if (!array || array->IsNull(row)) return result;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto values = std::static_pointer_cast<arrow::StringArray>(list_arr->value_slice(row));
    result.reserve(values->length());
    for (int64_t i = 0; i < values->length(); ++i)
    {
      result.emplace_back(values->GetString(i));
    }
    return result;
  }

  /// Read metavalues from a list<struct{name,value,value_type}> column at a given row.
  /// Sets them on the target MetaInfoInterface, excluding specified keys.
  void readMetaValues_(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    MetaInfoInterface& target,
    const std::unordered_set<std::string>& excluded_keys = {})
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto value_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));
    auto type_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(2));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      std::string name = name_arr->GetString(i);
      if (excluded_keys.count(name)) continue;

      std::string value_str = value_arr->GetString(i);
      std::string type_str = type_arr->GetString(i);

      if (type_str == "int")
      {
        try { target.setMetaValue(name, static_cast<int>(std::stol(value_str))); }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else if (type_str == "float")
      {
        try { target.setMetaValue(name, std::stod(value_str)); }
        catch (...) { target.setMetaValue(name, value_str); }
      }
      else
      {
        target.setMetaValue(name, value_str);
      }
    }
  }

  /// Build a map from run_identifier to index in the protein_identifications vector.
  std::unordered_map<std::string, size_t> buildRunIdMap_(
    const std::vector<ProteinIdentification>& protein_identifications)
  {
    std::unordered_map<std::string, size_t> map;
    for (size_t i = 0; i < protein_identifications.size(); ++i)
    {
      map[protein_identifications[i].getIdentifier()] = i;
    }
    return map;
  }

} // anonymous namespace


bool ProteinIdentificationArrowImport::importSearchParamsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Get all columns
  auto col_run_id = getColumn_(table, "run_identifier");
  auto col_search_engine = getColumn_(table, "search_engine");
  auto col_se_version = getColumn_(table, "search_engine_version", false);
  auto col_inf_engine = getColumn_(table, "inference_engine", false);
  auto col_inf_version = getColumn_(table, "inference_engine_version", false);
  auto col_date = getColumn_(table, "date", false);
  auto col_score_type = getColumn_(table, "score_type");
  auto col_higher_better = getColumn_(table, "higher_score_better");
  auto col_sig_threshold = getColumn_(table, "significance_threshold", false);
  auto col_db = getColumn_(table, "db", false);
  auto col_db_version = getColumn_(table, "db_version", false);
  auto col_taxonomy = getColumn_(table, "taxonomy", false);
  auto col_charges = getColumn_(table, "charges", false);
  auto col_mass_type = getColumn_(table, "mass_type");
  auto col_precursor_tol = getColumn_(table, "precursor_mass_tolerance");
  auto col_precursor_ppm = getColumn_(table, "precursor_mass_tolerance_ppm");
  auto col_fragment_tol = getColumn_(table, "fragment_mass_tolerance");
  auto col_fragment_ppm = getColumn_(table, "fragment_mass_tolerance_ppm");
  auto col_enzyme = getColumn_(table, "digestion_enzyme", false);
  auto col_enzyme_spec = getColumn_(table, "enzyme_term_specificity", false);
  auto col_missed_cleavages = getColumn_(table, "missed_cleavages");
  auto col_fixed_mods = getColumn_(table, "fixed_modifications", false);
  auto col_var_mods = getColumn_(table, "variable_modifications", false);
  auto col_ms_run_paths = getColumn_(table, "primary_ms_run_paths", false);
  auto col_metavalues = getColumn_(table, "metavalues", false);
  auto col_sp_metavalues = getColumn_(table, "sp_metavalues", false);

  if (!col_run_id || !col_search_engine || !col_score_type || !col_higher_better ||
      !col_mass_type || !col_precursor_tol || !col_precursor_ppm ||
      !col_fragment_tol || !col_fragment_ppm || !col_missed_cleavages)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Missing required columns in search_params table" << std::endl;
    return false;
  }

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    ProteinIdentification prot_id;

    // Run-level metadata
    prot_id.setIdentifier(getStringValue_(col_run_id, row));
    prot_id.setSearchEngine(getStringValue_(col_search_engine, row));
    prot_id.setSearchEngineVersion(getStringValue_(col_se_version, row));
    prot_id.setScoreType(getStringValue_(col_score_type, row));
    prot_id.setHigherScoreBetter(getBoolValue_(col_higher_better, row));
    prot_id.setSignificanceThreshold(getDoubleValue_(col_sig_threshold, row));

    // Date
    String date_str = getStringValue_(col_date, row);
    if (!date_str.empty())
    {
      DateTime dt;
      // Format: yyyy-MM-ddThh:mm:ss
      dt.set(date_str.substitute('T', ' '));
      prot_id.setDateTime(dt);
    }

    // Primary MS run paths
    auto ms_run_paths = readStringList_(col_ms_run_paths, row);
    if (!ms_run_paths.empty())
    {
      StringList sl(ms_run_paths.begin(), ms_run_paths.end());
      prot_id.setPrimaryMSRunPath(sl);
    }

    // SearchParameters
    ProteinIdentification::SearchParameters sp;
    sp.db = getStringValue_(col_db, row);
    sp.db_version = getStringValue_(col_db_version, row);
    sp.taxonomy = getStringValue_(col_taxonomy, row);
    sp.charges = getStringValue_(col_charges, row);

    String mass_type_str = getStringValue_(col_mass_type, row);
    sp.mass_type = (mass_type_str == "AVERAGE") ?
      ProteinIdentification::PeakMassType::AVERAGE :
      ProteinIdentification::PeakMassType::MONOISOTOPIC;

    sp.precursor_mass_tolerance = getDoubleValue_(col_precursor_tol, row);
    sp.precursor_mass_tolerance_ppm = getBoolValue_(col_precursor_ppm, row);
    sp.fragment_mass_tolerance = getDoubleValue_(col_fragment_tol, row);
    sp.fragment_mass_tolerance_ppm = getBoolValue_(col_fragment_ppm, row);
    sp.missed_cleavages = static_cast<UInt>(getInt32Value_(col_missed_cleavages, row));

    // Enzyme
    String enzyme_name = getStringValue_(col_enzyme, row);
    if (!enzyme_name.empty())
    {
      sp.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme_name);
    }

    // Enzyme term specificity
    String spec_str = getStringValue_(col_enzyme_spec, row);
    if (spec_str == "FULL") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_FULL;
    else if (spec_str == "SEMI") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_SEMI;
    else if (spec_str == "NONE") sp.enzyme_term_specificity = EnzymaticDigestion::SPEC_NONE;
    else if (!spec_str.empty())
    {
      OPENMS_LOG_WARN << "ProteinIdentificationArrowImport: Unknown enzyme_term_specificity '"
                      << spec_str << "'" << std::endl;
    }

    // Modifications
    auto fixed_mods = readStringList_(col_fixed_mods, row);
    sp.fixed_modifications.assign(fixed_mods.begin(), fixed_mods.end());
    auto var_mods = readStringList_(col_var_mods, row);
    sp.variable_modifications.assign(var_mods.begin(), var_mods.end());

    // SearchParameters metavalues (restore before setSearchParameters)
    readMetaValues_(col_sp_metavalues, row, sp);

    prot_id.setSearchParameters(sp);

    // Inference engine (must be set after setSearchParameters)
    String inf_engine = getStringValue_(col_inf_engine, row);
    if (!inf_engine.empty())
    {
      prot_id.setInferenceEngine(inf_engine);
      String inf_version = getStringValue_(col_inf_version, row);
      if (!inf_version.empty())
      {
        prot_id.setInferenceEngineVersion(inf_version);
      }
    }

    // ProteinIdentification metavalues (exclude keys handled via dedicated columns)
    static const std::unordered_set<std::string> excluded_prot_id_mvs = {
      "InferenceEngine", "InferenceEngineVersion",
      "spectra_data", "spectra_data_raw"
    };
    readMetaValues_(col_metavalues, row, prot_id, excluded_prot_id_mvs);

    protein_identifications.push_back(std::move(prot_id));
  }

  return true;
}


bool ProteinIdentificationArrowImport::importProteinsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Get all columns
  auto col_accession = getColumn_(table, "accession");
  auto col_score = getColumn_(table, "score");
  auto col_score_type = getColumn_(table, "score_type");
  auto col_higher_better = getColumn_(table, "higher_score_better");
  auto col_rank = getColumn_(table, "rank", false);
  auto col_coverage = getColumn_(table, "coverage", false);
  auto col_sequence = getColumn_(table, "sequence", false);
  auto col_description = getColumn_(table, "description", false);
  auto col_is_decoy = getColumn_(table, "is_decoy", false);
  auto col_run_id = getColumn_(table, "run_identifier");
  auto col_ref_file = getColumn_(table, "reference_file_name", false);
  auto col_search_engine = getColumn_(table, "search_engine");
  auto col_se_version = getColumn_(table, "search_engine_version", false);
  auto col_inf_engine = getColumn_(table, "inference_engine", false);
  auto col_inf_version = getColumn_(table, "inference_engine_version", false);
  auto col_sig_threshold = getColumn_(table, "significance_threshold", false);
  auto col_date = getColumn_(table, "date", false);
  auto col_modifications = getColumn_(table, "modifications", false);
  auto col_metavalues = getColumn_(table, "metavalues", false);

  if (!col_accession || !col_score || !col_score_type || !col_higher_better ||
      !col_run_id || !col_search_engine)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Missing required columns in proteins table" << std::endl;
    return false;
  }

  // Build map from run_identifier to index
  auto run_id_map = buildRunIdMap_(protein_identifications);

  // Metavalue keys excluded from hit metavalues (they have dedicated columns)
  static const std::unordered_set<std::string> excluded_hit_mvs = {
    "Description", "target_decoy"
  };

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    String run_id = getStringValue_(col_run_id, row);

    // Find or create matching ProteinIdentification
    auto it = run_id_map.find(run_id);
    size_t prot_id_idx;
    if (it != run_id_map.end())
    {
      prot_id_idx = it->second;
    }
    else
    {
      // Create new ProteinIdentification from run-level data in this row
      ProteinIdentification new_prot_id;
      new_prot_id.setIdentifier(run_id);
      new_prot_id.setSearchEngine(getStringValue_(col_search_engine, row));
      new_prot_id.setSearchEngineVersion(getStringValue_(col_se_version, row));
      new_prot_id.setScoreType(getStringValue_(col_score_type, row));
      new_prot_id.setHigherScoreBetter(getBoolValue_(col_higher_better, row));
      new_prot_id.setSignificanceThreshold(getDoubleValue_(col_sig_threshold, row));

      String date_str = getStringValue_(col_date, row);
      if (!date_str.empty())
      {
        DateTime dt;
        dt.set(date_str.substitute('T', ' '));
        new_prot_id.setDateTime(dt);
      }

      String ref_file = getStringValue_(col_ref_file, row);
      if (!ref_file.empty())
      {
        new_prot_id.setPrimaryMSRunPath(StringList{ref_file});
      }

      String inf_engine = getStringValue_(col_inf_engine, row);
      if (!inf_engine.empty())
      {
        new_prot_id.setInferenceEngine(inf_engine);
        String inf_version = getStringValue_(col_inf_version, row);
        if (!inf_version.empty())
        {
          new_prot_id.setInferenceEngineVersion(inf_version);
        }
      }

      prot_id_idx = protein_identifications.size();
      protein_identifications.push_back(std::move(new_prot_id));
      run_id_map[run_id] = prot_id_idx;
    }

    // Build ProteinHit
    ProteinHit hit;
    hit.setAccession(getStringValue_(col_accession, row));
    hit.setScore(getDoubleValue_(col_score, row));

    if (!isNull_(col_rank, row))
    {
      hit.setRank(static_cast<UInt>(getInt32Value_(col_rank, row)));
    }

    if (!isNull_(col_coverage, row))
    {
      hit.setCoverage(getDoubleValue_(col_coverage, row));
    }

    if (!isNull_(col_sequence, row))
    {
      hit.setSequence(getStringValue_(col_sequence, row));
    }

    if (!isNull_(col_description, row))
    {
      hit.setDescription(getStringValue_(col_description, row));
    }

    // is_decoy -> target_decoy metavalue
    if (!isNull_(col_is_decoy, row))
    {
      int32_t is_decoy = getInt32Value_(col_is_decoy, row);
      hit.setMetaValue("target_decoy", is_decoy == 1 ? "decoy" : "target");
    }

    // Modifications
    if (col_modifications && !col_modifications->IsNull(row))
    {
      auto list_arr = std::static_pointer_cast<arrow::ListArray>(col_modifications);
      auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
      if (struct_arr && struct_arr->length() > 0)
      {
        auto pos_arr = std::static_pointer_cast<arrow::Int32Array>(struct_arr->field(0));
        auto mod_name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));

        std::set<std::pair<Size, ResidueModification>> mods;
        for (int64_t i = 0; i < struct_arr->length(); ++i)
        {
          Size position = static_cast<Size>(pos_arr->Value(i));
          String mod_id = mod_name_arr->GetString(i);
          try
          {
            const ResidueModification* mod = ModificationsDB::getInstance()->getModification(mod_id);
            mods.insert(std::make_pair(position, *mod));
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "ProteinIdentificationArrowImport: Could not find modification '" << mod_id << "'" << std::endl;
          }
        }
        hit.setModifications(mods);
      }
    }

    // MetaValues
    readMetaValues_(col_metavalues, row, hit, excluded_hit_mvs);

    protein_identifications[prot_id_idx].insertHit(std::move(hit));
  }

  return true;
}


bool ProteinIdentificationArrowImport::importProteinGroupsFromArrow(
  const std::shared_ptr<arrow::Table>& table,
  std::vector<ProteinIdentification>& protein_identifications)
{
  if (!table || table->num_rows() == 0) return true;

  // Get all columns
  auto col_group_type = getColumn_(table, "group_type");
  auto col_probability = getColumn_(table, "probability");
  auto col_accessions = getColumn_(table, "accessions");
  auto col_run_id = getColumn_(table, "run_identifier");
  auto col_group_index = getColumn_(table, "group_index", false); // informational, not needed for reconstruction
  auto col_float_data = getColumn_(table, "float_data", false);
  auto col_string_data = getColumn_(table, "string_data", false);
  auto col_integer_data = getColumn_(table, "integer_data", false);

  if (!col_group_type || !col_probability || !col_accessions || !col_run_id)
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Missing required columns in protein_groups table" << std::endl;
    return false;
  }

  auto run_id_map = buildRunIdMap_(protein_identifications);

  // Helper lambda to read data arrays from list<struct{name,values: list<T>}>
  auto readFloatDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::FloatDataArray fda;
      fda.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::DoubleArray>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        fda.push_back(static_cast<float>(vals->Value(j)));
      }
      group.getFloatDataArrays().push_back(std::move(fda));
    }
  };

  auto readStringDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                 ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::StringDataArray sda;
      sda.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::StringArray>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        sda.push_back(vals->GetString(j));
      }
      group.getStringDataArrays().push_back(std::move(sda));
    }
  };

  auto readIntegerDataArrays = [](const std::shared_ptr<arrow::Array>& array, int64_t row,
                                  ProteinIdentification::ProteinGroup& group)
  {
    if (!array || array->IsNull(row)) return;
    auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
    auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
    if (!struct_arr || struct_arr->length() == 0) return;

    auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
    auto values_list_arr = std::static_pointer_cast<arrow::ListArray>(struct_arr->field(1));

    for (int64_t i = 0; i < struct_arr->length(); ++i)
    {
      ProteinIdentification::ProteinGroup::IntegerDataArray ida;
      ida.setName(name_arr->GetString(i));
      auto vals = std::static_pointer_cast<arrow::Int64Array>(values_list_arr->value_slice(i));
      for (int64_t j = 0; j < vals->length(); ++j)
      {
        ida.push_back(static_cast<Int>(vals->Value(j)));
      }
      group.getIntegerDataArrays().push_back(std::move(ida));
    }
  };

  for (int64_t row = 0; row < table->num_rows(); ++row)
  {
    String run_id = getStringValue_(col_run_id, row);

    // Find or create matching ProteinIdentification
    auto it = run_id_map.find(run_id);
    size_t prot_id_idx;
    if (it != run_id_map.end())
    {
      prot_id_idx = it->second;
    }
    else
    {
      ProteinIdentification new_prot_id;
      new_prot_id.setIdentifier(run_id);
      prot_id_idx = protein_identifications.size();
      protein_identifications.push_back(std::move(new_prot_id));
      run_id_map[run_id] = prot_id_idx;
    }

    // Build ProteinGroup
    ProteinIdentification::ProteinGroup group;
    group.probability = getDoubleValue_(col_probability, row);
    auto accessions = readStringList_(col_accessions, row);
    group.accessions.assign(accessions.begin(), accessions.end());

    // Data arrays
    readFloatDataArrays(col_float_data, row, group);
    readStringDataArrays(col_string_data, row, group);
    readIntegerDataArrays(col_integer_data, row, group);

    // Insert into appropriate vector
    String group_type = getStringValue_(col_group_type, row);
    if (group_type == "indistinguishable")
    {
      protein_identifications[prot_id_idx].insertIndistinguishableProteins(std::move(group));
    }
    else
    {
      protein_identifications[prot_id_idx].insertProteinGroup(std::move(group));
    }
  }

  return true;
}


bool ProteinIdentificationArrowImport::importSearchParamsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importSearchParamsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowImport::importProteinsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importProteinsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowImport::importProteinGroupsFromParquet(
  const String& filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  auto table = readParquetTable_(filename);
  if (!table) return false;
  return importProteinGroupsFromArrow(table, protein_identifications);
}


bool ProteinIdentificationArrowImport::importFromParquet(
  const String& proteins_filename,
  const String& protein_groups_filename,
  const String& search_params_filename,
  std::vector<ProteinIdentification>& protein_identifications)
{
  protein_identifications.clear();

  // 1. Import search params first (creates ProteinIdentification shells)
  if (!importSearchParamsFromParquet(search_params_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to import search params" << std::endl;
    return false;
  }

  // 2. Import proteins second (adds ProteinHits to matching shells)
  if (!importProteinsFromParquet(proteins_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to import proteins" << std::endl;
    return false;
  }

  // 3. Import protein groups third (adds groups to matching shells)
  if (!importProteinGroupsFromParquet(protein_groups_filename, protein_identifications))
  {
    OPENMS_LOG_ERROR << "ProteinIdentificationArrowImport: Failed to import protein groups" << std::endl;
    return false;
  }

  return true;
}


} // namespace OpenMS

#endif // WITH_PARQUET
