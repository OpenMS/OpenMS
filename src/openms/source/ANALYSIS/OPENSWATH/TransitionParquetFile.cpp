// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>

#include <filesystem>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <fstream>
#include <arrow/api.h>

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <map>
#include <sstream>
#include <cctype>

namespace
{
  using OpenMS::Size;

  std::string joinProteinAccessions_(const std::vector<std::string>& accessions)
  {
    OpenMS::String joined;
    for (OpenMS::Size i = 0; i < accessions.size(); ++i)
    {
      if (i > 0) joined += ";";
      joined += accessions[i];
    }
    return std::string(joined);
  }

  struct OpenMSLibraryStats
  {
    Size proteins_total = 0;
    Size proteins_decoy = 0;
    Size peptides_total = 0;
    Size peptides_decoy = 0;
    Size precursors_total = 0;
    Size precursors_decoy = 0;
    Size compounds_total = 0;
    Size compounds_decoy = 0;
    Size transitions_total = 0;
    Size transitions_decoy = 0;
    Size fragment_b_target = 0;
    Size fragment_b_decoy = 0;
    Size fragment_y_target = 0;
    Size fragment_y_decoy = 0;
    Size fragment_other_target = 0;
    Size fragment_other_decoy = 0;
    std::map<int, Size> precursor_charge_counts_target;
    std::map<int, Size> precursor_charge_counts_decoy;
    std::map<int, Size> transition_charge_counts_target;
    std::map<int, Size> transition_charge_counts_decoy;
  };

  std::string jsonEscape_(const OpenMS::String& input)
  {
    return OpenMS::ParquetFile::jsonEscape(input);
  }

  std::string jsonMapByClass_(const std::map<int, Size>& target, const std::map<int, Size>& decoy)
  {
    std::ostringstream oss;
    oss << "{";
    // target
    oss << "\"target\":{";
    bool first = true;
    for (const auto& p : target)
    {
      if (!first) oss << ",";
      first = false;
      oss << "\"" << p.first << "\":" << p.second;
    }
    oss << "},";
    // decoy
    oss << "\"decoy\":{";
    first = true;
    for (const auto& p : decoy)
    {
      if (!first) oss << ",";
      first = false;
      oss << "\"" << p.first << "\":" << p.second;
    }
    oss << "}";
    oss << "}";
    return oss.str();
  }

  void writeLibraryMetadata_(const OpenMS::String& library_dir, const OpenMS::String& library_name, const OpenMSLibraryStats& stats)
  {
    (void)library_name;
    const Size proteins_target = stats.proteins_total - stats.proteins_decoy;
    const Size peptides_target = stats.peptides_total - stats.peptides_decoy;
    const Size precursors_target = stats.precursors_total - stats.precursors_decoy;
    const Size compounds_target = stats.compounds_total - stats.compounds_decoy;
    const Size transitions_target = stats.transitions_total - stats.transitions_decoy;

    const OpenMS::String metadata_path = library_dir + "/metadata.json";
    std::ofstream out(metadata_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
      throw OpenMS::Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, metadata_path);
    }

    out << "{\n"
        << "  \"openms\": {\n"
        << "    \"schema_version\": 1,\n"
        << "    \"generator\": \"TransitionParquetFile\",\n"
      << "    \"openms_version\": \"" << jsonEscape_(OpenMS::VersionInfo::getVersion()) << "\",\n"
        << "    \"build_time\": \"" << jsonEscape_(OpenMS::VersionInfo::getTime()) << "\",\n"
        << "    \"tool\": {\"name\": \"OpenSwathWorkflow\", \"version\": \"" << jsonEscape_(OpenMS::VersionInfo::getVersion()) << "\"},\n"
        << "    \"counts\": {\n"
        << "      \"proteins\": {\"total\": " << stats.proteins_total << ", \"target\": " << proteins_target << ", \"decoy\": " << stats.proteins_decoy << "},\n"
        << "      \"peptides\": {\"total\": " << stats.peptides_total << ", \"target\": " << peptides_target << ", \"decoy\": " << stats.peptides_decoy << "},\n"
        << "      \"precursors\": {\"total\": " << stats.precursors_total << ", \"target\": " << precursors_target << ", \"decoy\": " << stats.precursors_decoy << "},\n"
        << "      \"compounds\": {\"total\": " << stats.compounds_total << ", \"target\": " << compounds_target << ", \"decoy\": " << stats.compounds_decoy << "},\n"
        << "      \"transitions\": {\"total\": " << stats.transitions_total << ", \"target\": " << transitions_target << ", \"decoy\": " << stats.transitions_decoy << "}\n"
        << "    },\n"
        << "    \"fragment_type_counts\": {\n"
        << "      \"target\": {\"b\": " << stats.fragment_b_target << ", \"y\": " << stats.fragment_y_target << ", \"other\": " << stats.fragment_other_target << "},\n"
        << "      \"decoy\": {\"b\": " << stats.fragment_b_decoy << ", \"y\": " << stats.fragment_y_decoy << ", \"other\": " << stats.fragment_other_decoy << "}\n"
        << "    },\n"
        << "    \"charge_counts\": {\n"
        << "      \"precursor\": " << jsonMapByClass_(stats.precursor_charge_counts_target, stats.precursor_charge_counts_decoy) << ",\n"
        << "      \"transition\": " << jsonMapByClass_(stats.transition_charge_counts_target, stats.transition_charge_counts_decoy) << "\n"
        << "    }\n"
        << "  }\n"
        << "}\n";
  }

  struct PrecursorInfo
  {
    double precursor_mz = 0.0;
    double library_rt = 0.0;
    double drift_time = -1.0;
    int charge = 0;
    bool decoy = false;
    std::string traml_id;
    std::string modified_sequence;
    std::string unmodified_sequence;
    std::vector<std::string> protein_accessions;
  };

  void addModifications_(const std::string& sequence, OpenSwath::LightCompound& compound)
  {
    if (sequence.empty()) return;
    try
    {
      OpenMS::AASequence aa_sequence = OpenMS::AASequence::fromString(sequence);
      if (aa_sequence.hasNTerminalModification())
      {
        OpenSwath::LightModification mod;
        mod.location = -1;
        mod.unimod_id = aa_sequence.getNTerminalModification()->getUniModRecordId();
        compound.modifications.push_back(mod);
      }
      if (aa_sequence.hasCTerminalModification())
      {
        OpenSwath::LightModification mod;
        mod.location = static_cast<int>(aa_sequence.size());
        mod.unimod_id = aa_sequence.getCTerminalModification()->getUniModRecordId();
        compound.modifications.push_back(mod);
      }
      for (OpenMS::Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          OpenSwath::LightModification mod;
          mod.location = static_cast<int>(i);
          mod.unimod_id = aa_sequence.getResidue(i).getModification()->getUniModRecordId();
          compound.modifications.push_back(mod);
        }
      }
    }
    catch (OpenMS::Exception::InvalidValue&)
    {
      OPENMS_LOG_DEBUG << "Could not parse modifications from sequence: " << sequence << '\n';
    }
  }
} // namespace

namespace OpenMS
{
  void TransitionParquetFile::convertParquetToTargetedExperiment(
    const String& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
    // Reset the output container to avoid appending to a caller-owned
    // object that may contain stale data from previous calls. The caller
    // expects this function to populate `targeted_exp` from the parquet
    // files, not to append to it.
    targeted_exp = OpenSwath::LightTargetedExperiment{};
    std::unique_ptr<File::TempDir> temp_dir;

    // Try to open parquet entries directly from the archive using a RandomAccessFile.
    // If that fails (e.g., compressed entry or libzip not available), fall back to
    // extracting the entry to a temporary file and reading from disk.
    auto open_table_from_entry = [&](const String& entry) -> std::shared_ptr<arrow::Table>
    {
      auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
      if (ra_res.ok())
      {
        const auto& raf = ra_res.ValueOrDie();
        return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
      }
      // Fallback to extract
      const String path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
      if (!File::exists(path))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Missing required parquet file '" + entry + "' in '" + oswpq_dir + "'");
      }
      return ParquetFile::readTable(path);
    };

    auto precursors_table = open_table_from_entry("library/precursors.parquet");
    auto transitions_table = open_table_from_entry("library/transitions.parquet");

    // Validate loaded tables against registry schemas (subset mode — file may have extra columns)
    auto prec_validation = ArrowSchemaValidation::validate(precursors_table, OSWPrecursorSchema::schema(), ArrowSchemaValidation::Mode::Subset);
    if (!prec_validation.valid)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Precursors table schema validation failed: " + prec_validation.toString(), "");
    }
    auto trans_validation = ArrowSchemaValidation::validate(transitions_table, OSWTransitionSchema::schema(), ArrowSchemaValidation::Mode::Subset);
    if (!trans_validation.valid)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Transitions table schema validation failed: " + trans_validation.toString(), "");
    }

    auto precursor_id_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::PRECURSOR_ID);
    auto precursor_mz_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::PRECURSOR_MZ);
    auto charge_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::CHARGE);
    auto library_rt_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::LIBRARY_RT);
    auto drift_time_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::LIBRARY_DRIFT_TIME);
    auto traml_id_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::TRAML_ID);
    auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::DECOY);
    auto modified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::MODIFIED_SEQUENCE);
    auto unmodified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::UNMODIFIED_SEQUENCE);
    auto protein_accessions_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::PROTEIN_ACCESSIONS);

    std::unordered_map<int64_t, PrecursorInfo> precursor_map;
    precursor_map.reserve(precursors_table->num_rows());

    for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
    {
      const int64_t precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      PrecursorInfo info;
      info.precursor_mz = ParquetFile::getDouble(precursor_mz_col, row, 0.0, false);
      info.charge = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
      info.library_rt = ParquetFile::getDouble(library_rt_col, row, 0.0, true);
      info.drift_time = ParquetFile::getDouble(drift_time_col, row, -1.0, true);
      info.decoy = ParquetFile::getBool(decoy_col, row, false, true);
      info.traml_id = ParquetFile::getString(traml_id_col, row);
      info.modified_sequence = ParquetFile::getString(modified_sequence_col, row);
      info.unmodified_sequence = ParquetFile::getString(unmodified_sequence_col, row);
      info.protein_accessions = ParquetFile::getStringList(protein_accessions_col, row);

      precursor_map.emplace(precursor_id, std::move(info));
    }

    std::unordered_map<std::string, int> compound_map;
    std::unordered_map<std::string, int> protein_map;
    compound_map.reserve(precursor_map.size());

    for (const auto& entry : precursor_map)
    {
      const int64_t precursor_id = entry.first;
      const PrecursorInfo& info = entry.second;
      const String precursor_id_str(precursor_id);

    OpenSwath::LightCompound compound;
    // Preserve source traml_id when available to maintain round-trip identity
    // fidelity. If traml_id is empty, fall back to the numeric precursor id.
    const String compound_id = info.traml_id.empty() ? precursor_id_str : String(info.traml_id);
    compound.id = compound_id;
      compound.drift_time = info.drift_time;
      compound.rt = info.library_rt;
      compound.charge = info.charge;
      compound.sequence = info.modified_sequence.empty() ? info.unmodified_sequence : info.modified_sequence;
      compound.protein_refs = info.protein_accessions;
      if (!compound.sequence.empty())
      {
        addModifications_(compound.sequence, compound);
      }

      targeted_exp.compounds.push_back(std::move(compound));
      compound_map[precursor_id_str] = 0;

      for (const auto& accession : info.protein_accessions)
      {
        if (protein_map.find(accession) == protein_map.end())
        {
          OpenSwath::LightProtein protein;
          protein.id = accession;
          protein.sequence = "";
          targeted_exp.proteins.push_back(std::move(protein));
          protein_map[accession] = 0;
        }
      }
    }

    auto transition_id_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::TRANSITION_ID);
    auto transition_traml_id_col = ParquetFile::getOptionalColumn(transitions_table, OSWTransitionSchema::TRAML_ID);
    auto transition_precursor_id_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::PRECURSOR_ID);
    auto product_mz_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::PRODUCT_MZ);
    auto fragment_charge_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::CHARGE);
    auto fragment_type_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::TYPE);
    auto fragment_annotation_col = ParquetFile::getOptionalColumn(transitions_table, OSWTransitionSchema::ANNOTATION);
    auto fragment_ordinal_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::ORDINAL);
    auto detecting_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::DETECTING);
    auto identifying_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::IDENTIFYING);
    auto quantifying_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::QUANTIFYING);
    auto transition_intensity_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::LIBRARY_INTENSITY);
    auto transition_decoy_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::DECOY);

    std::unordered_set<std::string> used_transition_names;
    used_transition_names.reserve(transitions_table->num_rows());
    bool warned_duplicate_transition = false;

    for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
    {
      const int64_t precursor_id = ParquetFile::getInt64(transition_precursor_id_col, row, 0, false);
      auto precursor_it = precursor_map.find(precursor_id);
      if (precursor_it == precursor_map.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown precursor_id");
      }

      const int64_t transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
      const std::string traml_id = ParquetFile::getString(transition_traml_id_col, row);
      std::string transition_name = traml_id.empty() ? String(transition_id) : String(traml_id);
      if (!used_transition_names.insert(transition_name).second)
      {
        if (!warned_duplicate_transition)
        {
          OPENMS_LOG_WARN << "Duplicate transition nativeID detected in Parquet library. "
                          << "Falling back to transition_id for uniqueness." << std::endl;
          warned_duplicate_transition = true;
        }
        transition_name = String(transition_id);
        if (!used_transition_names.insert(transition_name).second)
        {
          transition_name += "_" + std::to_string(row);
          used_transition_names.insert(transition_name);
        }
      }
      const std::string fragment_type = ParquetFile::getString(fragment_type_col, row);
      const std::string annotation = ParquetFile::getString(fragment_annotation_col, row);

      OpenSwath::LightTransition transition;
      transition.transition_name = transition_name;
      // Use precursor traml_id as peptide_ref when present to preserve source IDs
      const String peptide_ref = precursor_it->second.traml_id.empty() ? String(precursor_id) : String(precursor_it->second.traml_id);
      transition.peptide_ref = peptide_ref;
      transition.library_intensity = ParquetFile::getDouble(transition_intensity_col, row, 0.0, true);
      transition.precursor_mz = precursor_it->second.precursor_mz;
      transition.product_mz = ParquetFile::getDouble(product_mz_col, row, 0.0, false);
      transition.precursor_im = precursor_it->second.drift_time;
      transition.fragment_charge = static_cast<int8_t>(ParquetFile::getInt64(fragment_charge_col, row, 0, true));
      transition.fragment_nr = static_cast<int16_t>(ParquetFile::getInt64(fragment_ordinal_col, row, -1, true));
      transition.setFragmentType(fragment_type.empty() ? annotation : fragment_type);
      transition.setDecoy(ParquetFile::getBool(transition_decoy_col, row, false, true));
      transition.setDetectingTransition(ParquetFile::getBool(detecting_col, row, true, true));
      transition.setIdentifyingTransition(ParquetFile::getBool(identifying_col, row, false, true));
      transition.setQuantifyingTransition(ParquetFile::getBool(quantifying_col, row, true, true));

      targeted_exp.transitions.push_back(std::move(transition));
    }
  }

  void TransitionParquetFile::convertLightTargetedExperimentToParquet(
    const String& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
    const bool output_is_dir = File::isDirectory(oswpq_path);
    std::unique_ptr<File::TempDir> temp_dir;
    String base_dir = oswpq_path;
    if (!output_is_dir)
    {
      temp_dir = std::make_unique<File::TempDir>();
      base_dir = temp_dir->getPath() + "/oswpq_output";
      File::makeDir(base_dir);
    }

    const String library_dir = base_dir + "/library";
    File::makeDir(library_dir);
    String library_name = File::basename(oswpq_path);
    if (library_name.empty())
    {
      library_name = "openms_library";
    }

    std::unordered_map<std::string, int64_t> compound_to_precursor;
    compound_to_precursor.reserve(targeted_exp.compounds.size());

    std::unordered_map<std::string, bool> compound_decoy;
    compound_decoy.reserve(targeted_exp.compounds.size());
    for (const auto& transition : targeted_exp.transitions)
    {
      if (transition.getDecoy())
      {
        compound_decoy[transition.peptide_ref] = true;
      }
    }

    OpenMSLibraryStats stats;

    int64_t next_precursor_id = 1;
    std::unordered_set<int64_t> used_precursor_ids;
    for (const auto& compound : targeted_exp.compounds)
    {
      if (compound_to_precursor.find(compound.id) != compound_to_precursor.end())
      {
        continue;
      }

      int64_t precursor_id = 0;
      bool parsed_numeric = false;
      try
      {
        precursor_id = OpenMS::String(compound.id).toInt64();
        parsed_numeric = true;
      }
      catch (OpenMS::Exception::ConversionError&)
      {
        // will assign auto id below
      }

      if (parsed_numeric)
      {
        if (precursor_id <= 0 || used_precursor_ids.find(precursor_id) != used_precursor_ids.end())
        {
          precursor_id = next_precursor_id++;
        }
        else
        {
          if (precursor_id >= next_precursor_id)
          {
            next_precursor_id = precursor_id + 1;
          }
        }
      }
      else
      {
        precursor_id = next_precursor_id++;
      }

      used_precursor_ids.insert(precursor_id);
      compound_to_precursor[compound.id] = precursor_id;
    }

    std::unordered_map<std::string, double> precursor_mz;
    for (const auto& transition : targeted_exp.transitions)
    {
      if (precursor_mz.find(transition.peptide_ref) == precursor_mz.end())
      {
        precursor_mz[transition.peptide_ref] = transition.precursor_mz;
      }
    }

    arrow::Int64Builder precursor_id_builder;
    arrow::DoubleBuilder precursor_mz_builder;
    arrow::Int32Builder precursor_charge_builder;
    arrow::DoubleBuilder library_rt_builder;
    arrow::DoubleBuilder drift_time_builder;
    arrow::BooleanBuilder decoy_builder;
    arrow::StringBuilder traml_id_builder;
    arrow::StringBuilder modified_sequence_builder;
    arrow::StringBuilder unmodified_sequence_builder;
    arrow::StringBuilder protein_accessions_builder;

    for (const auto& compound : targeted_exp.compounds)
    {
      const int64_t precursor_id = compound_to_precursor[compound.id];
      const bool is_decoy = compound_decoy[compound.id] ||
        OpenMS::String(compound.id).hasPrefix("DECOY_");
      stats.compounds_total++;
      stats.precursors_total++;
      if (compound.isPeptide())
      {
        stats.peptides_total++;
        if (is_decoy) stats.peptides_decoy++;
      }
      if (is_decoy)
      {
        stats.compounds_decoy++;
        stats.precursors_decoy++;
      }
      if (is_decoy)
      {
        stats.precursor_charge_counts_decoy[compound.charge]++;
      }
      else
      {
        stats.precursor_charge_counts_target[compound.charge]++;
      }
      ParquetFile::appendOrThrow(precursor_id_builder.Append(precursor_id), "precursor_id");
      // Guard lookup in precursor_mz: using operator[] would insert a default
      // 0.0 value if the key is missing. Instead, ensure the value exists and
      // report a clear error if a compound has no matched transition.
      auto mz_it = precursor_mz.find(compound.id);
      if (mz_it == precursor_mz.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "No precursor_mz found for compound '" + String(compound.id) + "'");
      }
      ParquetFile::appendOrThrow(precursor_mz_builder.Append(mz_it->second), "precursor_mz");
      ParquetFile::appendOrThrow(precursor_charge_builder.Append(compound.charge), "charge");
      ParquetFile::appendOrThrow(library_rt_builder.Append(compound.rt), "library_rt");
      ParquetFile::appendOrThrow(drift_time_builder.Append(compound.drift_time), "library_drift_time");
      ParquetFile::appendOrThrow(decoy_builder.Append(is_decoy), "decoy");
      ParquetFile::appendOrThrow(traml_id_builder.Append(compound.id), "traml_id");
      ParquetFile::appendOrThrow(modified_sequence_builder.Append(compound.sequence), "modified_sequence");

      std::string unmodified_sequence;
      if (!compound.sequence.empty())
      {
        try
        {
          unmodified_sequence = OpenMS::AASequence::fromString(compound.sequence).toUnmodifiedString();
        }
        catch (OpenMS::Exception::InvalidValue&)
        {
          unmodified_sequence = "";
        }
      }
      ParquetFile::appendOrThrow(unmodified_sequence_builder.Append(unmodified_sequence), "unmodified_sequence");
      ParquetFile::appendOrThrow(protein_accessions_builder.Append(joinProteinAccessions_(compound.protein_refs)), "protein_accessions");
    }

    for (const auto& protein : targeted_exp.proteins)
    {
      stats.proteins_total++;
      if (OpenMS::String(protein.id).hasPrefix("DECOY_"))
      {
        stats.proteins_decoy++;
      }
    }

    auto precursor_schema = OSWPrecursorSchema::schema();

    auto precursors_table = arrow::Table::Make(
      precursor_schema,
      {ParquetFile::finishArray(precursor_id_builder, "precursor_id"),
       ParquetFile::finishArray(precursor_mz_builder, "precursor_mz"),
       ParquetFile::finishArray(precursor_charge_builder, "charge"),
       ParquetFile::finishArray(library_rt_builder, "library_rt"),
       ParquetFile::finishArray(drift_time_builder, "library_drift_time"),
       ParquetFile::finishArray(decoy_builder, "decoy"),
       ParquetFile::finishArray(traml_id_builder, "traml_id"),
       ParquetFile::finishArray(modified_sequence_builder, "modified_sequence"),
       ParquetFile::finishArray(unmodified_sequence_builder, "unmodified_sequence"),
       ParquetFile::finishArray(protein_accessions_builder, "protein_accessions")});

    // Validate precursors table against registry schema
    auto prec_validation = ArrowSchemaValidation::validate(precursors_table, OSWPrecursorSchema::schema());
    if (!prec_validation.valid)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Precursors table schema validation failed: " + prec_validation.toString(), "");
    }

    ParquetFile::writeTable(precursors_table, library_dir + "/precursors.parquet");

    arrow::Int64Builder transition_id_builder;
    arrow::Int64Builder transition_precursor_id_builder;
    arrow::StringBuilder transition_traml_id_builder;
    arrow::DoubleBuilder product_mz_builder;
    arrow::Int32Builder fragment_charge_builder;
    arrow::StringBuilder fragment_type_builder;
    arrow::StringBuilder annotation_builder;
    arrow::Int32Builder ordinal_builder;
    arrow::BooleanBuilder detecting_builder;
    arrow::BooleanBuilder identifying_builder;
    arrow::BooleanBuilder quantifying_builder;
    arrow::DoubleBuilder transition_intensity_builder;
    arrow::BooleanBuilder transition_decoy_builder;

    int64_t transition_id = 1;
    for (const auto& transition : targeted_exp.transitions)
    {
      auto precursor_it = compound_to_precursor.find(transition.peptide_ref);
      if (precursor_it == compound_to_precursor.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown peptide_ref '" + String(transition.peptide_ref) + "'");
      }
      const int64_t precursor_ref = precursor_it->second;
      ParquetFile::appendOrThrow(transition_id_builder.Append(transition_id++), "transition_id");
      ParquetFile::appendOrThrow(transition_precursor_id_builder.Append(precursor_ref), "precursor_id");
      ParquetFile::appendOrThrow(transition_traml_id_builder.Append(transition.transition_name), "traml_id");
      ParquetFile::appendOrThrow(product_mz_builder.Append(transition.product_mz), "product_mz");
      ParquetFile::appendOrThrow(fragment_charge_builder.Append(static_cast<int32_t>(transition.fragment_charge)), "charge");
      ParquetFile::appendOrThrow(fragment_type_builder.Append(transition.getFragmentType()), "type");
      ParquetFile::appendOrThrow(annotation_builder.Append(transition.getAnnotation()), "annotation");
      ParquetFile::appendOrThrow(ordinal_builder.Append(transition.fragment_nr), "ordinal");
      ParquetFile::appendOrThrow(detecting_builder.Append(transition.isDetectingTransition()), "detecting");
      ParquetFile::appendOrThrow(identifying_builder.Append(transition.isIdentifyingTransition()), "identifying");
      ParquetFile::appendOrThrow(quantifying_builder.Append(transition.isQuantifyingTransition()), "quantifying");
      ParquetFile::appendOrThrow(transition_intensity_builder.Append(transition.library_intensity), "library_intensity");
      ParquetFile::appendOrThrow(transition_decoy_builder.Append(transition.getDecoy()), "decoy");

      stats.transitions_total++;
      const bool transition_decoy = transition.getDecoy();
      if (transition_decoy) stats.transitions_decoy++;
      if (transition_decoy)
      {
        stats.transition_charge_counts_decoy[transition.fragment_charge]++;
      }
      else
      {
        stats.transition_charge_counts_target[transition.fragment_charge]++;
      }
      std::string fragment_type = transition.getFragmentType();
      for (auto& c : fragment_type)
      {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      if (fragment_type == "b")
      {
        if (transition_decoy) stats.fragment_b_decoy++;
        else stats.fragment_b_target++;
      }
      else if (fragment_type == "y")
      {
        if (transition_decoy) stats.fragment_y_decoy++;
        else stats.fragment_y_target++;
      }
      else
      {
        if (transition_decoy) stats.fragment_other_decoy++;
        else stats.fragment_other_target++;
      }
    }

    auto transition_schema = OSWTransitionSchema::schema();

    auto transitions_table = arrow::Table::Make(
      transition_schema,
      {ParquetFile::finishArray(transition_id_builder, "transition_id"),
       ParquetFile::finishArray(transition_precursor_id_builder, "precursor_id"),
       ParquetFile::finishArray(transition_traml_id_builder, "traml_id"),
       ParquetFile::finishArray(product_mz_builder, "product_mz"),
       ParquetFile::finishArray(fragment_charge_builder, "charge"),
       ParquetFile::finishArray(fragment_type_builder, "type"),
       ParquetFile::finishArray(annotation_builder, "annotation"),
       ParquetFile::finishArray(ordinal_builder, "ordinal"),
       ParquetFile::finishArray(detecting_builder, "detecting"),
       ParquetFile::finishArray(identifying_builder, "identifying"),
       ParquetFile::finishArray(quantifying_builder, "quantifying"),
       ParquetFile::finishArray(transition_intensity_builder, "library_intensity"),
       ParquetFile::finishArray(transition_decoy_builder, "decoy")});

    // Validate transitions table against registry schema
    auto trans_validation = ArrowSchemaValidation::validate(transitions_table, OSWTransitionSchema::schema());
    if (!trans_validation.valid)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Transitions table schema validation failed: " + trans_validation.toString(), "");
    }

    ParquetFile::writeTable(transitions_table, library_dir + "/transitions.parquet");

    writeLibraryMetadata_(library_dir, library_name, stats);

    if (!output_is_dir)
    {
      // Instead of zipping the whole directory in one go, add each file to the
      // archive individually. This allows streaming large parquet files into the
      // archive without unzipping/rezipping everything. Use a staging archive
      // and rename into place atomically to avoid destroying the existing
      // archive on partial failures.
      const std::filesystem::path dirpath = std::filesystem::u8path(std::string(base_dir));
      const std::filesystem::path outpath = std::filesystem::u8path(std::string(oswpq_path));
      const String output_zip_abs = File::absolutePath(oswpq_path);
      const String staging_zip = output_zip_abs + ".tmp";

      if (File::exists(staging_zip))
      {
        File::remove(staging_zip);
      }

      for (auto it = std::filesystem::recursive_directory_iterator(dirpath); it != std::filesystem::recursive_directory_iterator(); ++it)
      {
        if (it->is_directory()) continue;
        const auto full = it->path();
        std::string rel = std::filesystem::relative(full, dirpath).generic_string();
        // Ensure forward slashes (zip expects '/'). generic_string() already uses '/'.
        ZipArchiveFile::addOrReplaceFromFile(staging_zip, String(rel), String(full.string()));
      }
      // After adding/replacing files in the staging archive, write a sidecar index
      // and then atomically move the staging archive into place.
      ZipArchiveFile::writeSidecarIndex(staging_zip);
      if (File::exists(output_zip_abs))
      {
        File::remove(output_zip_abs);
      }
      std::filesystem::rename(std::filesystem::u8path(std::string(staging_zip)),
                              std::filesystem::u8path(std::string(output_zip_abs)));
    }
  }
} // namespace OpenMS
