// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/FORMAT/ParquetFile.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#ifdef WITH_PARQUET
#include <arrow/api.h>
#endif

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <map>
#include <sstream>
#include <cctype>

namespace
{
  using OpenMS::Size;
#ifdef WITH_PARQUET

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
    std::string out;
    out.reserve(input.size());
    for (OpenMS::Size i = 0; i < input.size(); ++i)
    {
      const char c = input[i];
      switch (c)
      {
        case '\\': out += "\\\\"; break;
        case '"': out += "\\\""; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
          if (static_cast<unsigned char>(c) < 0x20)
          {
            const char hex[] = "0123456789abcdef";
            out += "\\u00";
            out += hex[(c >> 4) & 0x0f];
            out += hex[c & 0x0f];
          }
          else
          {
            out += c;
          }
          break;
      }
    }
    return out;
  }

  std::string jsonMap_(const std::map<int, Size>& values)
  {
    std::ostringstream os;
    os << "{";
    bool first = true;
    for (const auto& entry : values)
    {
      if (!first) os << ", ";
      first = false;
      os << "\"" << entry.first << "\": " << entry.second;
    }
    os << "}";
    return os.str();
  }

  std::string jsonMapByClass_(const std::map<int, Size>& target_values,
                              const std::map<int, Size>& decoy_values)
  {
    std::ostringstream os;
    os << "{\"target\": " << jsonMap_(target_values)
       << ", \"decoy\": " << jsonMap_(decoy_values) << "}";
    return os.str();
  }

  void writeLibraryMetadata_(const OpenMS::String& library_dir, const OpenMS::String& library_name,
                             const OpenMSLibraryStats& stats)
  {
    const OpenMS::String metadata_path = library_dir + "/metadata.json";
    std::ofstream out(metadata_path.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open())
    {
      throw OpenMS::Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, metadata_path);
    }

    const Size proteins_target = stats.proteins_total - stats.proteins_decoy;
    const Size peptides_target = stats.peptides_total - stats.peptides_decoy;
    const Size precursors_target = stats.precursors_total - stats.precursors_decoy;
    const Size compounds_target = stats.compounds_total - stats.compounds_decoy;
    const Size transitions_target = stats.transitions_total - stats.transitions_decoy;

    out << "{\n"
        << "  \"mzspec_lib\": {\n"
        << "    \"format_version\": \"1.0\",\n"
        << "    \"attributes\": [\n"
        << "      {\"accession\": \"MS:1003186\", \"name\": \"library format version\", \"value\": \"1.0\"},\n"
        << "      {\"accession\": \"MS:1003188\", \"name\": \"library name\", \"value\": \"" << jsonEscape_(library_name) << "\"},\n"
        << "      {\"accession\": \"MS:1003207\", \"name\": \"library creation software\", \"value\": \"OpenMS\"}\n"
        << "    ]\n"
        << "  },\n"
        << "  \"openms\": {\n"
        << "    \"schema_version\": 1,\n"
        << "    \"generator\": \"OpenMS TransitionParquetFile\",\n"
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
#endif // WITH_PARQUET
} // namespace

namespace OpenMS
{
  void TransitionParquetFile::convertParquetToTargetedExperiment(
    const String& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
#ifdef WITH_PARQUET
    std::unique_ptr<File::TempDir> temp_dir;
    const String base_dir = ParquetFile::unzipDirectory(oswpq_dir, temp_dir);
    const String library_dir = base_dir + "/library";
    const String precursors_path = library_dir + "/precursors.parquet";
    const String transitions_path = library_dir + "/transitions.parquet";

    if (!File::exists(precursors_path) || !File::exists(transitions_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing required parquet files in '" + oswpq_dir + "'");
    }

    auto precursors_table = ParquetFile::readTable(precursors_path);
    auto transitions_table = ParquetFile::readTable(transitions_path);

    auto precursor_id_col = ParquetFile::getColumn(precursors_table, "precursor_id");
    auto precursor_mz_col = ParquetFile::getColumn(precursors_table, "precursor_mz");
    auto charge_col = ParquetFile::getColumn(precursors_table, "charge");
    auto library_rt_col = ParquetFile::getColumn(precursors_table, "library_rt");
    auto drift_time_col = ParquetFile::getOptionalColumn(precursors_table, "library_drift_time");
    auto traml_id_col = ParquetFile::getOptionalColumn(precursors_table, "traml_id");
    auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, "decoy");
    auto modified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, "modified_sequence");
    auto unmodified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, "unmodified_sequence");
    auto protein_accessions_col = ParquetFile::getOptionalColumn(precursors_table, "protein_accessions");

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
      compound.id = precursor_id_str;
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

    auto transition_id_col = ParquetFile::getColumn(transitions_table, "transition_id");
    auto transition_traml_id_col = ParquetFile::getOptionalColumn(transitions_table, "traml_id");
    auto transition_precursor_id_col = ParquetFile::getColumn(transitions_table, "precursor_id");
    auto product_mz_col = ParquetFile::getColumn(transitions_table, "product_mz");
    auto fragment_charge_col = ParquetFile::getColumn(transitions_table, "charge");
    auto fragment_type_col = ParquetFile::getColumn(transitions_table, "type");
    auto fragment_annotation_col = ParquetFile::getOptionalColumn(transitions_table, "annotation");
    auto fragment_ordinal_col = ParquetFile::getColumn(transitions_table, "ordinal");
    auto detecting_col = ParquetFile::getColumn(transitions_table, "detecting");
    auto identifying_col = ParquetFile::getColumn(transitions_table, "identifying");
    auto quantifying_col = ParquetFile::getColumn(transitions_table, "quantifying");
    auto transition_intensity_col = ParquetFile::getColumn(transitions_table, "library_intensity");
    auto transition_decoy_col = ParquetFile::getColumn(transitions_table, "decoy");

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
      transition.peptide_ref = String(precursor_id);
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
#else
    (void)oswpq_dir;
    (void)targeted_exp;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }

  void TransitionParquetFile::convertLightTargetedExperimentToParquet(
    const String& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const
  {
#ifdef WITH_PARQUET
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
    for (const auto& compound : targeted_exp.compounds)
    {
      if (compound_to_precursor.find(compound.id) != compound_to_precursor.end())
      {
        continue;
      }

      int64_t precursor_id = 0;
      try
      {
        precursor_id = OpenMS::String(compound.id).toInt64();
      }
      catch (OpenMS::Exception::ConversionError&)
      {
        precursor_id = next_precursor_id++;
      }

      if (precursor_id >= next_precursor_id)
      {
        next_precursor_id = precursor_id + 1;
      }

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
      ParquetFile::appendOrThrow(precursor_mz_builder.Append(precursor_mz[compound.id]), "precursor_mz");
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

    auto precursor_schema = arrow::schema({
      arrow::field("precursor_id", arrow::int64()),
      arrow::field("precursor_mz", arrow::float64()),
      arrow::field("charge", arrow::int32()),
      arrow::field("library_rt", arrow::float64()),
      arrow::field("library_drift_time", arrow::float64()),
      arrow::field("decoy", arrow::boolean()),
      arrow::field("traml_id", arrow::utf8()),
      arrow::field("modified_sequence", arrow::utf8()),
      arrow::field("unmodified_sequence", arrow::utf8()),
      arrow::field("protein_accessions", arrow::utf8())
    });

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
      const int64_t precursor_ref = compound_to_precursor[transition.peptide_ref];
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
      for (auto& c : fragment_type) c = static_cast<char>(std::tolower(c));
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

    auto transition_schema = arrow::schema({
      arrow::field("transition_id", arrow::int64()),
      arrow::field("precursor_id", arrow::int64()),
      arrow::field("traml_id", arrow::utf8()),
      arrow::field("product_mz", arrow::float64()),
      arrow::field("charge", arrow::int32()),
      arrow::field("type", arrow::utf8()),
      arrow::field("annotation", arrow::utf8()),
      arrow::field("ordinal", arrow::int32()),
      arrow::field("detecting", arrow::boolean()),
      arrow::field("identifying", arrow::boolean()),
      arrow::field("quantifying", arrow::boolean()),
      arrow::field("library_intensity", arrow::float64()),
      arrow::field("decoy", arrow::boolean())
    });

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

    ParquetFile::writeTable(transitions_table, library_dir + "/transitions.parquet");

    writeLibraryMetadata_(library_dir, library_name, stats);

    if (!output_is_dir)
    {
      ParquetFile::zipDirectory(base_dir, oswpq_path);
    }
#else
    (void)oswpq_path;
    (void)targeted_exp;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }
} // namespace OpenMS
