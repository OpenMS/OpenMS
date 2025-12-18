// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/StreamingTSVReader.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <sstream>
#include <algorithm>

namespace OpenMS
{

  const std::vector<std::string> StreamingTSVReader::header_names_ = {
    "PrecursorMz", "ProductMz", "PrecursorCharge", "ProductCharge",
    "LibraryIntensity", "NormalizedRetentionTime", "PeptideSequence",
    "ModifiedPeptideSequence", "ProteinId", "ProteinName", "GeneName",
    "FragmentType", "FragmentSeriesNumber", "Annotation", "CollisionEnergy",
    "PrecursorIonMobility", "TransitionGroupId", "TransitionId", "Decoy",
    "DetectingTransition", "IdentifyingTransition", "QuantifyingTransition",
    "Peptidoforms", "CompoundName", "SMILES", "SumFormula", "Adducts",
    "RetentionTime", "Tr_recalibrated", "iRT", "FullPeptideName",
    "FullUniModPeptideName", "UniprotId", "FragmentCharge", "PeptideGroupLabel",
    "LabelType", "Charge", "Sequence", "transition_name", "transition_group_id"
  };

  StreamingTSVReader::StreamingTSVReader() :
    ProgressLogger()
  {
  }

  StreamingTSVReader::~StreamingTSVReader() = default;

  void StreamingTSVReader::parseHeader_(const std::string& line, char& delimiter,
                                         std::map<std::string, int>& header_dict) const
  {
    // Determine delimiter (tab or comma)
    delimiter = (line.find('\t') != std::string::npos) ? '\t' : ',';

    std::stringstream ss(line);
    std::string field;
    int col_idx = 0;

    while (std::getline(ss, field, delimiter))
    {
      // Remove quotes if present
      if (!field.empty() && field.front() == '"' && field.back() == '"')
      {
        field = field.substr(1, field.size() - 2);
      }
      header_dict[field] = col_idx++;
    }
  }

  void StreamingTSVReader::parseLine_(const std::string& line, char delimiter,
                                       const std::map<std::string, int>& header_dict,
                                       OpenSwath::LightTransition& transition,
                                       std::string& group_id,
                                       std::string& sequence,
                                       std::string& full_peptide_name,
                                       int& charge,
                                       std::vector<std::string>& protein_names,
                                       std::string& compound_name,
                                       double& rt) const
  {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, delimiter))
    {
      // Remove quotes if present
      if (!field.empty() && field.front() == '"' && field.back() == '"')
      {
        field = field.substr(1, field.size() - 2);
      }
      fields.push_back(field);
    }

    // Helper to get field value by name
    auto getField = [&](const std::string& name) -> std::string {
      auto it = header_dict.find(name);
      if (it != header_dict.end() && it->second < (int)fields.size())
      {
        return fields[it->second];
      }
      return "";
    };

    auto getFieldDouble = [&](const std::string& name, double default_val = 0.0) -> double {
      std::string val = getField(name);
      if (val.empty() || val == "NA") return default_val;
      try { return std::stod(val); }
      catch (...) { return default_val; }
    };

    auto getFieldInt = [&](const std::string& name, int default_val = 0) -> int {
      std::string val = getField(name);
      if (val.empty() || val == "NA") return default_val;
      try { return std::stoi(val); }
      catch (...) { return default_val; }
    };

    auto getFieldBool = [&](const std::string& name, bool default_val = false) -> bool {
      std::string val = getField(name);
      if (val.empty()) return default_val;
      return val == "1" || val == "TRUE" || val == "true" || val == "True";
    };

    // Parse transition fields
    transition.precursor_mz = getFieldDouble("PrecursorMz");
    transition.product_mz = getFieldDouble("ProductMz");
    if (transition.product_mz == 0.0)
    {
      transition.product_mz = getFieldDouble("FragmentMz");
    }
    transition.library_intensity = getFieldDouble("LibraryIntensity");
    if (transition.library_intensity == 0.0)
    {
      transition.library_intensity = getFieldDouble("RelativeFragmentIntensity");
    }

    // Transition identifiers
    transition.transition_name = getField("TransitionId");
    if (transition.transition_name.empty())
    {
      transition.transition_name = getField("TransitionName");
    }
    if (transition.transition_name.empty())
    {
      transition.transition_name = getField("transition_name");
    }

    // Group ID
    group_id = getField("TransitionGroupId");
    if (group_id.empty())
    {
      group_id = getField("TransitionGroupName");
    }
    if (group_id.empty())
    {
      group_id = getField("transition_group_id");
    }
    transition.peptide_ref = group_id;

    // Charge
    transition.fragment_charge = getFieldInt("ProductCharge");
    if (transition.fragment_charge == 0)
    {
      transition.fragment_charge = getFieldInt("FragmentCharge");
    }

    // Flags
    transition.decoy = getFieldBool("Decoy") || getFieldBool("IsDecoy") || getFieldBool("decoy");
    transition.detecting_transition = getFieldBool("DetectingTransition", true);
    transition.quantifying_transition = getFieldBool("QuantifyingTransition", true);
    transition.identifying_transition = getFieldBool("IdentifyingTransition", false);

    // Ion mobility
    transition.precursor_im = getFieldDouble("PrecursorIonMobility", -1);

    // Peptide/compound data
    sequence = getField("PeptideSequence");
    if (sequence.empty())
    {
      sequence = getField("Sequence");
    }
    if (sequence.empty())
    {
      sequence = getField("StrippedSequence");
    }

    full_peptide_name = getField("ModifiedPeptideSequence");
    if (full_peptide_name.empty())
    {
      full_peptide_name = getField("FullPeptideName");
    }
    if (full_peptide_name.empty())
    {
      full_peptide_name = getField("FullUniModPeptideName");
    }
    if (full_peptide_name.empty())
    {
      full_peptide_name = getField("ModifiedSequence");
    }

    charge = getFieldInt("PrecursorCharge");
    if (charge == 0)
    {
      charge = getFieldInt("Charge");
    }

    // Proteins
    protein_names.clear();
    std::string protein_str = getField("ProteinId");
    if (protein_str.empty())
    {
      protein_str = getField("ProteinName");
    }
    if (!protein_str.empty())
    {
      // Split by semicolon or slash
      std::stringstream pss(protein_str);
      std::string prot;
      char pdelim = (protein_str.find(';') != std::string::npos) ? ';' : '/';
      while (std::getline(pss, prot, pdelim))
      {
        if (!prot.empty())
        {
          protein_names.push_back(prot);
        }
      }
    }

    // Compound name (for metabolomics)
    compound_name = getField("CompoundName");
    if (compound_name.empty())
    {
      compound_name = getField("CompoundId");
    }

    // Retention time
    rt = getFieldDouble("NormalizedRetentionTime");
    if (rt == 0.0)
    {
      rt = getFieldDouble("RetentionTime");
    }
    if (rt == 0.0)
    {
      rt = getFieldDouble("Tr_recalibrated");
    }
    if (rt == 0.0)
    {
      rt = getFieldDouble("iRT");
    }
  }

  std::string StreamingTSVReader::buildModifiedSequence_(const OpenSwath::LightCompound& compound)
  {
    if (compound.modifications.empty())
    {
      return compound.sequence;
    }

    // Sort modifications by location
    std::vector<OpenSwath::LightModification> sorted_mods = compound.modifications;
    std::sort(sorted_mods.begin(), sorted_mods.end(),
              [](const auto& a, const auto& b) { return a.location < b.location; });

    std::string result;
    int seq_len = static_cast<int>(compound.sequence.size());

    // N-terminal modifications (location -1)
    for (const auto& mod : sorted_mods)
    {
      if (mod.location == -1)
      {
        result += "(UniMod:" + std::to_string(mod.unimod_id) + ")";
      }
    }

    // Sequence with internal modifications
    for (int i = 0; i < seq_len; ++i)
    {
      result += compound.sequence[i];
      for (const auto& mod : sorted_mods)
      {
        if (mod.location == i)
        {
          result += "(UniMod:" + std::to_string(mod.unimod_id) + ")";
        }
      }
    }

    // C-terminal modifications (location == sequence length)
    for (const auto& mod : sorted_mods)
    {
      if (mod.location == seq_len)
      {
        result += "(UniMod:" + std::to_string(mod.unimod_id) + ")";
      }
    }

    return result;
  }

  std::string StreamingTSVReader::buildCollisionKey(const OpenSwath::LightCompound& compound)
  {
    return buildModifiedSequence_(compound) + "+" + std::to_string(compound.charge);
  }

  StreamingTSVReader::IndexResult StreamingTSVReader::buildIndex(const String& filename)
  {
    IndexResult result;

    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // Read header
    std::string line;
    if (!std::getline(ifs, line))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   "Empty TSV file", filename);
    }

    char delimiter;
    std::map<std::string, int> header_dict;
    parseHeader_(line, delimiter, header_dict);

    // Track group ordering
    std::string current_group;
    std::string previous_group;
    std::unordered_set<std::string> seen_groups;

    // Temporary storage for current compound
    OpenSwath::LightCompound current_compound;
    std::vector<std::string> current_protein_ids;

    OPENMS_LOG_INFO << "Pass 1: Building index and validating group order..." << std::endl;

    size_t line_num = 1;
    while (std::getline(ifs, line))
    {
      ++line_num;
      if (line.empty()) continue;

      OpenSwath::LightTransition transition;
      std::string group_id, sequence, full_peptide_name, compound_name;
      int charge;
      double rt;
      std::vector<std::string> protein_names;

      parseLine_(line, delimiter, header_dict, transition, group_id,
                 sequence, full_peptide_name, charge, protein_names, compound_name, rt);

      ++result.total_transitions;

      // Check group ordering
      if (group_id != current_group)
      {
        // Finish previous group if exists
        if (!current_group.empty())
        {
          // Add collision key for completed compound
          std::string collision_key = buildCollisionKey(current_compound);
          result.collision_index[collision_key] = current_compound.id;
        }

        // Check if we've seen this group before (split group detection)
        if (seen_groups.count(group_id) > 0)
        {
          ifs.close();
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Streaming mode requires all transitions of a peptide group to appear "
            "consecutively in the TSV file. Found split group: '" + group_id + "' at line " +
            std::to_string(line_num) + ". Pre-sort input with: "
            "sort -t$'\\t' -k<group_col> -s input.tsv > sorted.tsv",
            filename);
        }

        seen_groups.insert(current_group);
        previous_group = current_group;
        current_group = group_id;
        ++result.total_groups;

        // Start new compound
        current_compound = OpenSwath::LightCompound();
        current_compound.id = group_id;
        current_compound.sequence = sequence;
        current_compound.charge = charge;
        current_compound.rt = rt;
        current_compound.compound_name = compound_name;
        current_protein_ids = protein_names;

        // Parse modifications from full peptide name if available
        // (simplified - in production would use AASequence parser)

        // Add to compound map if not exists
        if (result.compound_map.find(group_id) == result.compound_map.end())
        {
          result.compound_map[group_id] = result.compounds.size();
          current_compound.protein_refs = protein_names;
          result.compounds.push_back(current_compound);
        }

        // Add proteins
        for (const auto& prot_id : protein_names)
        {
          if (result.protein_map.find(prot_id) == result.protein_map.end())
          {
            result.protein_map[prot_id] = result.proteins.size();
            OpenSwath::LightProtein protein;
            protein.id = prot_id;
            result.proteins.push_back(protein);
          }
        }
      }

      if (line_num % 1000000 == 0)
      {
        OPENMS_LOG_INFO << "  Processed " << line_num << " lines, "
                        << result.total_groups << " groups..." << std::endl;
      }
    }

    // Handle last group
    if (!current_group.empty())
    {
      std::string collision_key = buildCollisionKey(current_compound);
      result.collision_index[collision_key] = current_compound.id;
    }

    ifs.close();

    OPENMS_LOG_INFO << "Pass 1 complete: " << result.total_transitions << " transitions, "
                    << result.total_groups << " groups, "
                    << result.collision_index.size() << " collision keys" << std::endl;

    return result;
  }

  void StreamingTSVReader::processGroups(
    const String& filename,
    const IndexResult& index,
    GroupCallback callback)
  {
    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // Read header
    std::string line;
    std::getline(ifs, line);

    char delimiter;
    std::map<std::string, int> header_dict;
    parseHeader_(line, delimiter, header_dict);

    std::string current_group;
    std::vector<OpenSwath::LightTransition> current_transitions;
    std::vector<std::string> current_protein_ids;
    OpenSwath::LightCompound current_compound;

    OPENMS_LOG_INFO << "Pass 2: Processing groups..." << std::endl;

    size_t groups_processed = 0;

    auto flushGroup = [&]() {
      if (!current_group.empty() && !current_transitions.empty())
      {
        // Look up compound from index
        auto comp_it = index.compound_map.find(current_group);
        if (comp_it != index.compound_map.end())
        {
          current_compound = index.compounds[comp_it->second];
        }

        callback(current_group, current_transitions, current_compound, current_protein_ids);

        ++groups_processed;
        if (groups_processed % 10000 == 0)
        {
          OPENMS_LOG_INFO << "  Processed " << groups_processed << " groups..." << std::endl;
        }
      }
      current_transitions.clear();
      current_protein_ids.clear();
    };

    while (std::getline(ifs, line))
    {
      if (line.empty()) continue;

      OpenSwath::LightTransition transition;
      std::string group_id, sequence, full_peptide_name, compound_name;
      int charge;
      double rt;
      std::vector<std::string> protein_names;

      parseLine_(line, delimiter, header_dict, transition, group_id,
                 sequence, full_peptide_name, charge, protein_names, compound_name, rt);

      if (group_id != current_group)
      {
        flushGroup();
        current_group = group_id;
        current_protein_ids = protein_names;
      }

      current_transitions.push_back(std::move(transition));
    }

    // Flush last group
    flushGroup();

    ifs.close();

    OPENMS_LOG_INFO << "Pass 2 complete: processed " << groups_processed << " groups" << std::endl;
  }

} // namespace OpenMS
