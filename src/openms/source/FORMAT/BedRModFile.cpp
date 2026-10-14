// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/BedRModFile.h>

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/ModomicsJSONDataProvider.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <limits>
#include <optional>
#include <set>
#include <tuple>
#include <vector>

namespace OpenMS
{
  namespace
  {
    struct BedRow
    {
      std::string chrom;
      Size chrom_start{0};
      Size chrom_end{0};
      const Ribonucleotide* ribonucleotide{nullptr};
      std::string modification_name;
      std::string modomics_short_name;
      double score{std::numeric_limits<double>::infinity()};
      Size target_mapping_count{0};
      std::set<std::pair<Size, Size>> fragments;
      std::set<std::pair<std::string, std::string>> observations;
    };

    struct ModificationMapping
    {
      Int chebi_id{0};
      std::string modomics_short_name;
    };

    const std::map<std::string, Ribonucleotide>& modomicsCodes_()
    {
      static const auto codes = []
      {
        std::map<std::string, Ribonucleotide> result;
        ModomicsJSONDataProvider provider("CHEMISTRY/Modomics.json");
        for (const auto& entry : provider.loadRibonucleotides())
        {
          result.emplace(entry.ribo->getCode(), *entry.ribo);
        }
        return result;
      }();
      return codes;
    }

    std::string modomicsName_(const Ribonucleotide& ribo, const ModificationMapping& mapping)
    {
      const auto& codes = modomicsCodes_();
      const std::string name = mapping.modomics_short_name.empty() ? ribo.getCode() : mapping.modomics_short_name;
      auto known = codes.find(name);
      if (known == codes.end() && mapping.modomics_short_name.empty())
      {
        known = codes.find(ribo.getNewCode());
      }
      if (known == codes.end() || known->first.find_first_of(",:\t\r\n") != std::string::npos ||
          known->second.getOrigin() != ribo.getOrigin() || known->second.getFormula() != ribo.getFormula())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Cannot represent RNA modification '" + ribo.getCode() +
          "' as a bedRMod MODOMICS short name. Supply a chemically compatible, registered "
          "modomics_short_name in the mapping CSV; comma-containing dictionary names are unsupported.");
      }
      return known->first;
    }

    std::string normalizeHeader_(std::string value)
    {
      value = StringUtils::toLowered(StringUtils::trimmed(value));
      std::replace(value.begin(), value.end(), ' ', '_');
      return value;
    }

    bool toInt_(const std::string& value, Int& result)
    {
      try
      {
        result = StringUtils::toInt32(StringUtils::trimmed(value));
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
      return true;
    }

    std::map<std::string, ModificationMapping> readChebiMapping_(const std::string& chebi_mapping_file)
    {
      std::map<std::string, ModificationMapping> mapping;
      if (chebi_mapping_file.empty())
      {
        return mapping;
      }

      std::string full_path = File::find(chebi_mapping_file);
      TextFile input(full_path, true, -1, true, "");
      std::vector<std::string> lines(input.begin(), input.end());
      if (lines.empty())
      {
        return mapping;
      }

      std::vector<std::string> header;
      if (!StringUtils::split(lines[0], ',', header, true))
      {
        return mapping;
      }

      Size mod_col = Size(-1);
      Size chebi_col = Size(-1);
      Size modomics_col = Size(-1);
      for (Size i = 0; i < header.size(); ++i)
      {
        std::string key = normalizeHeader_(header[i]);
        if ((key == "mod") || (key == "name"))
        {
          mod_col = i;
        }
        else if (key == "chebi_id")
        {
          chebi_col = i;
        }
        else if (key == "modomics_short_name")
        {
          modomics_col = i;
        }
      }

      if ((mod_col == Size(-1)) || ((chebi_col == Size(-1)) && (modomics_col == Size(-1))))
      {
        OPENMS_LOG_WARN << "Warning: ChEBI mapping file '" << full_path
                        << "' requires 'mod'/'name' and 'chebi_id' or 'modomics_short_name' columns."
                        << std::endl;
        return mapping;
      }

      for (Size row = 1; row < lines.size(); ++row)
      {
        if (StringUtils::trimmed(lines[row]).empty())
        {
          continue;
        }
        std::vector<std::string> values;
        StringUtils::split(lines[row], ',', values, true);
        if (mod_col >= values.size())
        {
          continue;
        }

        std::string mod = StringUtils::trimmed(values[mod_col]);
        if (mod.empty())
        {
          continue;
        }

        ModificationMapping entry;
        if (chebi_col < values.size())
        {
          toInt_(values[chebi_col], entry.chebi_id);
        }
        if (modomics_col < values.size())
        {
          entry.modomics_short_name = StringUtils::trimmed(values[modomics_col]);
        }
        mapping[mod] = std::move(entry);
      }

      return mapping;
    }

    using SequenceKey = std::pair<UInt32, std::string>; // database and accession

    /// The database sequences of @p run that are decoys; evidence for other (or uncatalogued) sequences counts as target.
    std::set<SequenceKey> decoySequences_(const IdentificationData::Run& run)
    {
      std::set<SequenceKey> decoys;
      if (run.getDatabaseSequences())
      {
        for (const auto& sequence : *run.getDatabaseSequences())
        {
          if (sequence.target_decoy == IdentificationData::TargetDecoy::DECOY) decoys.emplace(sequence.database.value, sequence.accession);
        }
      }
      return decoys;
    }

    std::string joinFragmentPositions_(const std::set<std::pair<Size, Size>>& fragments,
                                       const bool use_start)
    {
      std::string result;
      for (const auto& fragment : fragments)
      {
        if (!result.empty()) result += ",";
        result += std::to_string((use_start ? fragment.first : fragment.second) + 1);
      }
      return result;
    }

  }

  void BedRModFile::store(const std::string& out_file,
                          const IdentificationData& id_data,
                          const std::string& chebi_mapping_file,
                          const std::map<std::string, std::string>& reference_metadata) const
  {
    const auto metadata = [&](const std::string& key)
    {
      const auto value = reference_metadata.find(key);
      if (value == reference_metadata.end() || value->second.empty()) return std::string("NA");
      if (value->second.find_first_of("\r\n") != std::string::npos)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "bedRMod reference metadata must fit on one line: " + key);
      }
      return value->second;
    };
    // Validate all header values before opening (and potentially truncating) output.
    const std::string organism = metadata("organism");
    const std::string assembly = metadata("assembly");
    const std::string annotation_source = metadata("annotation_source");
    const std::string annotation_version = metadata("annotation_version");
    const std::string sequencing_platform = metadata("sequencing_platform");
    const std::string basecalling = metadata("basecalling");
    const std::string bioinformatics_workflow = metadata("bioinformatics_workflow");
    const std::string experiment = metadata("experiment");
    const std::string external_source = metadata("external_source");
    const auto chebi_mapping = readChebiMapping_(chebi_mapping_file);

    // Keep aggregation private until another consumer needs this representation.
    // Coverage counts distinct observations of the modification's canonical base.
    // An observation is a spectrum (data ID) of a file; spectra of an unknown file are told apart by their source.
    using SiteKey = std::tuple<std::string, Size, std::string>;
    using ObservationKey = std::pair<std::string, std::string>;
    std::map<SiteKey, BedRow> rows;
    std::map<std::tuple<std::string, Size, char>, std::set<ObservationKey>> obs_per_position;
    Size missing_qvalues = 0;

    for (const auto& run : id_data.getRuns())
    {
      if (run.getMoleculeKind() != IdentificationData::MoleculeKind::OLIGONUCLEOTIDE)
      {
        continue;
      }
      std::optional<IdentificationData::ScoreId> qvalue_score;
      for (Size i = 0; i < run.getScoreDefinitions().size(); ++i)
      {
        if (run.getScoreDefinitions()[i].name == "PSM-level q-value")
        {
          qvalue_score = run.getScoreId(UInt32(i));
          break;
        }
      }
      const auto decoys = decoySequences_(run);
      for (Size source_index = 0; source_index < run.getSources().size(); ++source_index)
      {
        const auto& source = run.getSources()[source_index];
        const std::string file = source.file.path.empty() ? "run " + run.getUuid() + " source " + std::to_string(source_index) : source.file.path;
        for (const auto& query : source.identifications)
        {
          const ObservationKey observation {file, query.data_id};
          for (const auto& match : query.getMatches())
          {
            if (match.sequence_evidence.empty())
            {
              continue;
            }

            // Only output matches for which a q-value has been determined (i.e. FDR
            // was run). A value of -1 is the sentinel for a match without an
            // assigned q-value, so it must not contribute rows or frequencies.
            double score = std::numeric_limits<double>::quiet_NaN();
            if (qvalue_score)
            {
              const auto qval = run.getScore(match.getId(), *qvalue_score);
              if (qval && std::isfinite(*qval) && *qval >= 0.0)
              {
                score = *qval;
              }
            }
            if (std::isnan(score))
            {
              ++missing_qvalues;
              continue; // no valid q-value — skip this match
            }
            const NASequence sequence = NASequence::fromString(match.representation);
            // Only count mappings to target (non-decoy) sequences
            const Size target_mapping_count = static_cast<Size>(std::count_if(match.sequence_evidence.begin(), match.sequence_evidence.end(),
              [&](const auto& evidence) { return ! decoys.contains({evidence.database.value, evidence.accession}); }));

            for (const auto& evidence : match.sequence_evidence)
            {
              if (decoys.contains({evidence.database.value, evidence.accession})) continue;
              if (! evidence.start || ! evidence.end || *evidence.end < *evidence.start)
              {
                continue;
              }
              const std::string& chrom = evidence.accession;
              const Size start_pos = *evidence.start, end_pos = *evidence.end;
              for (Size i = 0; i < sequence.size(); ++i)
              {
                const auto* ribo = sequence[i];
                if (ribo == nullptr || i > end_pos - start_pos)
                {
                  continue;
                }
                const Size position = start_pos + i;
                obs_per_position[{chrom, position, ribo->getOrigin()}].insert(observation);

                // Canonical residues contribute coverage above, but do not describe
                // positive modification calls. Terminal groups are also excluded.
                if (!ribo->isModified() || ribo->getTermSpecificity() != Ribonucleotide::ANYWHERE)
                {
                  continue;
                }

                const std::string mod_code = ribo->getCode();
                const auto mapping = chebi_mapping.find(mod_code);
                const ModificationMapping entry = mapping == chebi_mapping.end() ? ModificationMapping{} : mapping->second;
                const std::string short_name = modomicsName_(*ribo, entry);
                // Explicit aliases of the same MODOMICS modification share a site.
                BedRow& row = rows[{chrom, position, short_name}];
                row.chrom = chrom;
                row.chrom_start = position;
                row.chrom_end = row.chrom_start + 1;
                row.ribonucleotide = ribo;
                row.score = std::min(row.score, score); // best supporting q-value; lower is better
                row.target_mapping_count = std::max(row.target_mapping_count, target_mapping_count);
                row.fragments.emplace(start_pos, end_pos);
                row.observations.insert(observation);

                row.modomics_short_name = short_name;
                if (entry.chebi_id > 0)
                {
                  const std::string name = std::to_string(entry.chebi_id);
                  if (!row.modification_name.empty() && row.modification_name != short_name && row.modification_name != name)
                  {
                    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "Conflicting ChEBI identifiers for MODOMICS modification '" + short_name + "'");
                  }
                  row.modification_name = name;
                }
                else if (row.modification_name.empty())
                {
                  row.modification_name = short_name;
                }
              }
            }
          }
        }
      }
    }

    if (missing_qvalues > 0)
    {
      OPENMS_LOG_WARN << "bedRMod: skipped " << missing_qvalues
                      << " RNA matches without an assigned PSM-level q-value; run FDR before export.\n";
    }

    std::vector<std::string> modification_names;
    std::set<std::string> seen_mod_names;
    std::map<std::string, std::string> dictionary;
    for (const auto& [site, row] : rows)
    {
      if (!row.ribonucleotide)
      {
        continue;  // Skip rows with null ribonucleotide
      }
      // Use actual origin from ribonucleotide instead of inferring from mod code
      const std::string base_origin(1, row.ribonucleotide->getOrigin());
      std::string mod_name = row.modification_name + ":" + row.modomics_short_name + ":" + base_origin;
      const auto [name, inserted] = dictionary.emplace(row.modification_name, mod_name);
      if (!inserted && name->second != mod_name)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Conflicting bedRMod dictionary entries for identifier '" + row.modification_name + "'");
      }
      if (seen_mod_names.insert(mod_name).second)
      {
        modification_names.push_back(mod_name);
      }
    }

    std::ofstream out(out_file.c_str());
    if (!out.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file);
    }

    out << "#fileformat=bedRModv2\n";
    out << "#organism=" << organism << "\n";
    out << "#modification_type=RNA\n";
    out << "#modification_names=";
    for (Size i = 0; i < modification_names.size(); ++i)
    {
      if (i != 0) out << ",";
      out << modification_names[i];
    }
    out << "\n";
    out << "#assembly=" << assembly << "\n";
    out << "#annotation_source=" << annotation_source << "\n";
    out << "#annotation_version=" << annotation_version << "\n";
    out << "#sequencing_platform=" << sequencing_platform << '\n';
    out << "#basecalling=" << basecalling << '\n';
    out << "#bioinformatics_workflow=" << bioinformatics_workflow << '\n';
    out << "#experiment=" << experiment << '\n';
    out << "#external_source=" << external_source << '\n';
    out << "#score=minimum supporting PSM-level q-value (lower is better)\n";
    out << "#coverage=distinct observations of the canonical base at the site\n";
    out << "#frequency=100 * distinct modified observations / coverage\n";
    out << "#chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb coverage frequency max_target_mapping_count frag_start frag_end\n";

    for (const auto& [site, row] : rows)
    {
      if (!row.ribonucleotide)
      {
        continue;  // Skip rows with null ribonucleotide
      }
      
      // Calculate frequency: count of this modification / total observations at this position
      const Size total = obs_per_position.at({row.chrom, row.chrom_start, row.ribonucleotide->getOrigin()}).size();
      const double frequency = 100.0 * row.observations.size() / total;

      out << row.chrom << '\t'
          << row.chrom_start << '\t'
          << row.chrom_end << '\t'
          << row.modification_name << '\t'
          << row.score << '\t'
          << "+" << '\t'
          << row.chrom_start << '\t'
          << row.chrom_end << '\t'
          << "0,0,0" << '\t'
          << total << '\t'
          << frequency << '\t'
          << row.target_mapping_count << '\t'
          << joinFragmentPositions_(row.fragments, true) << '\t'
          << joinFragmentPositions_(row.fragments, false) << '\n';
    }
    // Buffered writes may fail after a successful open (e.g. a full device).
    out.flush();
    if (!out)
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file);
    }
  }
}
