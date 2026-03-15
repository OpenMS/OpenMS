// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/BedRModFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
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
      String chrom;
      Int chrom_start{0};
      Int chrom_end{0};
      String mod;
      Int chebi_id{0};
      Int score{0};
      Int coverage{1};

      bool operator<(const BedRow& rhs) const
      {
        return std::tie(chrom, chrom_start, chrom_end, mod, score, coverage) <
               std::tie(rhs.chrom, rhs.chrom_start, rhs.chrom_end, rhs.mod, rhs.score, rhs.coverage);
      }
    };

    String normalizeHeader_(String value)
    {
      value = value.trim().toLower();
      value.substitute(" ", "_");
      return value;
    }

    String inferBase_(const String& mod)
    {
      if (mod.empty())
      {
        return "N";
      }
      const char last = mod.back();
      if ((last == 'A') || (last == 'C') || (last == 'G') || (last == 'U') ||
          (last == 'I') || (last == 'Y'))
      {
        return String(last);
      }
      return "N";
    }

    bool toInt_(const String& value, Int& result)
    {
      try
      {
        String tmp = value;
        tmp.trim();
        result = tmp.toInt();
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
      return true;
    }

    std::map<String, Int> readChebiMapping_(const String& chebi_mapping_file)
    {
      std::map<String, Int> mapping;
      if (chebi_mapping_file.empty())
      {
        return mapping;
      }

      String full_path = File::find(chebi_mapping_file);
      TextFile input(full_path, true, -1, true, "");
      StringList lines(input.begin(), input.end());
      if (lines.empty())
      {
        return mapping;
      }

      StringList header;
      if (!lines[0].split(',', header, true))
      {
        return mapping;
      }

      Size mod_col = Size(-1);
      Size chebi_col = Size(-1);
      for (Size i = 0; i < header.size(); ++i)
      {
        String key = normalizeHeader_(header[i]);
        if ((key == "mod") || (key == "name"))
        {
          mod_col = i;
        }
        else if (key == "chebi_id")
        {
          chebi_col = i;
        }
      }

      if ((mod_col == Size(-1)) || (chebi_col == Size(-1)))
      {
        OPENMS_LOG_WARN << "Warning: ChEBI mapping file '" << full_path
                        << "' misses required columns ('mod'/'name' and 'chebi_id'/'chebi id')."
                        << std::endl;
        return mapping;
      }

      for (Size row = 1; row < lines.size(); ++row)
      {
        if (lines[row].trim().empty())
        {
          continue;
        }
        StringList values;
        lines[row].split(',', values, true);
        if ((mod_col >= values.size()) || (chebi_col >= values.size()))
        {
          continue;
        }

        String mod = values[mod_col].trim();
        if (mod.empty())
        {
          continue;
        }

        Int chebi = 0;
        if (!toInt_(values[chebi_col], chebi))
        {
          continue;
        }
        mapping[mod] = chebi;
      }

      return mapping;
    }

    Int getRoundedScore_(const IdentificationData::ObservationMatch& match,
                         const IdentificationData::ScoreTypeRef& score_ref,
                         const IdentificationData::ScoreTypes& all_score_types)
    {
      double score = 0.0;
      bool found = false;
      if (score_ref != all_score_types.end())
      {
        std::tie(score, found) = match.getScore(score_ref);
      }
      if (!found)
      {
        std::optional<IdentificationData::ScoreTypeRef> any_ref;
        std::tie(score, any_ref, found) = match.getMostRecentScore();
      }
      if (!found || !std::isfinite(score))
      {
        return 0;
      }
      return Int(std::lround(score));
    }

    Int getCoverage_(const IdentificationData::ObservationRef& observation_ref)
    {
      if (!observation_ref->metaValueExists("precursor_intensity"))
      {
        return 1;
      }
      double intensity = 1.0;
      try
      {
        intensity = static_cast<double>(observation_ref->getMetaValue("precursor_intensity"));
      }
      catch (Exception::ConversionError&)
      {
        return 1;
      }

      if (!std::isfinite(intensity))
      {
        return 1;
      }
      Int cov = Int(std::lround(intensity));
      return std::max<Int>(1, cov);
    }
  }

  void BedRModFile::store(const String& out_file,
                          const IdentificationData& id_data,
                          const String& chebi_mapping_file) const
  {
    const auto chebi_mapping = readChebiMapping_(chebi_mapping_file);
    const auto hyperscore_ref = id_data.findScoreType("hyperscore");
    const auto& score_types = id_data.getScoreTypes();

    std::vector<BedRow> rows;
    std::set<String> missing_mods;

    for (const IdentificationData::ObservationMatch& match : id_data.getObservationMatches())
    {
      if (match.identified_molecule_var.getMoleculeType() != IdentificationData::MoleculeType::RNA)
      {
        continue;
      }
      auto oligo_ref = match.identified_molecule_var.getIdentifiedOligoRef();
      const auto& oligo = *oligo_ref;
      const auto& sequence = oligo.sequence;

      if (oligo.parent_matches.empty())
      {
        continue;
      }

      const Int score = getRoundedScore_(match, hyperscore_ref, score_types);
      const Int coverage = getCoverage_(match.observation_ref);

      for (const auto& parent_pair : oligo.parent_matches)
      {
        const String& chrom = parent_pair.first->accession;
        for (const auto& parent_match : parent_pair.second)
        {
          if (!parent_match.hasValidPositions())
          {
            continue;
          }
          for (Size i = 0; i < sequence.size(); ++i)
          {
            const auto* ribo = sequence[i];
            if ((ribo == nullptr) || !ribo->isModified())
            {
              continue;
            }

            BedRow row;
            row.chrom = chrom;
            row.chrom_start = Int(parent_match.start_pos + i);
            row.chrom_end = row.chrom_start + 1;
            row.mod = ribo->getCode();
            row.score = score;
            row.coverage = coverage;

            auto pos = chebi_mapping.find(row.mod);
            if (pos != chebi_mapping.end())
            {
              row.chebi_id = pos->second;
            }
            else
            {
              row.chebi_id = 0;
              missing_mods.insert(row.mod);
            }

            if (row.chrom_end > row.chrom_start)
            {
              rows.push_back(std::move(row));
            }
          }
        }
      }
    }

    if (!missing_mods.empty())
    {
      OPENMS_LOG_WARN << "Warning: Missing ChEBI ids for modifications (using 0): ";
      bool first = true;
      for (const auto& mod : missing_mods)
      {
        if (!first) OPENMS_LOG_WARN << ", ";
        OPENMS_LOG_WARN << mod;
        first = false;
      }
      OPENMS_LOG_WARN << std::endl;
    }

    std::sort(rows.begin(), rows.end());

    std::vector<String> modification_names;
    std::set<String> seen_mod_names;
    for (const auto& row : rows)
    {
      String mod_name = String(row.chebi_id) + ":" + row.mod + ":" + inferBase_(row.mod);
      if (seen_mod_names.insert(mod_name).second)
      {
        modification_names.push_back(mod_name);
      }
    }

    std::ofstream out(std::string(out_file).c_str());
    if (!out.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out_file);
    }

    out << "#fileformat=bedRModv2\n";
    out << "#organism=9606\n";
    out << "#modification_type=RNA\n";
    out << "#modification_names=";
    for (Size i = 0; i < modification_names.size(); ++i)
    {
      if (i != 0) out << ",";
      out << modification_names[i];
    }
    out << "\n";
    out << "#assembly=GRCh38\n";
    out << "#annotation_source=gtrnadb\n";
    out << "#annotation_version=2.0\n";
    out << "#sequencing_platform=ddMS2\n";
    out << "#basecalling=NASE\n";
    out << "#bioinformatics_workflow=NA\n";
    out << "#experiment=NA\n";
    out << "#external_source=NA\n";
    out << "#chrom chromStart chromEnd name score strand thickStart thickEnd itemRgb coverage frequency\n";

    for (const auto& row : rows)
    {
      out << row.chrom << '\t'
          << row.chrom_start << '\t'
          << row.chrom_end << '\t'
          << row.chebi_id << '\t'
          << row.score << '\t'
          << "+" << '\t'
          << row.chrom_start << '\t'
          << row.chrom_end << '\t'
          << "0,0,0" << '\t'
          << row.coverage << '\t'
          << 100 << '\n';
    }
  }
}
