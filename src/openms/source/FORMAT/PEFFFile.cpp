// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PEFFFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>

#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>

namespace OpenMS
{
  using namespace std;

  namespace
  {
    /// Sentinel value for "annotation ID not set"
    constexpr UInt ANNOT_ID_NOT_SET = std::numeric_limits<UInt>::max();

    /// Extract optional annotation ID prefix from a field (format: "id:value" or just "value")
    /// Returns pair of (annotation_id, remaining_value). annotation_id is ANNOT_ID_NOT_SET if not present.
    std::pair<UInt, String> extractAnnotationId(const String& field)
    {
      size_t colon_pos = field.find(':');
      if (colon_pos != String::npos)
      {
        String id_str = field.substr(0, colon_pos);
        // Check if the prefix is a non-negative integer (annotation ID)
        bool is_numeric = !id_str.empty();
        for (size_t i = 0; i < id_str.size() && is_numeric; ++i)
        {
          if (!std::isdigit(id_str[i])) is_numeric = false;
        }
        if (is_numeric)
        {
          return {static_cast<UInt>(id_str.toInt()), field.substr(colon_pos + 1)};
        }
      }
      return {ANNOT_ID_NOT_SET, field};
    }

    /// Split a string by pipe character, respecting backslash escapes.
    /// Handles \| (escaped pipe), \\ (escaped backslash).
    /// Unescapes the result strings.
    std::vector<String> splitByPipeEscapeAware(const String& str)
    {
      std::vector<String> result;
      String current;
      for (size_t i = 0; i < str.size(); ++i)
      {
        if (str[i] == '\\' && i + 1 < str.size())
        {
          char next = str[i + 1];
          if (next == '|' || next == '\\' || next == '(' || next == ')')
          {
            current += next;
            ++i; // skip escaped char
            continue;
          }
        }
        if (str[i] == '|')
        {
          result.push_back(current);
          current.clear();
        }
        else
        {
          current += str[i];
        }
      }
      result.push_back(current);
      return result;
    }
    /// Escape PEFF reserved characters: \, |, (, )
    /// This is the inverse of the unescape done in splitByPipeEscapeAware.
    String escape_peff(const String& str)
    {
      String result;
      result.reserve(str.size());
      for (size_t i = 0; i < str.size(); ++i)
      {
        char c = str[i];
        if (c == '\\' || c == '|' || c == '(' || c == ')')
        {
          result += '\\';
        }
        result += c;
      }
      return result;
    }

    bool isXsdInteger_(const String& s)
    {
      if (s.empty()) return false;
      size_t i = 0;
      if (s[0] == '+' || s[0] == '-')
      {
        i = 1;
      }
      if (i >= s.size()) return false;
      for (; i < s.size(); ++i)
      {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
      }
      return true;
    }

    bool isUnsignedInteger_(const String& s)
    {
      if (s.empty()) return false;
      for (Size i = 0; i < s.size(); ++i)
      {
        if (!std::isdigit(static_cast<unsigned char>(s[i]))) return false;
      }
      return true;
    }

    bool isXsdDecimal_(const String& s)
    {
      if (s.empty()) return false;
      size_t i = 0;
      if (s[0] == '+' || s[0] == '-')
      {
        i = 1;
      }
      if (i >= s.size()) return false;

      bool saw_digit = false;
      bool saw_dot = false;
      for (; i < s.size(); ++i)
      {
        const char c = s[i];
        if (std::isdigit(static_cast<unsigned char>(c)))
        {
          saw_digit = true;
          continue;
        }
        if (c == '.' && !saw_dot)
        {
          saw_dot = true;
          continue;
        }
        return false;
      }
      return saw_digit;
    }

    bool isXsdBoolean_(const String& s)
    {
      String lower(s);
      lower.toLower();
      return lower == "true" || lower == "false" || lower == "1" || lower == "0";
    }

    bool isXsdDate_(const String& s)
    {
      // xsd:date is ISO 8601-ish. We accept YYYY-MM-DD with optional timezone.
      // We deliberately do not validate month/day ranges here.
      if (s.size() < 10) return false;
      size_t i = 0;
      if (s[0] == '-')
      {
        i = 1;
      }
      if (i + 10 > s.size()) return false;

      auto is_digit = [&s](const size_t p) -> bool
      {
        return std::isdigit(static_cast<unsigned char>(s[p])) != 0;
      };

      if (!(is_digit(i + 0) && is_digit(i + 1) && is_digit(i + 2) && is_digit(i + 3))) return false;
      if (s[i + 4] != '-') return false;
      if (!(is_digit(i + 5) && is_digit(i + 6))) return false;
      if (s[i + 7] != '-') return false;
      if (!(is_digit(i + 8) && is_digit(i + 9))) return false;

      // Optional timezone
      if (i + 10 == s.size()) return true;
      const size_t tz = i + 10;
      if (s[tz] == 'Z') return (tz + 1 == s.size());
      if (s[tz] != '+' && s[tz] != '-') return false;
      if (tz + 6 != s.size()) return false;
      if (!(is_digit(tz + 1) && is_digit(tz + 2))) return false;
      if (s[tz + 3] != ':') return false;
      if (!(is_digit(tz + 4) && is_digit(tz + 5))) return false;
      return true;
    }

    bool isXsdTime_(const String& s)
    {
      // xsd:time is ISO 8601-ish. We accept hh:mm:ss(.fff)? with optional timezone.
      // We deliberately do not validate hour/min/sec ranges here.
      if (s.size() < 8) return false;

      auto is_digit = [&s](const size_t p) -> bool
      {
        return std::isdigit(static_cast<unsigned char>(s[p])) != 0;
      };

      if (!(is_digit(0) && is_digit(1))) return false;
      if (s[2] != ':') return false;
      if (!(is_digit(3) && is_digit(4))) return false;
      if (s[5] != ':') return false;
      if (!(is_digit(6) && is_digit(7))) return false;

      size_t pos = 8;
      // Optional fractional seconds
      if (pos < s.size() && s[pos] == '.')
      {
        ++pos;
        if (pos >= s.size()) return false;
        bool saw_digit = false;
        while (pos < s.size() && std::isdigit(static_cast<unsigned char>(s[pos])))
        {
          saw_digit = true;
          ++pos;
        }
        if (!saw_digit) return false;
      }

      if (pos == s.size()) return true;

      // Optional timezone
      if (s[pos] == 'Z') return (pos + 1 == s.size());
      if (s[pos] != '+' && s[pos] != '-') return false;
      if (pos + 6 != s.size()) return false;
      if (!(is_digit(pos + 1) && is_digit(pos + 2))) return false;
      if (s[pos + 3] != ':') return false;
      if (!(is_digit(pos + 4) && is_digit(pos + 5))) return false;
      return true;
    }

    bool validateCustomKeyFieldType_(const String& field_type, const String& value)
    {
      String type(field_type);
      type.trim();
      const String type_lower = type.toLower();

      if (type_lower.empty()) return true;
      if (type_lower == "string") return true;
      if (type_lower == "integer") return isXsdInteger_(value);
      if (type_lower == "decimal") return isXsdDecimal_(value);
      if (type_lower == "boolean") return isXsdBoolean_(value);
      if (type_lower == "date") return isXsdDate_(value);
      if (type_lower == "time") return isXsdTime_(value);

      if (type_lower.hasPrefix("enumeration(") && type.hasSuffix(")"))
      {
        const size_t lparen = type.find('(');
        const size_t rparen = type.rfind(')');
        if (lparen == String::npos || rparen == String::npos || rparen <= lparen + 1) return false;

        const String inside = type.substr(lparen + 1, rparen - lparen - 1);
        std::vector<String> options;
        inside.split('|', options);
        for (String opt : options)
        {
          opt.trim();
          if (opt == value) return true;
        }
        return false;
      }

      // Unknown type: do not validate.
      return true;
    }

    /// Split a CustomKeyDef "(defKey=value|...)" string into defKey=value segments.
    /// Respects quoted values and nested parentheses (e.g., enumeration(...) in FieldTypes).
    std::vector<String> splitCustomKeyDefSegments_(const String& str)
    {
      std::vector<String> result;
      String current;
      bool in_quotes = false;
      int paren_depth = 0;

      for (size_t i = 0; i < str.size(); ++i)
      {
        const char c = str[i];

        // Keep escapes verbatim, but ensure we don't treat escaped characters as syntax.
        if (c == '\\' && i + 1 < str.size())
        {
          current += c;
          current += str[i + 1];
          ++i;
          continue;
        }

        if (c == '"')
        {
          in_quotes = !in_quotes;
          current += c;
          continue;
        }

        if (!in_quotes)
        {
          if (c == '(')
          {
            ++paren_depth;
          }
          else if (c == ')' && paren_depth > 0)
          {
            --paren_depth;
          }
          else if (c == '|' && paren_depth == 0)
          {
            result.push_back(current);
            current.clear();
            continue;
          }
        }

        current += c;
      }

      result.push_back(current);
      return result;
    }

    /// Unquote a CustomKeyDef string value like "text", unescaping \" inside.
    String unquoteCustomKeyDefValue_(const String& value)
    {
      String trimmed(value);
      trimmed.trim();
      if (trimmed.size() >= 2 && trimmed[0] == '"' && trimmed[trimmed.size() - 1] == '"')
      {
        String inner = trimmed.substr(1, trimmed.size() - 2);
        String result;
        result.reserve(inner.size());
        for (size_t i = 0; i < inner.size(); ++i)
        {
          if (inner[i] == '\\' && i + 1 < inner.size() && inner[i + 1] == '"')
          {
            result += '"';
            ++i;
          }
          else
          {
            result += inner[i];
          }
        }
        return result;
      }
      return trimmed;
    }

    std::vector<String> splitCommaSeparatedList_(const String& value)
    {
      std::vector<String> items;
      std::vector<String> parts;
      value.split(',', parts);
      for (String& p : parts)
      {
        p.trim();
        if (!p.empty())
        {
          items.push_back(p);
        }
      }
      return items;
    }

    PEFFCustomKeyDef parseCustomKeyDef_(const String& value)
    {
      PEFFCustomKeyDef ckd;
      String trimmed(value);
      trimmed.trim();

      // Spec format: CustomKeyDef=(defKey=value|defKey=value|...)
      if (trimmed.size() >= 2 && trimmed[0] == '(' && trimmed[trimmed.size() - 1] == ')')
      {
        String inner = trimmed.substr(1, trimmed.size() - 2);
        std::vector<String> segments = splitCustomKeyDefSegments_(inner);
        for (String& segment : segments)
        {
          segment.trim();
          if (segment.empty())
          {
            continue;
          }

          size_t eq_pos = segment.find('=');
          if (eq_pos == String::npos)
          {
            OPENMS_LOG_WARN << "PEFF: CustomKeyDef segment is missing '=': '" << segment << "'.\n";
            continue;
          }

          String def_key = segment.substr(0, eq_pos);
          def_key.trim();
          String def_value = segment.substr(eq_pos + 1);
          def_value.trim();

          if (def_key == "KeyName")
          {
            ckd.key_name = def_value;
          }
          else if (def_key == "Description")
          {
            ckd.description = unquoteCustomKeyDefValue_(def_value);
          }
          else if (def_key == "ConceptCURIE")
          {
            ckd.concept_curie = def_value;
          }
          else if (def_key == "RegExp")
          {
            ckd.regexp = unquoteCustomKeyDefValue_(def_value);
          }
          else if (def_key == "FieldNames")
          {
            ckd.field_names = splitCommaSeparatedList_(def_value);
          }
          else if (def_key == "FieldTypes")
          {
            ckd.field_types = splitCommaSeparatedList_(def_value);
          }
          else
          {
            OPENMS_LOG_WARN << "PEFF: Unknown CustomKeyDef attribute '" << def_key << "'.\n";
          }
        }
      }
      else
      {
        // Legacy/simple format (historical): CustomKeyDef=KeyName:description
        size_t colon_pos = trimmed.find(':');
        if (colon_pos != String::npos)
        {
          ckd.key_name = trimmed.substr(0, colon_pos);
          ckd.description = trimmed.substr(colon_pos + 1);
        }
        else
        {
          ckd.key_name = trimmed;
        }
      }

      if (ckd.key_name.empty())
      {
        OPENMS_LOG_WARN << "PEFF: CustomKeyDef is missing required KeyName.\n";
      }
      return ckd;
    }

    String escapeCustomKeyDefQuotedValue_(const String& value)
    {
      String result;
      result.reserve(value.size());
      for (size_t i = 0; i < value.size(); ++i)
      {
        if (value[i] == '"')
        {
          result += '\\';
        }
        result += value[i];
      }
      return result;
    }
  } // anonymous namespace

  FASTAFile::FASTAEntry PEFFEntry::toFASTAEntry() const
  {
    // Build description with protein names
    String desc;
    if (!protein_names.empty())
    {
      desc = protein_names[0];
    }
    return FASTAFile::FASTAEntry(identifier, desc, sequence);
  }

  PEFFEntry PEFFEntry::fromFASTAEntry(const FASTAFile::FASTAEntry& fasta)
  {
    PEFFEntry entry;
    entry.identifier = fasta.identifier;
    entry.sequence = fasta.sequence;
    if (!fasta.description.empty())
    {
      entry.protein_names.push_back(fasta.description);
    }
    entry.sequence_length = fasta.sequence.size();
    return entry;
  }

  AASequence PEFFEntry::getSequence() const
  {
    return AASequence::fromString(sequence);
  }

  AASequence PEFFEntry::getModifiedSequence() const
  {
    // Start with the base sequence
    AASequence seq = AASequence::fromString(sequence);

    // Apply modifications
    for (const auto& mod : modifications)
    {
      // Skip modifications with unknown position
      if (mod.position == 0)
      {
        continue;
      }

      // Convert 1-based position to 0-based index
      Size idx = mod.position - 1;
      if (idx >= seq.size())
      {
        continue; // Skip if position is out of range
      }

      try
      {
        // Try to find the modification - accession first, then name
        const ResidueModification* res_mod = nullptr;
        String residue = seq[idx].getOneLetterCode();

        // Try accession first (UNIMOD:xx, MOD:xxxxx, etc.)
        if (!mod.accession.empty())
        {
          try
          {
            res_mod = ModificationsDB::getInstance()->getModification(mod.accession, residue, ResidueModification::ANYWHERE);
          }
          catch (const Exception::BaseException&)
          {
            res_mod = nullptr; // Fall through to try name
          }
        }

        // Fall back to name if accession lookup failed
        if (res_mod == nullptr && !mod.name.empty())
        {
          res_mod = ModificationsDB::getInstance()->getModification(mod.name, residue, ResidueModification::ANYWHERE);
        }

        if (res_mod != nullptr)
        {
          seq.setModification(idx, res_mod->getFullId());
        }
        else
        {
          // Neither accession nor name lookup succeeded
          OPENMS_LOG_WARN << "Could not resolve modification: " << mod.accession << " (" << mod.name << ")" << "\n";
        }
      }
      catch (const Exception::BaseException&)
      {
        // Skip if modification cannot be resolved (exception during lookup)
        OPENMS_LOG_WARN << "Could not resolve modification: " << mod.accession << " (" << mod.name << ")" << "\n";
      }
    }

    return seq;
  }

  void PEFFEntry::getVariantSequences(
    std::vector<std::string>& descriptions,
    std::vector<AASequence>& sequences,
    bool include_complex) const
  {
    descriptions.clear();
    sequences.clear();

    // Simple variants (single AA substitution)
    for (const auto& var : simple_variants)
    {
      if (var.position == 0 || var.position > sequence.size() || var.variant_aa == '\0')
      {
        continue;
      }

      // Create a copy of the sequence with the variant
      String variant_seq = sequence;
      variant_seq[var.position - 1] = var.variant_aa; // Convert 1-based to 0-based

      // Create description
      String desc = String(sequence[var.position - 1]) + String(var.position) + String(var.variant_aa);
      if (!var.optional_tag.empty())
      {
        desc += " (" + var.optional_tag + ")";
      }

      try
      {
        descriptions.push_back(desc);
        sequences.push_back(AASequence::fromString(variant_seq));
      }
      catch (const Exception::BaseException&)
      {
        // Skip if sequence cannot be parsed
      }
    }

    // Complex variants (insertion, deletion, multi-residue substitution)
    if (include_complex)
    {
      for (const auto& var : complex_variants)
      {
        if (var.start_position == 0 || var.end_position > sequence.size() ||
            var.start_position > var.end_position)
        {
          continue;
        }

        // Replace region [start, end] with replacement sequence
        String variant_seq = sequence.substr(0, var.start_position - 1)
                           + var.replacement
                           + sequence.substr(var.end_position);

        // Create description
        String desc = String(var.start_position) + "-" + String(var.end_position) + ">" + var.replacement;
        if (!var.optional_tag.empty())
        {
          desc += " (" + var.optional_tag + ")";
        }

        try
        {
          descriptions.push_back(desc);
          sequences.push_back(AASequence::fromString(variant_seq));
        }
        catch (const Exception::BaseException&)
        {
          // Skip if sequence cannot be parsed
        }
      }
    }
  }

  AASequence PEFFEntry::getProcessedSequence(const String& region_accession) const
  {
    // Find the first processed region of the given accession
    for (const auto& region : processed_regions)
    {
      if (region.accession == region_accession)
      {
        // For signal peptides (PEFF:0001021), return the mature protein (after the signal)
        if (region_accession == "PEFF:0001021" || region.name == "signal peptide")
        {
          if (region.end_position < sequence.size())
          {
            String mature_seq = sequence.substr(region.end_position);
            return AASequence::fromString(mature_seq);
          }
        }
        // For other region types, return the region itself
        else if (region.start_position > 0 && region.end_position >= region.start_position &&
                 region.end_position <= sequence.size())
        {
          String region_seq = sequence.substr(region.start_position - 1, region.end_position - region.start_position + 1);
          return AASequence::fromString(region_seq);
        }
      }
    }

    // Return empty sequence if region not found
    return AASequence();
  }

  void PEFFEntry::digestWithVariants(
    const ProteaseDigestion& digestor,
    std::vector<std::string>& descriptions,
    std::vector<AASequence>& sequences,
    Size min_length,
    Size max_length,
    bool include_reference,
    bool include_variants,
    bool include_modifications) const
  {
    descriptions.clear();
    sequences.clear();

    if (sequence.empty())
    {
      return;
    }

    // Step 1: Digest the reference sequence to get peptide boundaries
    AASequence protein = AASequence::fromString(sequence);
    std::vector<std::pair<size_t, size_t>> peptide_ranges;
    digestor.digest(protein, peptide_ranges, min_length, max_length);

    // Step 2: For each peptide, find variants/modifications and generate combinations
    for (const auto& [start, end] : peptide_ranges)
    {
      // Collect simple variants within this peptide range (if requested)
      // PEFF positions are 1-based, ranges are 0-based
      std::vector<const PEFFVariantSimple*> local_variants;
      if (include_variants)
      {
        for (const auto& var : simple_variants)
        {
          // Skip stop codons and invalid positions
          if (var.variant_aa == '*' || var.variant_aa == '\0' || var.position == 0)
          {
            continue;
          }

          // Convert 1-based PEFF position to 0-based
          size_t var_pos_0based = var.position - 1;
          if (var_pos_0based >= start && var_pos_0based < end)
          {
            local_variants.push_back(&var);
          }
        }
      }

      // Collect modifications within this peptide range (if requested)
      std::vector<const PEFFModification*> local_mods;
      if (include_modifications)
      {
        for (const auto& mod : modifications)
        {
          // Skip modifications with unknown position
          if (mod.position == 0)
          {
            continue;
          }

          // Convert 1-based PEFF position to 0-based
          size_t mod_pos_0based = mod.position - 1;
          if (mod_pos_0based >= start && mod_pos_0based < end)
          {
            local_mods.push_back(&mod);
          }
        }
      }

      // Extract reference peptide sequence
      String ref_peptide = sequence.substr(start, end - start);

      // Total number of combinatorial elements
      size_t num_variants = local_variants.size();
      size_t num_mods = local_mods.size();
      size_t total_elements = num_variants + num_mods;

      // Safeguard against combinatorial explosion
      // Limit to 20 elements (2^20 = ~1 million combinations per peptide)
      constexpr size_t MAX_COMBINATORIAL_ELEMENTS = 20;
      if (total_elements > MAX_COMBINATORIAL_ELEMENTS)
      {
        OPENMS_LOG_WARN << "Peptide at position " << start << "-" << end
                        << " has " << total_elements << " variants/modifications. "
                        << "Limiting to first " << MAX_COMBINATORIAL_ELEMENTS
                        << " to avoid combinatorial explosion." << "\n";
        // Prioritize variants over modifications, truncate excess
        if (num_variants > MAX_COMBINATORIAL_ELEMENTS)
        {
          local_variants.resize(MAX_COMBINATORIAL_ELEMENTS);
          local_mods.clear();
        }
        else
        {
          local_mods.resize(MAX_COMBINATORIAL_ELEMENTS - num_variants);
        }
        num_variants = local_variants.size();
        num_mods = local_mods.size();
        total_elements = num_variants + num_mods;
      }

      // Include reference peptide if requested (combo = 0)
      if (include_reference)
      {
        try
        {
          descriptions.push_back("");
          sequences.push_back(AASequence::fromString(ref_peptide));
        }
        catch (const Exception::BaseException&)
        {
          // Skip if sequence cannot be parsed
        }
      }

      // Generate all 2^(k+m) combinations for k variants and m modifications
      if (total_elements > 0)
      {
        size_t num_combinations = 1ULL << total_elements;

        // Start from 1 to skip the all-reference combination (already added above if include_reference)
        for (size_t combo = 1; combo < num_combinations; ++combo)
        {
          String variant_peptide = ref_peptide;
          String description;

          // Apply variants (bits 0 to num_variants-1), skip if multiple at same position
          std::set<size_t> modified_positions;
          bool has_conflict = false;
          for (size_t i = 0; i < num_variants && !has_conflict; ++i)
          {
            if (combo & (1ULL << i))  // Apply this variant
            {
              const PEFFVariantSimple* var = local_variants[i];
              size_t local_pos = var->position - 1 - start;  // Position within peptide

              // Check for conflicting variants at same position
              if (modified_positions.count(local_pos))
              {
                has_conflict = true;
                break;
              }
              modified_positions.insert(local_pos);

              // Build description
              if (!description.empty())
              {
                description += ";";
              }
              description += String(ref_peptide[local_pos]) + String(var->position) + String(var->variant_aa);
              if (!var->optional_tag.empty())
              {
                description += "(" + var->optional_tag + ")";
              }

              // Apply the substitution
              variant_peptide[local_pos] = var->variant_aa;
            }
          }

          // Skip combinations with conflicting variants
          if (has_conflict)
          {
            continue;
          }

          // Convert to AASequence for modification application
          AASequence pep_seq;
          try
          {
            pep_seq = AASequence::fromString(variant_peptide);
          }
          catch (const Exception::BaseException&)
          {
            continue;  // Skip if sequence cannot be parsed
          }

          // Apply modifications (bits num_variants to num_variants+num_mods-1)
          bool mod_failed = false;
          for (size_t i = 0; i < num_mods; ++i)
          {
            if (combo & (1ULL << (num_variants + i)))  // Apply this modification
            {
              const PEFFModification* mod = local_mods[i];
              size_t local_pos = mod->position - 1 - start;  // Position within peptide

              // Build description
              if (!description.empty())
              {
                description += ";";
              }
              if (!mod->accession.empty())
              {
                description += mod->accession + "@" + String(mod->position);
              }
              else if (!mod->name.empty())
              {
                description += mod->name + "@" + String(mod->position);
              }

              // Try to apply the modification
              try
              {
                // Try accession first (UNIMOD:xx, MOD:xxxxx, etc.)
                const ResidueModification* res_mod = nullptr;
                String residue = pep_seq[local_pos].getOneLetterCode();

                if (!mod->accession.empty())
                {
                  try
                  {
                    res_mod = ModificationsDB::getInstance()->getModification(mod->accession, residue, ResidueModification::ANYWHERE);
                  }
                  catch (const Exception::BaseException&)
                  {
                    res_mod = nullptr;
                  }
                }

                // Fall back to name
                if (res_mod == nullptr && !mod->name.empty())
                {
                  try
                  {
                    res_mod = ModificationsDB::getInstance()->getModification(mod->name, residue, ResidueModification::ANYWHERE);
                  }
                  catch (const Exception::BaseException&)
                  {
                    res_mod = nullptr;
                  }
                }

                if (res_mod != nullptr)
                {
                  pep_seq.setModification(local_pos, res_mod->getFullId());
                }
                else
                {
                  mod_failed = true;  // Skip this combination if mod not found
                }
              }
              catch (const Exception::BaseException&)
              {
                mod_failed = true;
              }
            }
          }

          if (!mod_failed)
          {
            descriptions.push_back(description);
            sequences.push_back(pep_seq);
          }
        }
      }
    }
  }

  std::vector<std::pair<String, AASequence>> PEFFEntry::enumeratePEFFModifications_(
    const AASequence& peptide,
    const std::vector<std::pair<Size, const PEFFModification*>>& peff_mods,
    const String& base_description)
  {
    std::vector<std::pair<String, AASequence>> result;

    if (peff_mods.empty())
    {
      result.emplace_back(base_description, peptide);
      return result;
    }

    // Limit combinatorial explosion
    constexpr size_t MAX_MODS = 20;
    size_t num_mods = std::min(peff_mods.size(), MAX_MODS);
    if (peff_mods.size() > MAX_MODS)
    {
      OPENMS_LOG_WARN << "Peptide has " << peff_mods.size() << " PEFF modifications. "
                      << "Limiting to first " << MAX_MODS << " to avoid combinatorial explosion." << "\n";
    }

    size_t num_combinations = 1ULL << num_mods;

    for (size_t combo = 0; combo < num_combinations; ++combo)
    {
      AASequence mod_peptide = peptide;
      String description = base_description;
      bool mod_failed = false;

      for (size_t i = 0; i < num_mods && !mod_failed; ++i)
      {
        if (combo & (1ULL << i))
        {
          const auto& [local_pos, mod] = peff_mods[i];

          // Build description
          if (!description.empty())
          {
            description += ";";
          }
          if (!mod->accession.empty())
          {
            description += mod->accession + "@" + String(local_pos + 1);  // 1-based for display
          }
          else if (!mod->name.empty())
          {
            description += mod->name + "@" + String(local_pos + 1);
          }

          // Apply modification
          try
          {
            const ResidueModification* res_mod = nullptr;
            String residue = mod_peptide[local_pos].getOneLetterCode();

            if (!mod->accession.empty())
            {
              try
              {
                res_mod = ModificationsDB::getInstance()->getModification(mod->accession, residue, ResidueModification::ANYWHERE);
              }
              catch (const Exception::BaseException&) { res_mod = nullptr; }
            }

            if (res_mod == nullptr && !mod->name.empty())
            {
              try
              {
                res_mod = ModificationsDB::getInstance()->getModification(mod->name, residue, ResidueModification::ANYWHERE);
              }
              catch (const Exception::BaseException&) { res_mod = nullptr; }
            }

            if (res_mod != nullptr)
            {
              mod_peptide.setModification(local_pos, res_mod->getFullId());
            }
            else
            {
              mod_failed = true;
            }
          }
          catch (const Exception::BaseException&)
          {
            mod_failed = true;
          }
        }
      }

      if (!mod_failed)
      {
        result.emplace_back(description, mod_peptide);
      }
    }

    return result;
  }

  void PEFFEntry::generatePeptides(
    const ProteaseDigestion& digestor,
    std::vector<std::string>& descriptions,
    std::vector<AASequence>& sequences,
    const std::vector<std::string>& fixed_mods,
    const std::vector<std::string>& variable_mods,
    Size max_variable_mods_per_peptide,
    Size min_length,
    Size max_length,
    bool include_reference,
    bool include_peff_variants,
    bool include_peff_modifications) const
  {
    descriptions.clear();
    sequences.clear();

    if (sequence.empty())
    {
      return;
    }

    // Step 1: Digest the reference sequence to get peptide boundaries
    AASequence protein = AASequence::fromString(sequence);
    std::vector<std::pair<size_t, size_t>> peptide_ranges;
    digestor.digest(protein, peptide_ranges, min_length, max_length);

    // Prepare fixed modifications lookup (applied to all compatible residues)
    ModifiedPeptideGenerator::MapToResidueType fixed_mods_map;
    if (!fixed_mods.empty())
    {
      StringList mod_list(fixed_mods.begin(), fixed_mods.end());
      fixed_mods_map = ModifiedPeptideGenerator::getModifications(mod_list);
    }

    // Prepare variable modifications lookup (for combinatorial generation)
    ModifiedPeptideGenerator::MapToResidueType var_mods_map;
    if (!variable_mods.empty())
    {
      StringList mod_list(variable_mods.begin(), variable_mods.end());
      var_mods_map = ModifiedPeptideGenerator::getModifications(mod_list);
    }

    // Step 2: Process each peptide
    for (const auto& [start, end] : peptide_ranges)
    {
      // Collect PEFF variants within this peptide range
      std::vector<const PEFFVariantSimple*> local_variants;
      if (include_peff_variants)
      {
        for (const auto& var : simple_variants)
        {
          if (var.variant_aa == '*' || var.variant_aa == '\0' || var.position == 0)
          {
            continue;
          }
          size_t var_pos_0based = var.position - 1;
          if (var_pos_0based >= start && var_pos_0based < end)
          {
            local_variants.push_back(&var);
          }
        }
      }

      // Collect PEFF modifications within this peptide range
      std::vector<std::pair<Size, const PEFFModification*>> local_mods;
      if (include_peff_modifications)
      {
        for (const auto& mod : modifications)
        {
          if (mod.position == 0) continue;
          size_t mod_pos_0based = mod.position - 1;
          if (mod_pos_0based >= start && mod_pos_0based < end)
          {
            local_mods.emplace_back(mod_pos_0based - start, &mod);  // Store local position
          }
        }
      }

      // Limit combined variant + modification combinations to avoid explosion
      // (variants and mods each contribute 2^n combinations, combined is 2^(v+m))
      constexpr size_t MAX_COMBINATORIAL_ELEMENTS = 20;
      size_t total_elements = local_variants.size() + local_mods.size();
      if (total_elements > MAX_COMBINATORIAL_ELEMENTS)
      {
        OPENMS_LOG_WARN << "Peptide at " << start << "-" << end << " has " << total_elements
                        << " variants+modifications. Limiting to " << MAX_COMBINATORIAL_ELEMENTS
                        << " to avoid combinatorial explosion." << "\n";
        // Prioritize variants over PEFF modifications
        if (local_variants.size() > MAX_COMBINATORIAL_ELEMENTS)
        {
          local_variants.resize(MAX_COMBINATORIAL_ELEMENTS);
          local_mods.clear();
        }
        else
        {
          local_mods.resize(MAX_COMBINATORIAL_ELEMENTS - local_variants.size());
        }
      }

      String ref_peptide_str = sequence.substr(start, end - start);

      // Step 3: Generate all variant combinations (2^k for k variants)
      size_t num_variant_combos = 1ULL << local_variants.size();

      for (size_t var_combo = 0; var_combo < num_variant_combos; ++var_combo)
      {
        // Skip non-reference if no variants selected and we don't want reference
        if (var_combo == 0 && !include_reference && local_variants.empty() && local_mods.empty())
        {
          continue;
        }

        String variant_peptide_str = ref_peptide_str;
        String variant_desc;

        // Apply selected variants (skip if multiple variants at same position)
        std::set<size_t> modified_positions;
        bool has_conflict = false;
        for (size_t i = 0; i < local_variants.size() && !has_conflict; ++i)
        {
          if (var_combo & (1ULL << i))
          {
            const PEFFVariantSimple* var = local_variants[i];
            size_t local_pos = var->position - 1 - start;

            // Check for conflicting variants at same position
            if (modified_positions.count(local_pos))
            {
              has_conflict = true;
              break;
            }
            modified_positions.insert(local_pos);

            if (!variant_desc.empty()) variant_desc += ";";
            variant_desc += String(ref_peptide_str[local_pos]) + String(var->position) + String(var->variant_aa);
            if (!var->optional_tag.empty())
            {
              variant_desc += "(" + var->optional_tag + ")";
            }

            variant_peptide_str[local_pos] = var->variant_aa;
          }
        }

        // Skip combinations with conflicting variants
        if (has_conflict)
        {
          continue;
        }

        // Skip reference combination if not requested (but allow mods-only peptides)
        if (var_combo == 0 && !include_reference && local_mods.empty())
        {
          continue;
        }

        // Parse variant peptide
        AASequence variant_peptide;
        try
        {
          variant_peptide = AASequence::fromString(variant_peptide_str);
        }
        catch (const Exception::BaseException&)
        {
          continue;  // Skip invalid sequences
        }

        // Step 4: Enumerate PEFF modification combinations
        auto peff_mod_peptides = enumeratePEFFModifications_(variant_peptide, local_mods, variant_desc);

        // Filter out the base (no-PEFF-mod) entry when reference not wanted
        // This handles the case where var_combo == 0 but local_mods is not empty:
        // enumeratePEFFModifications_ returns all combinations including combo=0 (no mods),
        // which would be the unmodified reference peptide that should be excluded.
        if (var_combo == 0 && !include_reference && !local_mods.empty())
        {
          String variant_peptide_string = variant_peptide.toString();
          peff_mod_peptides.erase(
            std::remove_if(peff_mod_peptides.begin(), peff_mod_peptides.end(),
              [&](const std::pair<String, AASequence>& p) {
                // Remove entries where no PEFF mods were applied (base peptide)
                return p.second.toString() == variant_peptide_string;
              }),
            peff_mod_peptides.end());
        }

        // Step 5: Apply fixed modifications (e.g., Carbamidomethyl on all Cys)
        // and then variable modifications (like Oxidation)
        for (auto& [peff_desc, peff_peptide] : peff_mod_peptides)
        {
          // Apply fixed modifications to all compatible residues
          if (!fixed_mods.empty())
          {
            ModifiedPeptideGenerator::applyFixedModifications(fixed_mods_map, peff_peptide);
          }

          if (variable_mods.empty())
          {
            // No variable mods - add directly (with fixed mods already applied)
            descriptions.push_back(peff_desc);
            sequences.push_back(peff_peptide);
          }
          else
          {
            // Apply variable modifications using ModifiedPeptideGenerator
            std::vector<AASequence> var_mod_peptides;
            ModifiedPeptideGenerator::applyVariableModifications(
              var_mods_map,
              peff_peptide,
              max_variable_mods_per_peptide,
              var_mod_peptides,
              true  // keep_original
            );

            for (auto& var_mod_pep : var_mod_peptides)
            {
              String full_desc = peff_desc;
              // If additional mods were applied, note it in description
              if (var_mod_pep.isModified() && var_mod_pep.toString() != peff_peptide.toString())
              {
                if (!full_desc.empty()) full_desc += ";";
                full_desc += "+varmod";
              }
              descriptions.push_back(full_desc);
              sequences.push_back(std::move(var_mod_pep));
            }
          }
        }
      }
    }
  }

  void PEFFFile::load(const String& filename,
                      std::vector<PEFFEntry>& entries,
                      std::vector<PEFFDatabaseMetadata>& headers) const
  {
    startProgress(0, 1, "Loading PEFF file");
    entries.clear();
    headers.clear();

    PEFFEntry entry;
    PEFFFile f;
    f.readStart(filename);
    headers = f.getHeaders();

    while (f.readNext(entry))
    {
      entries.push_back(std::move(entry));
    }
    endProgress();
  }

  void PEFFFile::store(const String& filename,
                       const std::vector<PEFFEntry>& entries,
                       const PEFFDatabaseMetadata& header) const
  {
    startProgress(0, entries.size(), "Writing PEFF file");
    PEFFFile f;
    f.writeStart(filename, header);
    for (const PEFFEntry& entry : entries)
    {
      f.writeNext(entry);
      nextProgress();
    }
    f.writeEnd();
    endProgress();
  }

  void PEFFFile::store(const String& filename,
                       const std::vector<PEFFEntry>& entries,
                       const std::vector<PEFFDatabaseMetadata>& headers) const
  {
    if (headers.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "PEFF: store() requires at least one database header");
    }

    startProgress(0, entries.size(), "Writing PEFF file");
    PEFFFile f;
    f.writeStart(filename, headers);
    for (const PEFFEntry& entry : entries)
    {
      f.writeNext(entry);
      nextProgress();
    }
    f.writeEnd();
    endProgress();
  }

  void PEFFFile::readStart(const String& filename)
  {
    if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    if (infile_.is_open())
    {
      infile_.close();
    }

    infile_.open(filename.c_str(), std::ios::binary | std::ios::in);
    if (!infile_.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    infile_.seekg(0, infile_.end);
    fileSize_ = infile_.tellg();
    infile_.seekg(0, infile_.beg);

    headers_.clear();
    entries_read_ = 0;

    // Parse header section
    // PEFF has two block types:
    // 1. File Description Block: # PEFF N.N, optional GeneralComment lines, ends with # //
    // 2. Database Description Block: starts with # DbName=, has mandatory keys, ends with # //
    // When we encounter # // for a block with no DbName, we treat it as a file description block
    // and carry its data (version, general_comments) into the next database block.
    std::streambuf* sb = infile_.rdbuf();
    std::string line;
    bool in_header = true;
    PEFFDatabaseMetadata current_header;
    bool has_header = false;
    // File-level data from the file description block (no DbName)
    String file_version;
    std::vector<String> file_comments;

    while (in_header && sb->sgetc() != std::streambuf::traits_type::eof())
    {
      int c = sb->sgetc();

      // Skip leading whitespace
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
      {
        sb->sbumpc();
        continue;
      }

      if (c == '#')
      {
        // Read the header line
        line.clear();
        while ((c = sb->sbumpc()) != std::streambuf::traits_type::eof() && c != '\n')
        {
          if (c != '\r')
          {
            line += static_cast<char>(c);
          }
        }

        // Parse the header line (line starts with # or is just #)
        bool is_db_separator = false;
        parseHeaderLine_(String(line), current_header, is_db_separator);

        if (is_db_separator)
        {
          // The # // line marks end of a block
          if (has_header)
          {
            if (current_header.db_name.empty() && headers_.empty())
            {
              // This is the file description block (no DbName) — save its data for the next block
              file_version = current_header.version;
              file_comments = current_header.general_comments;
              current_header = PEFFDatabaseMetadata();
              has_header = false;
            }
            else
            {
              // This is a database description block — merge file-level data if available
              if (!file_version.empty())
              {
                current_header.version = file_version;
              }
              if (!file_comments.empty() && current_header.general_comments.empty())
              {
                current_header.general_comments = file_comments;
              }
              headers_.push_back(current_header);
              current_header = PEFFDatabaseMetadata();
              has_header = false;
            }
          }
        }
        else
        {
          has_header = true;
        }
      }
      else
      {
        // Not a header line, we're done with headers
        in_header = false;
        // Don't push again if we already pushed on # //
        if (has_header)
        {
          // Merge file-level data
          if (!file_version.empty())
          {
            current_header.version = file_version;
          }
          if (!file_comments.empty() && current_header.general_comments.empty())
          {
            current_header.general_comments = file_comments;
          }
          headers_.push_back(current_header);
        }
      }
    }

    // If no headers found, create a default one
    if (headers_.empty())
    {
      headers_.emplace_back();
    }

    // Validate required header keys (spec section 3.3.1)
    for (Size i = 0; i < headers_.size(); ++i)
    {
      const auto& h = headers_[i];
      String db_label = h.db_name.empty() ? String("database block ") + String(i + 1) : h.db_name;

      // Fix #16: Database blocks must start with DbName
      if (h.db_name.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Database block should start with 'DbName' key in " << db_label << ".\n";
      }
      if (h.prefix.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'Prefix' missing in " << db_label << ".\n";
      }
      if (h.db_version.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'DbVersion' missing in " << db_label << ".\n";
      }
      if (h.db_sources.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'DbSource' missing in " << db_label << ".\n";
      }
      if (h.number_of_entries == 0)
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'NumberOfEntries' missing or zero in " << db_label << ".\n";
      }
      // SequenceType always has a default (AA), so no warning needed

      // Fix #6: ProteoformDb and HasAnnotationIdentifiers mutual exclusion (spec section 3.4.2)
      if (h.is_proteoform_db && h.has_annotation_identifiers)
      {
        OPENMS_LOG_WARN << "PEFF: ProteoformDb and HasAnnotationIdentifiers MUST NOT both be true in " << db_label << ".\n";
      }
    }
  }

  bool PEFFFile::readNext(PEFFEntry& entry)
  {
    if (infile_.eof() || atEnd())
    {
      return false;
    }

    seq_.clear();
    id_.clear();
    description_.clear();

    if (!readEntry_(id_, description_, seq_))
    {
      if (entries_read_ == 0)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "", "The first PEFF entry could not be read!");
      }
      else
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "", "Only " + String(entries_read_) + " entries could be read. Parsing next record failed.");
      }
    }
    ++entries_read_;

    // Initialize entry
    entry = PEFFEntry();
    entry.identifier = id_;
    entry.sequence = seq_;
    entry.sequence_length = seq_.size();

    // Extract prefix from identifier (format: Prefix:DbUniqueId)
    size_t colon_pos = id_.find(':');
    if (colon_pos != std::string::npos)
    {
      entry.prefix = id_.substr(0, colon_pos);
    }

    // Parse annotations from description
    parseAnnotations_(String(description_), entry);

    return true;
  }

  bool PEFFFile::readEntry_(std::string& id, std::string& description, std::string& seq)
  {
    std::streambuf* sb = infile_.rdbuf();
    bool keep_reading = true;
    bool description_exists = true;

    // Skip leading whitespace (including newlines, tabs, spaces)
    int c;
    while ((c = sb->sgetc()) != std::streambuf::traits_type::eof() &&
           (c == ' ' || c == '\t' || c == '\n' || c == '\r'))
    {
      sb->sbumpc();
    }

    // Fix #12: Lines beginning with ";" are not permitted in PEFF (spec section 3.3.3)
    if (sb->sgetc() == ';')
    {
      OPENMS_LOG_WARN << "PEFF: Lines beginning with ';' are not permitted.\n";
      int ch;
      while ((ch = sb->sbumpc()) != std::streambuf::traits_type::eof() && ch != '\n') {}
      return false;
    }

    if (sb->sbumpc() != '>')
    {
      return false;
    }

    while (keep_reading) // reading the ID
    {
      c = sb->sbumpc(); // get and advance to next char
      switch (c)
      {
        case ' ':
        case '\t':
          if (!id.empty())
          {
            keep_reading = false; // ID finished
          }
          break;
        case '\n': // ID finished and no description available
          keep_reading = false;
          description_exists = false;
          break;
        case '\r':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return false;
        default:
          id += static_cast<char>(c);
      }
    }

    if (id.empty())
    {
      return false;
    }

    if (description_exists)
    {
      keep_reading = true;
    }

    // reading the description
    while (keep_reading)
    {
      c = sb->sbumpc(); // get and advance to next char
      switch (c)
      {
        case '\n': // description finished
          keep_reading = false;
          break;
        case '\r':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return false;
        default:
          description += static_cast<char>(c);
      }
    }

    // reading the sequence
    keep_reading = true;
    while (keep_reading)
    {
      c = sb->sbumpc(); // get and advance to next char
      switch (c)
      {
        case '\n':
          if (sb->sgetc() == '>') // reaching the beginning of the next protein-entry
          {
            keep_reading = false;
          }
          break;
        case '\r': // not saving white spaces
        case ' ':
        case '\t':
          break;
        case std::streambuf::traits_type::eof():
          infile_.setstate(std::ios::eofbit);
          return !seq.empty();
        default:
          seq += static_cast<char>(c);
      }
    }
    return !seq.empty();
  }

  const std::vector<PEFFDatabaseMetadata>& PEFFFile::getHeaders() const
  {
    return headers_;
  }

  bool PEFFFile::atEnd() const
  {
    return const_cast<std::fstream&>(infile_).peek() == std::streambuf::traits_type::eof();
  }

  void PEFFFile::writeStart(const String& filename, const PEFFDatabaseMetadata& header)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::PEFF))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "invalid file extension; expected '" +
                                          FileTypes::typeToName(FileTypes::PEFF) + "'");
    }

    outfile_.open(filename.c_str(), ofstream::out | ofstream::binary);

    if (!outfile_.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // Write header
    outfile_ << formatHeader_(header);
  }

  void PEFFFile::writeStart(const String& filename, const std::vector<PEFFDatabaseMetadata>& headers)
  {
    if (headers.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "PEFF: writeStart() requires at least one database header");
    }

    if (!FileHandler::hasValidExtension(filename, FileTypes::PEFF))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                          "invalid file extension; expected '" +
                                          FileTypes::typeToName(FileTypes::PEFF) + "'");
    }

    outfile_.open(filename.c_str(), ofstream::out | ofstream::binary);

    if (!outfile_.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // Write header
    outfile_ << formatHeader_(headers);
  }

  void PEFFFile::writeNext(const PEFFEntry& entry)
  {
    outfile_ << formatEntry_(entry);
  }

  void PEFFFile::writeEnd()
  {
    outfile_.close();
  }

  bool PEFFFile::isPEFFFile(const String& filename)
  {
    if (!File::exists(filename) || !File::readable(filename))
    {
      return false;
    }

    std::ifstream file(filename.c_str());
    std::string line;

    // Skip whitespace and check for PEFF header
    while (std::getline(file, line))
    {
      // Trim whitespace
      size_t start = line.find_first_not_of(" \t\r\n");
      if (start == std::string::npos)
      {
        continue; // Empty line
      }

      // Check if it's a comment line
      if (line[start] == '#')
      {
        // Check for PEFF version marker - must be "# PEFF x.x" format at start
        // Skip # and whitespace
        size_t content_start = start + 1;
        while (content_start < line.size() && (line[content_start] == ' ' || line[content_start] == '\t'))
        {
          content_start++;
        }
        // Check if remaining content starts with "PEFF " followed by version number
        if (content_start + 5 < line.size() &&
            line.substr(content_start, 5) == "PEFF ")
        {
          // Check for version number pattern after "PEFF "
          size_t version_start = content_start + 5;
          if (std::isdigit(line[version_start]))
          {
            return true;
          }
        }
      }
      else if (line[start] == '>')
      {
        // Check if description contains PEFF annotations (backslash keywords)
        if (line.find('\\') != std::string::npos)
        {
          return true;
        }
        return false;
      }
      else
      {
        return false;
      }
    }
    return false;
  }

  String PEFFFile::toProForma(const PEFFEntry& entry)
  {
    // If no modifications, return plain sequence
    if (entry.modifications.empty())
    {
      return entry.sequence;
    }

    // Separate modifications by position
    // Position 0 = unlocalised (unknown position in PEFF)
    std::vector<const PEFFModification*> unlocalised_mods;
    std::map<Size, std::vector<const PEFFModification*>> mods_by_pos;

    for (const auto& mod : entry.modifications)
    {
      if (mod.position == 0)
      {
        unlocalised_mods.push_back(&mod);
      }
      else
      {
        mods_by_pos[mod.position].push_back(&mod);
      }
    }

    // Helper to format a single modification in ProForma notation
    auto formatMod = [](const PEFFModification& mod) -> String
    {
      // Prefer accession if available (UNIMOD: or MOD:)
      if (!mod.accession.empty())
      {
        return "[" + mod.accession + "]";
      }
      // Fall back to name
      if (!mod.name.empty())
      {
        return "[" + mod.name + "]";
      }
      return "";
    };

    String result;

    // Add unlocalised modifications as prefix: <?[mod]>
    // In ProForma, unlocalised mods are written as: <?[mod1][mod2]>SEQUENCE
    if (!unlocalised_mods.empty())
    {
      result += "<?";
      for (const auto* mod : unlocalised_mods)
      {
        result += formatMod(*mod);
      }
      result += ">";
    }

    // Build sequence with localised modifications
    // Modifications in ProForma appear right after the amino acid they modify
    for (Size i = 0; i < entry.sequence.size(); ++i)
    {
      // Add the amino acid (1-based position is i+1)
      result += entry.sequence[i];

      // Add any modifications at this position (1-based)
      Size pos = i + 1;
      auto it = mods_by_pos.find(pos);
      if (it != mods_by_pos.end())
      {
        for (const auto* mod : it->second)
        {
          result += formatMod(*mod);
        }
      }
    }

    return result;
  }

  void PEFFFile::parseHeaderLine_(const String& line, PEFFDatabaseMetadata& header, bool& new_db)
  {
    new_db = false;

    // Line includes the leading # character
    String trimmed(line);
    trimmed.trim();

    // Check for database separator (# // or #//, but not bare # which is just a blank comment)
    if (trimmed == "# //" || trimmed == "#//")
    {
      new_db = true;
      return;
    }

    // Skip blank comment lines (just #)
    if (trimmed == "#")
    {
      return;
    }

    // Skip the # and any leading whitespace
    if (!trimmed.hasPrefix("#"))
    {
      return;
    }
    String content = trimmed.substr(1);
    content.trim();

    // Check for PEFF version line
    if (content.hasPrefix("PEFF"))
    {
      // Extract version (e.g., "PEFF 1.0" or "PEFF1.0")
      String version_str = content.substr(4);
      version_str.trim();
      if (!version_str.empty())
      {
        header.version = version_str;
      }
      return;
    }

    // Parse Key=Value format
    size_t eq_pos = content.find('=');
    if (eq_pos == String::npos)
    {
      return;
    }

    String key = content.substr(0, eq_pos);
    key.trim();
    String value = content.substr(eq_pos + 1);
    value.trim();

    if (key == "DbName")
    {
      header.db_name = value;
    }
    else if (key == "Prefix")
    {
      header.prefix = value;
    }
    else if (key == "DbDescription")
    {
      header.db_description = value;
    }
    else if (key == "Decoy")
    {
      header.is_decoy = (value.toLower() == "true" || value == "1");
    }
    else if (key == "DbSource")
    {
      header.db_sources.push_back(value);
    }
    else if (key == "DbVersion")
    {
      header.db_version = value;
    }
    else if (key == "NumberOfEntries")
    {
      header.number_of_entries = value.toInt();
    }
    else if (key == "SequenceType")
    {
      header.sequence_type = (value == "NA") ? PEFFDatabaseMetadata::SequenceType::NA
                                              : PEFFDatabaseMetadata::SequenceType::AA;
    }
    else if (key == "GeneralComment")
    {
      if (value.empty())
      {
        OPENMS_LOG_WARN << "PEFF: GeneralComment MUST NOT be empty if present (spec section 3.3.1).\n";
      }
      else
      {
        header.general_comments.push_back(value);
      }
    }
    else if (key == "Conversion")
    {
      header.conversion = value;
    }
    else if (key == "SpecificKey")
    {
      // Format: SpecificKey=KeyName:description
      size_t colon_pos = value.find(':');
      if (colon_pos != String::npos)
      {
        String sk_key = value.substr(0, colon_pos);
        String sk_desc = value.substr(colon_pos + 1);
        header.specific_keys[sk_key] = sk_desc;
      }
    }
    else if (key == "SpecificValue")
    {
      // Format: SpecificValue=KeyName:(type_definition)
      size_t colon_pos = value.find(':');
      if (colon_pos != String::npos)
      {
        String sv_key = value.substr(0, colon_pos);
        String sv_type = value.substr(colon_pos + 1);
        header.specific_values[sv_key] = sv_type;
      }
    }
    else if (key == "OptionalTagDef")
    {
      header.optional_tag_defs.push_back(value);
    }
    else if (key == "CustomKeyDef")
    {
      header.custom_key_defs.push_back(parseCustomKeyDef_(value));
    }
    else if (key == "HasAnnotationIdentifiers")
    {
      header.has_annotation_identifiers = (value.toLower() == "true" || value == "1");
    }
    else if (key == "ProteoformDb" || key == "ProteoformDB" || key == "IsProteoformDB")
    {
      header.is_proteoform_db = (value.toLower() == "true" || value == "1");
    }
    else
    {
      // Store unrecognized keys for round-trip preservation (e.g., DbDate from UniProt)
      header.unrecognized_keys[key] = value;
    }
  }

  void PEFFFile::parseAnnotations_(const String& description, PEFFEntry& entry)
  {
    if (description.empty())
    {
      return;
    }

    // Split by backslash to find annotation blocks
    // We need to find unescaped backslashes (annotation separators)
    // In PEFF, \\ \| \( \) are escapes within a value, not separators
    size_t pos = 0;
    size_t prev = 0;

    // First part before any backslash is not an annotation
    // Find first unescaped backslash
    auto is_key_start = [](const char c) -> bool
    {
      return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
    };
    auto is_escape_sequence = [&description](const size_t i) -> bool
    {
      if (i + 1 >= description.size()) return false;
      const char next = description[i + 1];
      return next == '\\' || next == '|' || next == '(' || next == ')';
    };
    auto find_next_key_sep = [&description, &is_key_start, &is_escape_sequence](const size_t from) -> size_t
    {
      for (size_t i = from; i < description.size(); ++i)
      {
        if (description[i] != '\\' || i + 1 >= description.size())
        {
          continue;
        }
        if (is_escape_sequence(i))
        {
          ++i; // skip escaped char
          continue;
        }
        if (is_key_start(description[i + 1]))
        {
          return i;
        }
      }
      return String::npos;
    };

    pos = String::npos;
    pos = find_next_key_sep(0);
    if (pos == String::npos)
    {
      // No annotations, entire description could be protein name
      if (!description.empty())
      {
        String trimmed_desc(description);
        trimmed_desc.trim();
        entry.protein_names.push_back(trimmed_desc);
      }
      return;
    }

    // Handle text before first annotation (if any)
    String before_annotations = description.substr(0, pos);
    before_annotations.trim();
    if (!before_annotations.empty())
    {
      entry.protein_names.push_back(before_annotations);
    }

    // Determine corresponding database header (by Prefix) for CustomKeyDef validation.
    const PEFFDatabaseMetadata* header_for_entry = nullptr;
    if (!headers_.empty())
    {
      if (!entry.prefix.empty())
      {
        for (const auto& h : headers_)
        {
          if (h.prefix == entry.prefix)
          {
            header_for_entry = &h;
            break;
          }
        }
      }
      if (header_for_entry == nullptr && headers_.size() == 1)
      {
        header_for_entry = &headers_.front();
      }
    }

    const Size seq_len = entry.sequence.size();
    const bool is_nucleotide_db = (header_for_entry != nullptr &&
                                   header_for_entry->sequence_type == PEFFDatabaseMetadata::SequenceType::NA);

    auto is_valid_variant_code = [is_nucleotide_db](const char c) -> bool
    {
      if (c == '*') return true;
      if (!std::isalpha(static_cast<unsigned char>(c))) return false;
      if (!is_nucleotide_db) return true;

      const char up = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      // IUPAC nucleotide ambiguity codes.
      const std::string allowed = "ACGTURYSWKMBDHVN";
      return allowed.find(up) != std::string::npos;
    };

    // Parse each annotation
    std::set<String> seen_keys; // Fix #7: track seen keys for duplicate detection
    prev = pos;
    while (prev < description.size())
    {
      // Skip the backslash
      prev++;

      // Find the next annotation backslash (key separator) or end
      pos = String::npos;
      pos = find_next_key_sep(prev);
      if (pos == String::npos)
      {
        pos = description.size();
      }

      String annotation = description.substr(prev, pos - prev);
      annotation.trim();
      prev = pos;

      if (annotation.empty())
      {
        continue;
      }

      // Parse Key=Value or Key=(Value) format
      size_t eq_pos = annotation.find('=');
      if (eq_pos == String::npos)
      {
        continue;
      }

      String key = annotation.substr(0, eq_pos);
      String value = annotation.substr(eq_pos + 1);

      // Fix #7: Detect duplicate keys (spec section 3.3.3: same key MUST NOT appear twice)
      if (seen_keys.count(key))
      {
        OPENMS_LOG_WARN << "PEFF: Duplicate key '" << key << "' in entry '" << entry.identifier << "'.\n";
      }
      seen_keys.insert(key);

      if (key == "PName")
      {
        // Protein names can be multiple: (Name1)(Name2)
        std::vector<String> names = parseParenList_(value);
        if (names.empty() && !value.empty())
        {
          entry.protein_names.push_back(value);
        }
        else
        {
          for (const String& name : names)
          {
            entry.protein_names.push_back(name);
          }
        }
      }
      else if (key == "GName")
      {
        entry.gene_name = value;
      }
      else if (key == "NcbiTaxId")
      {
        entry.ncbi_tax_id = value.toInt();
      }
      else if (key == "TaxName")
      {
        entry.taxonomy_name = value;
      }
      else if (key == "Length")
      {
        entry.sequence_length = value.toInt();
      }
      else if (key == "SV")
      {
        entry.sequence_version = value;
      }
      else if (key == "EV")
      {
        entry.entry_version = value;
      }
      else if (key == "PE")
      {
        entry.protein_existence = value.toInt();
      }
      else if (key == "DbUniqueId")
      {
        entry.db_unique_id = value;
      }
      else if (key == "ID")
      {
        entry.entry_id = value;
      }
      else if (key == "OX")
      {
        // OX is a shorthand for NcbiTaxId
        entry.ncbi_tax_id = value.toInt();
      }
      else if (key == "AltAC")
      {
        // May be parenthesized list: (AC1)(AC2) or single value
        std::vector<String> acs = parseParenList_(value);
        for (const String& ac : acs)
        {
          entry.alt_accessions.push_back(ac);
        }
      }
      else if (key == "ModRes" || key == "ModResPsi" || key == "ModResUnimod")
      {
        if (is_nucleotide_db)
        {
          OPENMS_LOG_WARN << "PEFF: " << key << " annotations are not applicable to nucleotide databases (entry '"
                          << entry.identifier << "').\n";
        }

        std::vector<String> mods = parseParenList_(value);
        for (const String& mod_str : mods)
        {
          // Basic structural checks per spec sections 3.3.10-3.3.12.
          const std::vector<String> parts_check = splitByPipeEscapeAware(mod_str);
          if (parts_check.size() < 3)
          {
            OPENMS_LOG_WARN << "PEFF: " << key << " tuple must have at least 3 fields (position|accession|name) in entry '"
                            << entry.identifier << "'.\n";
          }
          else if (parts_check[2].empty())
          {
            OPENMS_LOG_WARN << "PEFF: " << key << " tuple is missing required name field in entry '"
                            << entry.identifier << "'.\n";
          }
          if ((key == "ModResPsi" || key == "ModResUnimod") && parts_check.size() >= 2 && parts_check[1].empty())
          {
            OPENMS_LOG_WARN << "PEFF: " << key << " tuple is missing required accession field in entry '"
                            << entry.identifier << "'.\n";
          }

          if (parts_check.size() >= 2 && !parts_check[1].empty())
          {
            if (key == "ModResPsi" && !parts_check[1].hasPrefix("MOD:"))
            {
              OPENMS_LOG_WARN << "PEFF: " << key << " accession should start with 'MOD:' in entry '"
                              << entry.identifier << "'.\n";
            }
            else if (key == "ModResUnimod" && !parts_check[1].hasPrefix("UNIMOD:"))
            {
              OPENMS_LOG_WARN << "PEFF: " << key << " accession should start with 'UNIMOD:' in entry '"
                              << entry.identifier << "'.\n";
            }
          }

          // Validate position(s) if present and known.
          if (!parts_check.empty())
          {
            auto [annot_id, pos_field] = extractAnnotationId(parts_check[0]);
            (void)annot_id;

            std::vector<String> pos_tokens;
            if (pos_field.find(',') != String::npos)
            {
              pos_field.split(',', pos_tokens);
            }
            else
            {
              pos_tokens.push_back(pos_field);
            }

            for (const String& tok : pos_tokens)
            {
              if (tok == "?" || tok.empty())
              {
                continue; // unknown position allowed
              }
              if (!isUnsignedInteger_(tok) || tok == "0")
              {
                OPENMS_LOG_WARN << "PEFF: " << key << " position must be a positive integer or '?' in entry '"
                                << entry.identifier << "'.\n";
                continue;
              }
              const Size p = static_cast<Size>(tok.toInt());
              if (seq_len > 0 && p > seq_len)
              {
                OPENMS_LOG_WARN << "PEFF: " << key << " position " << p << " exceeds sequence length "
                                << seq_len << " in entry '" << entry.identifier << "'.\n";
              }
            }
          }

          PEFFModification mod = parseModification_(mod_str);
          // Override type based on key name if accession didn't determine it
          if (key == "ModResPsi")
          {
            mod.type = PEFFModification::Type::PSI_MOD;
          }
          else if (key == "ModResUnimod")
          {
            mod.type = PEFFModification::Type::UNIMOD;
          }

          // Warn if name is missing for PSI_MOD or UNIMOD types (spec requires both accession and name)
          if ((mod.type == PEFFModification::Type::PSI_MOD || mod.type == PEFFModification::Type::UNIMOD) && mod.name.empty())
          {
            OPENMS_LOG_WARN << "PEFF: Modification with accession '" << mod.accession << "' is missing required name field." << "\n";
          }

	          // Expand comma-separated positions (spec sections 3.3.10-3.3.12)
	          // Check the original tuple for comma-separated positions
	          if (!parts_check.empty())
	          {
	            auto [annot_id, pos_field] = extractAnnotationId(parts_check[0]);
	            if (pos_field.find(',') != String::npos)
            {
              std::vector<String> positions;
              pos_field.split(',', positions);
              for (const String& pos_str : positions)
              {
                if (!(pos_str == "?" || pos_str.empty()) && (!isUnsignedInteger_(pos_str) || pos_str == "0"))
                {
                  OPENMS_LOG_WARN << "PEFF: " << key << " position must be a positive integer or '?' in entry '"
                                  << entry.identifier << "'.\n";
                }
                PEFFModification expanded_mod = mod;
                expanded_mod.annotation_id = annot_id;
                if (pos_str == "?" || pos_str.empty())
                {
                  expanded_mod.position = 0;
                }
                else
                {
                  expanded_mod.position = pos_str.toInt();
                }
                if (expanded_mod.position > 0 && seq_len > 0 && expanded_mod.position > seq_len)
                {
                  OPENMS_LOG_WARN << "PEFF: " << key << " position " << expanded_mod.position
                                  << " exceeds sequence length " << seq_len
                                  << " in entry '" << entry.identifier << "'.\n";
                }
                entry.modifications.push_back(expanded_mod);
              }
              continue; // Skip the default push_back below
            }
          }

          if (mod.position > 0 && seq_len > 0 && mod.position > seq_len)
          {
            OPENMS_LOG_WARN << "PEFF: " << key << " position " << mod.position << " exceeds sequence length "
                            << seq_len << " in entry '" << entry.identifier << "'.\n";
          }
          entry.modifications.push_back(mod);
        }
      }
      else if (key == "VariantSimple")
      {
        std::vector<String> variants = parseParenList_(value);
        for (const String& var : variants)
        {
          const std::vector<String> parts_check = splitByPipeEscapeAware(var);
          if (parts_check.size() >= 2 && parts_check[1].size() != 1)
          {
            OPENMS_LOG_WARN << "PEFF: VariantSimple newAminoAcid must be exactly 1 character in entry '"
                            << entry.identifier << "'.\n";
          }
          if (!parts_check.empty())
          {
            auto [annot_id, pos_field] = extractAnnotationId(parts_check[0]);
            (void)annot_id;
            if (!isUnsignedInteger_(pos_field) || pos_field == "0")
            {
              OPENMS_LOG_WARN << "PEFF: VariantSimple position must be a positive integer in entry '"
                              << entry.identifier << "'.\n";
            }
          }

          PEFFVariantSimple parsed = parseVariantSimple_(var);
          // Fix #9: Validate VariantSimple (spec section 3.3.8)
          if (parsed.position == 0)
          {
            OPENMS_LOG_WARN << "PEFF: VariantSimple position must be > 0 in entry '" << entry.identifier << "'.\n";
          }
          if (seq_len > 0 && parsed.position > seq_len)
          {
            OPENMS_LOG_WARN << "PEFF: VariantSimple position " << parsed.position << " exceeds sequence length "
                            << seq_len << " in entry '" << entry.identifier << "'.\n";
          }
          if (parsed.variant_aa == '\0')
          {
            OPENMS_LOG_WARN << "PEFF: VariantSimple variant_aa must not be empty in entry '" << entry.identifier << "'.\n";
          }
          else if (!is_valid_variant_code(parsed.variant_aa))
          {
            OPENMS_LOG_WARN << "PEFF: VariantSimple newAminoAcid '" << parsed.variant_aa
                            << "' is not valid in entry '" << entry.identifier << "'.\n";
          }
          entry.simple_variants.push_back(parsed);
        }
      }
      else if (key == "VariantComplex")
      {
        std::vector<String> variants = parseParenList_(value);
        for (const String& var : variants)
        {
          const std::vector<String> parts_check = splitByPipeEscapeAware(var);
          if (parts_check.size() < 3)
          {
            OPENMS_LOG_WARN << "PEFF: VariantComplex tuple must have at least 3 fields (start|end|newSequence) in entry '"
                            << entry.identifier << "'.\n";
          }
          if (parts_check.size() >= 2)
          {
            auto [annot_id, start_field] = extractAnnotationId(parts_check[0]);
            (void)annot_id;
            if (!isUnsignedInteger_(start_field) || start_field == "0" || !isUnsignedInteger_(parts_check[1]) || parts_check[1] == "0")
            {
              OPENMS_LOG_WARN << "PEFF: VariantComplex positions must be positive integers in entry '"
                              << entry.identifier << "'.\n";
            }
          }

          PEFFVariantComplex parsed = parseVariantComplex_(var);
          // Fix #10: Validate VariantComplex (spec section 3.3.9)
          // Single AA substitution is illegal as VariantComplex — must use VariantSimple
          if (parsed.start_position == parsed.end_position && parsed.replacement.size() == 1)
          {
            OPENMS_LOG_WARN << "PEFF: VariantComplex at position " << parsed.start_position
                            << " with single AA replacement must be encoded as VariantSimple in entry '"
                            << entry.identifier << "'.\n";
          }
          if (parsed.start_position == 0 || parsed.end_position == 0)
          {
            OPENMS_LOG_WARN << "PEFF: VariantComplex positions must be > 0 in entry '" << entry.identifier << "'.\n";
          }
          if (parsed.start_position > parsed.end_position)
          {
            OPENMS_LOG_WARN << "PEFF: VariantComplex startPosition must be <= endPosition in entry '"
                            << entry.identifier << "'.\n";
          }
          if (seq_len > 0 && (parsed.start_position > seq_len || parsed.end_position > seq_len))
          {
            OPENMS_LOG_WARN << "PEFF: VariantComplex positions exceed sequence length " << seq_len
                            << " in entry '" << entry.identifier << "'.\n";
          }
          if (!parsed.replacement.empty())
          {
            for (Size i = 0; i < parsed.replacement.size(); ++i)
            {
              const char c = parsed.replacement[i];
              if (!is_valid_variant_code(c))
              {
                OPENMS_LOG_WARN << "PEFF: VariantComplex newSequence contains invalid character '" << c
                                << "' in entry '" << entry.identifier << "'.\n";
                break;
              }
            }
          }
          entry.complex_variants.push_back(parsed);
        }
      }
      else if (key == "Variant")
      {
        // Deprecated in favor of VariantSimple and VariantComplex (spec section 3.3.7).
        OPENMS_LOG_WARN << "PEFF: Deprecated key 'Variant' encountered in entry '" << entry.identifier
                        << "'. Readers should use VariantSimple/VariantComplex.\n";
        entry.custom_annotations[key] = value;
      }
      else if (key == "Processed")
      {
        std::vector<String> regions = parseParenList_(value);
        for (const String& reg : regions)
        {
          const std::vector<String> parts_check = splitByPipeEscapeAware(reg);
          if (parts_check.size() >= 2)
          {
            auto [annot_id, start_field] = extractAnnotationId(parts_check[0]);
            (void)annot_id;
            if (!isUnsignedInteger_(start_field) || start_field == "0" || !isUnsignedInteger_(parts_check[1]) || parts_check[1] == "0")
            {
              OPENMS_LOG_WARN << "PEFF: Processed positions must be positive integers in entry '"
                              << entry.identifier << "'.\n";
            }
          }
          if (parts_check.size() < 4)
          {
            OPENMS_LOG_WARN << "PEFF: Processed tuple must have at least 4 fields (start|end|accession|name) in entry '"
                            << entry.identifier << "'.\n";
          }
          else
          {
            if (parts_check[2].empty())
            {
              OPENMS_LOG_WARN << "PEFF: Processed tuple is missing required accession field in entry '"
                              << entry.identifier << "'.\n";
            }
            if (parts_check[3].empty())
            {
              OPENMS_LOG_WARN << "PEFF: Processed tuple is missing required name field in entry '"
                              << entry.identifier << "'.\n";
            }
          }

          PEFFProcessedRegion parsed = parseProcessedRegion_(reg);
          if (parsed.start_position == 0 || parsed.end_position == 0)
          {
            OPENMS_LOG_WARN << "PEFF: Processed positions must be > 0 in entry '" << entry.identifier << "'.\n";
          }
          if (parsed.start_position > parsed.end_position)
          {
            OPENMS_LOG_WARN << "PEFF: Processed startPosition must be <= endPosition in entry '"
                            << entry.identifier << "'.\n";
          }
          if (seq_len > 0 && (parsed.start_position > seq_len || parsed.end_position > seq_len))
          {
            OPENMS_LOG_WARN << "PEFF: Processed positions exceed sequence length " << seq_len
                            << " in entry '" << entry.identifier << "'.\n";
          }
          if (!parsed.accession.empty() && !parsed.accession.hasPrefix("PEFF:"))
          {
            OPENMS_LOG_WARN << "PEFF: Processed accession should start with 'PEFF:' in entry '"
                            << entry.identifier << "'.\n";
          }
          entry.processed_regions.push_back(parsed);
        }
      }
      else if (key == "DisulfideBond")
      {
        std::vector<String> bonds = parseParenList_(value);
        for (const String& bond : bonds)
        {
          entry.disulfide_bonds.push_back(parseDisulfideBond_(bond));
        }
      }
      else if (key == "Proteoform")
      {
        std::vector<String> pforms = parseParenList_(value);
        for (const String& pf : pforms)
        {
          entry.proteoforms.push_back(pf);
        }
      }
      else
      {
        // Custom annotation
        if (header_for_entry != nullptr)
        {
          // Validate value against CustomKeyDef if present.
          const PEFFCustomKeyDef* ckd = nullptr;
          for (const auto& def : header_for_entry->custom_key_defs)
          {
            if (def.key_name == key)
            {
              ckd = &def;
              break;
            }
          }

          if (ckd != nullptr)
          {
            const Size expected_fields = !ckd->field_names.empty() ? ckd->field_names.size() : ckd->field_types.size();
            if (!ckd->field_names.empty() && !ckd->field_types.empty() && ckd->field_names.size() != ckd->field_types.size())
            {
              OPENMS_LOG_WARN << "PEFF: CustomKeyDef for '" << key
                              << "' has mismatched FieldNames/FieldTypes lengths.\n";
            }

            bool last_is_optional_tag = false;
            if (!ckd->field_names.empty() && !ckd->field_names.back().empty())
            {
              last_is_optional_tag = (ckd->field_names.back() == "OptionalTag");
            }

            bool has_regex = !ckd->regexp.empty();
            bool regex_valid = false;
            std::regex re;
            if (has_regex)
            {
              try
              {
                re = std::regex(ckd->regexp.c_str());
                regex_valid = true;
              }
              catch (const std::regex_error&)
              {
                OPENMS_LOG_WARN << "PEFF: CustomKeyDef RegExp for '" << key << "' is not a valid C++ regex.\n";
              }
            }

            const std::vector<String> items = parseParenList_(value);
            for (const String& item : items)
            {
              if (has_regex && regex_valid && !std::regex_match(item.c_str(), re))
              {
                OPENMS_LOG_WARN << "PEFF: Value of custom key '" << key << "' does not match RegExp in entry '"
                                << entry.identifier << "'.\n";
              }

              const std::vector<String> components = splitByPipeEscapeAware(item);
              if (expected_fields > 0)
              {
                const bool ok_count = (components.size() == expected_fields) ||
                                      (last_is_optional_tag && components.size() + 1 == expected_fields);
                if (!ok_count)
                {
                  OPENMS_LOG_WARN << "PEFF: Custom key '" << key << "' has " << components.size()
                                  << " component(s), expected " << expected_fields
                                  << " in entry '" << entry.identifier << "'.\n";
                }
              }

              const Size check_fields = std::min(components.size(), ckd->field_types.size());
              for (Size i = 0; i < check_fields; ++i)
              {
                if (!validateCustomKeyFieldType_(ckd->field_types[i], components[i]))
                {
                  OPENMS_LOG_WARN << "PEFF: Custom key '" << key << "' field " << (i + 1)
                                  << " does not match type '" << ckd->field_types[i]
                                  << "' in entry '" << entry.identifier << "'.\n";
                }
              }
            }
          }
        }
        entry.custom_annotations[key] = value;
      }
    }
  }

  std::vector<String> PEFFFile::parseParenList_(const String& value)
  {
    std::vector<String> result;

    if (value.empty())
    {
      return result;
    }

    // Handle format: (value1)(value2)(value3) or single value
    // Respects backslash escapes: \( \) \\ are not structural characters
    size_t pos = 0;
    while (pos < value.size())
    {
      if (value[pos] == '(')
      {
        // Find matching closing parenthesis, handling escapes and nested paired parens
        int depth = 1;
        size_t start = pos + 1;
        size_t end = start;

        while (end < value.size() && depth > 0)
        {
          if (value[end] == '\\' && end + 1 < value.size())
          {
            // Skip escaped character
            end += 2;
            continue;
          }
          if (value[end] == '(')
          {
            depth++;
          }
          else if (value[end] == ')')
          {
            depth--;
          }
          if (depth > 0)
          {
            end++;
          }
        }

        if (end > start)
        {
          result.push_back(value.substr(start, end - start));
        }
        pos = end + 1;
      }
      else
      {
        pos++;
      }
    }

    // If no parentheses found, treat entire value as single item
    if (result.empty() && !value.empty())
    {
      result.push_back(value);
    }

    return result;
  }

  PEFFModification PEFFFile::parseModification_(const String& tuple)
  {
    // Format: pos|accession|name|OptionalTag  (OptionalTag optional)
    // With annotation IDs: annotationID:pos|accession|name|OptionalTag
    // Position can be ? for unknown position
    PEFFModification mod;

    std::vector<String> parts = splitByPipeEscapeAware(tuple);

    if (parts.size() >= 2)
    {
      auto [annotation_id, pos_field] = extractAnnotationId(parts[0]);
      mod.annotation_id = annotation_id;

      // Position - may be comma-separated (e.g., "100,157,214") per spec sections 3.3.10-3.3.12
      // We store the first position here; additional positions are handled by the caller
      if (pos_field == "?" || pos_field.empty())
      {
        mod.position = 0;
      }
      else if (pos_field.find(',') != String::npos)
      {
        // Take first position; caller will handle expansion
        std::vector<String> positions;
        pos_field.split(',', positions);
        if (!positions.empty() && !positions[0].empty())
        {
          mod.position = positions[0].toInt();
        }
      }
      else
      {
        mod.position = pos_field.toInt();
      }

      // Accession
      mod.accession = parts[1];

      // Determine type from accession
      if (mod.accession.hasPrefix("MOD:"))
      {
        mod.type = PEFFModification::Type::PSI_MOD;
      }
      else if (mod.accession.hasPrefix("UNIMOD:"))
      {
        mod.type = PEFFModification::Type::UNIMOD;
      }
      else
      {
        mod.type = PEFFModification::Type::GENERIC;
      }

      // Name
      if (parts.size() >= 3)
      {
        mod.name = parts[2];
      }

      // OptionalTag
      if (parts.size() >= 4)
      {
        mod.optional_tag = parts[3];
      }
      if (parts.size() > 4)
      {
        OPENMS_LOG_WARN << "PEFF: Modification tuple has more than 4 fields; treating last field as OptionalTag.\n";
        mod.optional_tag = parts.back();
      }
    }

    return mod;
  }

  PEFFVariantSimple PEFFFile::parseVariantSimple_(const String& tuple)
  {
    // Format: pos|aa|OptionalTag (OptionalTag optional)
    // With annotation IDs: annotationID:pos|aa|OptionalTag
    PEFFVariantSimple var;

    std::vector<String> parts = splitByPipeEscapeAware(tuple);

    if (parts.size() >= 2)
    {
      auto [annotation_id, pos_field] = extractAnnotationId(parts[0]);
      var.annotation_id = annotation_id;
      var.position = pos_field.toInt();
      var.variant_aa = parts[1].empty() ? '\0' : parts[1][0];

      if (parts.size() >= 3)
      {
        var.optional_tag = parts[2];
      }
    }

    return var;
  }

  PEFFVariantComplex PEFFFile::parseVariantComplex_(const String& tuple)
  {
    // Format: start|end|replacement|OptionalTag (OptionalTag optional)
    // With annotation IDs: annotationID:start|end|replacement|OptionalTag
    PEFFVariantComplex var;

    std::vector<String> parts = splitByPipeEscapeAware(tuple);

    if (parts.size() >= 3)
    {
      auto [annotation_id, start_field] = extractAnnotationId(parts[0]);
      var.annotation_id = annotation_id;
      var.start_position = start_field.toInt();
      var.end_position = parts[1].toInt();
      var.replacement = parts[2];

      if (parts.size() >= 4)
      {
        var.optional_tag = parts[3];
      }
    }

    return var;
  }

  PEFFProcessedRegion PEFFFile::parseProcessedRegion_(const String& tuple)
  {
    // Format: start|end|accession|name|OptionalTag (name and OptionalTag optional)
    // With annotation IDs: annotationID:start|end|accession|name|OptionalTag
    PEFFProcessedRegion reg;

    std::vector<String> parts = splitByPipeEscapeAware(tuple);

    if (parts.size() >= 3)
    {
      auto [annotation_id, start_field] = extractAnnotationId(parts[0]);
      reg.annotation_id = annotation_id;
      reg.start_position = start_field.toInt();
      reg.end_position = parts[1].toInt();
      reg.accession = parts[2];

      if (parts.size() >= 4)
      {
        reg.name = parts[3];
      }

      if (parts.size() >= 5)
      {
        reg.optional_tag = parts[4];
      }
    }

    return reg;
  }

  PEFFDisulfideBond PEFFFile::parseDisulfideBond_(const String& tuple)
  {
    // Format: id1,id2|description  (description is optional)
    // With annotation IDs: annotationID:id1,id2|description
    // IDs can be annotation identifiers or positions
    PEFFDisulfideBond bond;

    // Split by | first to separate IDs from description
    size_t pipe_pos = tuple.find('|');
    String ids_part = (pipe_pos != String::npos) ? tuple.substr(0, pipe_pos) : tuple;

    if (pipe_pos != String::npos)
    {
      bond.optional_tag = tuple.substr(pipe_pos + 1);
    }

    // Extract annotation ID prefix from the IDs portion
    // The annotation ID prefix appears before the first ID: "annotationID:id1,id2"
    auto [annotation_id, remaining] = extractAnnotationId(ids_part);
    bond.annotation_id = annotation_id;

    // Split the IDs by comma
    size_t comma_pos = remaining.find(',');
    if (comma_pos != String::npos)
    {
      bond.id1 = remaining.substr(0, comma_pos);
      bond.id2 = remaining.substr(comma_pos + 1);
    }

    return bond;
  }

  namespace
  {
    void appendDatabaseDescriptionBlock_(std::ostringstream& out, const PEFFDatabaseMetadata& header)
    {
      // Validate mandatory fields on write (spec section 3.3.1)
      if (header.db_name.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'DbName' is empty on write.\n";
      }
      if (header.prefix.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'Prefix' is empty on write.\n";
      }
      if (header.db_version.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'DbVersion' is empty on write.\n";
      }
      if (header.db_sources.empty())
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'DbSource' is missing on write.\n";
      }
      if (header.number_of_entries == 0)
      {
        OPENMS_LOG_WARN << "PEFF: Required key 'NumberOfEntries' is zero on write.\n";
      }

      // Database Description Block
      if (!header.db_name.empty())
      {
        out << "# DbName=" << header.db_name << "\n";
      }

      if (!header.prefix.empty())
      {
        out << "# Prefix=" << header.prefix << "\n";
      }

      if (!header.db_description.empty())
      {
        out << "# DbDescription=" << header.db_description << "\n";
      }

      out << "# Decoy=" << (header.is_decoy ? "true" : "false") << "\n";

      for (const String& source : header.db_sources)
      {
        out << "# DbSource=" << source << "\n";
      }

      if (!header.db_version.empty())
      {
        out << "# DbVersion=" << header.db_version << "\n";
      }

      if (header.number_of_entries > 0)
      {
        out << "# NumberOfEntries=" << header.number_of_entries << "\n";
      }

      out << "# SequenceType=" << (header.sequence_type == PEFFDatabaseMetadata::SequenceType::AA ? "AA" : "NA") << "\n";

      if (!header.conversion.empty())
      {
        out << "# Conversion=" << header.conversion << "\n";
      }

      for (const auto& [key, desc] : header.specific_keys)
      {
        out << "# SpecificKey=" << key << ":" << desc << "\n";
      }

      for (const auto& [key, type_def] : header.specific_values)
      {
        out << "# SpecificValue=" << key << ":" << type_def << "\n";
      }

      for (const String& tag_def : header.optional_tag_defs)
      {
        out << "# OptionalTagDef=" << tag_def << "\n";
      }

      for (const auto& ckd : header.custom_key_defs)
      {
        out << "# CustomKeyDef=(";
        out << "KeyName=" << ckd.key_name;
        if (!ckd.description.empty())
        {
          out << "|Description=\"" << escapeCustomKeyDefQuotedValue_(ckd.description) << "\"";
        }
        if (!ckd.concept_curie.empty())
        {
          out << "|ConceptCURIE=" << ckd.concept_curie;
        }
        if (!ckd.regexp.empty())
        {
          out << "|RegExp=\"" << escapeCustomKeyDefQuotedValue_(ckd.regexp) << "\"";
        }
        if (!ckd.field_names.empty())
        {
          out << "|FieldNames=";
          for (Size i = 0; i < ckd.field_names.size(); ++i)
          {
            if (i > 0)
            {
              out << ",";
            }
            out << ckd.field_names[i];
          }
        }
        if (!ckd.field_types.empty())
        {
          out << "|FieldTypes=";
          for (Size i = 0; i < ckd.field_types.size(); ++i)
          {
            if (i > 0)
            {
              out << ",";
            }
            out << ckd.field_types[i];
          }
        }
        out << ")\n";
      }

      if (header.is_proteoform_db)
      {
        out << "# ProteoformDb=true\n";
      }

      if (header.has_annotation_identifiers)
      {
        out << "# HasAnnotationIdentifiers=true\n";
      }

      // Write unrecognized keys (round-trip preservation)
      for (const auto& [key, value] : header.unrecognized_keys)
      {
        out << "# " << key << "=" << value << "\n";
      }

      out << "# //\n";
    }
  } // anonymous namespace

  String PEFFFile::formatHeader_(const PEFFDatabaseMetadata& header) const
  {
    std::ostringstream out;

    // File Description Block
    String version = header.version;
    if (version.empty())
    {
      OPENMS_LOG_WARN << "PEFF: Missing PEFF version on write; defaulting to 1.0.\n";
      version = "1.0";
    }
    out << "# PEFF " << version << "\n";

    for (const String& comment : header.general_comments)
    {
      if (!comment.empty())
      {
        out << "# GeneralComment=" << comment << "\n";
      }
    }

    out << "# //\n";

    appendDatabaseDescriptionBlock_(out, header);

    return out.str();
  }

  String PEFFFile::formatHeader_(const std::vector<PEFFDatabaseMetadata>& headers) const
  {
    std::ostringstream out;
    if (headers.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "PEFF: formatHeader() requires at least one database header");
    }

    String version;
    for (const auto& h : headers)
    {
      if (!h.version.empty())
      {
        version = h.version;
        break;
      }
    }
    if (version.empty())
    {
      OPENMS_LOG_WARN << "PEFF: Missing PEFF version on write; defaulting to 1.0.\n";
      version = "1.0";
    }

    for (const auto& h : headers)
    {
      if (!h.version.empty() && h.version != version)
      {
        OPENMS_LOG_WARN << "PEFF: Multiple PEFF version values provided; writing '" << version << "'.\n";
        break;
      }
    }

    // Merge GeneralComment lines (file description block only).
    std::set<String> seen_comments;
    std::vector<String> comments;
    for (const auto& h : headers)
    {
      for (const String& c : h.general_comments)
      {
        if (c.empty()) continue;
        if (seen_comments.insert(c).second)
        {
          comments.push_back(c);
        }
      }
    }

    out << "# PEFF " << version << "\n";
    for (const String& comment : comments)
    {
      out << "# GeneralComment=" << comment << "\n";
    }
    out << "# //\n";

    for (const auto& h : headers)
    {
      appendDatabaseDescriptionBlock_(out, h);
    }

    return out.str();
  }

  String PEFFFile::formatEntry_(const PEFFEntry& entry) const
  {
    std::ostringstream out;

    // Start with >Prefix:DbUniqueId (or just >identifier if no prefix)
    out << ">" << entry.identifier;

    // Build description with annotations
    std::ostringstream desc;

    // Protein names
    if (!entry.protein_names.empty())
    {
      desc << " \\PName=";
      if (entry.protein_names.size() == 1)
      {
        desc << "(" << escape_peff(entry.protein_names[0]) << ")";
      }
      else
      {
        for (const String& name : entry.protein_names)
        {
          desc << "(" << escape_peff(name) << ")";
        }
      }
    }

    // Gene name
    if (!entry.gene_name.empty())
    {
      desc << " \\GName=" << escape_peff(entry.gene_name);
    }

    // NCBI Tax ID
    if (entry.ncbi_tax_id != 0)
    {
      desc << " \\NcbiTaxId=" << entry.ncbi_tax_id;
    }

    // Taxonomy name
    if (!entry.taxonomy_name.empty())
    {
      desc << " \\TaxName=" << escape_peff(entry.taxonomy_name);
    }

    // Length
    if (entry.sequence_length > 0)
    {
      desc << " \\Length=" << entry.sequence_length;
    }

    // Sequence version
    if (!entry.sequence_version.empty())
    {
      desc << " \\SV=" << escape_peff(entry.sequence_version);
    }

    // Entry version
    if (!entry.entry_version.empty())
    {
      desc << " \\EV=" << escape_peff(entry.entry_version);
    }

    // Protein existence
    if (entry.protein_existence > 0)
    {
      desc << " \\PE=" << entry.protein_existence;
    }

    // Database unique ID
    if (!entry.db_unique_id.empty())
    {
      desc << " \\DbUniqueId=" << escape_peff(entry.db_unique_id);
    }

    // Entry ID
    if (!entry.entry_id.empty())
    {
      desc << " \\ID=" << escape_peff(entry.entry_id);
    }

    // Alternative accessions - single key with parenthesized list
    if (!entry.alt_accessions.empty())
    {
      desc << " \\AltAC=";
      for (const String& alt_ac : entry.alt_accessions)
      {
        desc << "(" << escape_peff(alt_ac) << ")";
      }
    }

    // Modifications - group by type into separate keys
    if (!entry.modifications.empty())
    {
      // Helper lambda to write a list of modifications
      auto writeMods = [&desc](const std::vector<const PEFFModification*>& mods, const String& key)
      {
        if (mods.empty()) return;
        desc << " \\" << key << "=";
        for (const PEFFModification* mod : mods)
        {
          desc << "(";
          if (mod->annotation_id != std::numeric_limits<UInt>::max())
          {
            desc << mod->annotation_id << ":";
          }
          if (mod->position == 0)
          {
            desc << "?";
          }
          else
          {
            desc << mod->position;
          }
          desc << "|" << escape_peff(mod->accession);
          if (!mod->name.empty() || !mod->optional_tag.empty())
          {
            desc << "|" << escape_peff(mod->name);
          }
          if (!mod->optional_tag.empty())
          {
            desc << "|" << escape_peff(mod->optional_tag);
          }
          desc << ")";
        }
      };

      std::vector<const PEFFModification*> psi_mods, unimod_mods, generic_mods;
      for (const PEFFModification& mod : entry.modifications)
      {
        switch (mod.type)
        {
          case PEFFModification::Type::PSI_MOD:  psi_mods.push_back(&mod); break;
          case PEFFModification::Type::UNIMOD:   unimod_mods.push_back(&mod); break;
          case PEFFModification::Type::GENERIC:   generic_mods.push_back(&mod); break;
        }
      }

      writeMods(psi_mods, "ModResPsi");
      writeMods(unimod_mods, "ModResUnimod");
      writeMods(generic_mods, "ModRes");
    }

    // Simple variants
    if (!entry.simple_variants.empty())
    {
      desc << " \\VariantSimple=";
      for (const PEFFVariantSimple& var : entry.simple_variants)
      {
        desc << "(";
        if (var.annotation_id != std::numeric_limits<UInt>::max())
        {
          desc << var.annotation_id << ":";
        }
        desc << var.position << "|" << var.variant_aa;
        if (!var.optional_tag.empty())
        {
          desc << "|" << escape_peff(var.optional_tag);
        }
        desc << ")";
      }
    }

    // Complex variants
    if (!entry.complex_variants.empty())
    {
      desc << " \\VariantComplex=";
      for (const PEFFVariantComplex& var : entry.complex_variants)
      {
        desc << "(";
        if (var.annotation_id != std::numeric_limits<UInt>::max())
        {
          desc << var.annotation_id << ":";
        }
        desc << var.start_position << "|" << var.end_position << "|" << escape_peff(var.replacement);
        if (!var.optional_tag.empty())
        {
          desc << "|" << escape_peff(var.optional_tag);
        }
        desc << ")";
      }
    }

    // Processed regions
    if (!entry.processed_regions.empty())
    {
      desc << " \\Processed=";
      for (const PEFFProcessedRegion& reg : entry.processed_regions)
      {
        desc << "(";
        if (reg.annotation_id != std::numeric_limits<UInt>::max())
        {
          desc << reg.annotation_id << ":";
        }
        desc << reg.start_position << "|" << reg.end_position << "|" << escape_peff(reg.accession);
        // Always emit name separator if name or optional_tag is present (preserve field positions)
        if (!reg.name.empty() || !reg.optional_tag.empty())
        {
          desc << "|" << escape_peff(reg.name);  // May be empty placeholder
        }
        if (!reg.optional_tag.empty())
        {
          desc << "|" << escape_peff(reg.optional_tag);
        }
        desc << ")";
      }
    }

    // Disulfide bonds
    if (!entry.disulfide_bonds.empty())
    {
      desc << " \\DisulfideBond=";
      for (const PEFFDisulfideBond& bond : entry.disulfide_bonds)
      {
        desc << "(";
        if (bond.annotation_id != std::numeric_limits<UInt>::max())
        {
          desc << bond.annotation_id << ":";
        }
        desc << escape_peff(bond.id1) << "," << escape_peff(bond.id2);
        if (!bond.optional_tag.empty())
        {
          desc << "|" << escape_peff(bond.optional_tag);
        }
        desc << ")";
      }
    }

    // Proteoforms
    if (!entry.proteoforms.empty())
    {
      desc << " \\Proteoform=";
      for (const String& pf : entry.proteoforms)
      {
        desc << "(" << escape_peff(pf) << ")";
      }
    }

    // Custom annotations
    for (const auto& [key, value] : entry.custom_annotations)
    {
      // Values are stored in PEFF-encoded form (including escapes and item parentheses),
      // so we write them verbatim to preserve semantics and round-trip.
      desc << " \\" << key << "=" << value;
    }

    out << desc.str() << "\n";

    // Write sequence (80 characters per line)
    const String& seq = entry.sequence;
    Size chunks = seq.size() / 80;
    Size chunk_pos = 0;

    for (Size i = 0; i < chunks; ++i)
    {
      out.write(&seq[chunk_pos], 80);
      out << "\n";
      chunk_pos += 80;
    }

    if (seq.size() > chunk_pos)
    {
      out.write(&seq[chunk_pos], seq.size() - chunk_pos);
      out << "\n";
    }

    return out.str();
  }

} // namespace OpenMS
