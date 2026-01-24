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

#include <limits>
#include <regex>
#include <sstream>

namespace OpenMS
{
  using namespace std;

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
        // Try to find the modification by accession
        const ResidueModification* res_mod = nullptr;

        if (mod.accession.hasPrefix("UNIMOD:"))
        {
          // Extract UNIMOD ID
          String unimod_id = mod.accession.substr(7);
          res_mod = ModificationsDB::getInstance()->getModification(mod.name, seq[idx].getOneLetterCode(), ResidueModification::ANYWHERE);
        }
        else if (mod.accession.hasPrefix("MOD:"))
        {
          // Try to find by PSI-MOD accession
          res_mod = ModificationsDB::getInstance()->getModification(mod.name, seq[idx].getOneLetterCode(), ResidueModification::ANYWHERE);
        }
        else if (!mod.name.empty())
        {
          // Try to find by name
          res_mod = ModificationsDB::getInstance()->getModification(mod.name, seq[idx].getOneLetterCode(), ResidueModification::ANYWHERE);
        }

        if (res_mod != nullptr)
        {
          seq.setModification(idx, res_mod->getFullId());
        }
      }
      catch (const Exception::BaseException&)
      {
        // Skip if modification cannot be resolved
        OPENMS_LOG_WARN << "Could not resolve modification: " << mod.accession << " (" << mod.name << ")" << std::endl;
      }
    }

    return seq;
  }

  std::vector<std::pair<String, AASequence>> PEFFEntry::getVariantSequences() const
  {
    std::vector<std::pair<String, AASequence>> variants;

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
      if (!var.sources.empty())
      {
        desc += " (" + var.sources + ")";
      }

      try
      {
        variants.push_back(std::make_pair(desc, AASequence::fromString(variant_seq)));
      }
      catch (const Exception::BaseException&)
      {
        // Skip if sequence cannot be parsed
      }
    }

    return variants;
  }

  AASequence PEFFEntry::getProcessedSequence(const String& region_type) const
  {
    // Find the first processed region of the given type
    for (const auto& region : processed_regions)
    {
      if (region.type == region_type)
      {
        // For signal peptides (PEFF:0001021), return the mature protein (after the signal)
        if (region_type == "PEFF:0001021" || region.name == "signal peptide")
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
    infile_.seekg(0, infile_.end);
    fileSize_ = infile_.tellg();
    infile_.seekg(0, infile_.beg);

    headers_.clear();
    entries_read_ = 0;

    // Parse header section
    std::streambuf* sb = infile_.rdbuf();
    std::string line;
    bool in_header = true;
    PEFFDatabaseMetadata current_header;
    bool has_header = false;

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
          // The # // line marks end of header section for current database
          // Push the current header if we have one
          if (has_header)
          {
            headers_.push_back(current_header);
            current_header = PEFFDatabaseMetadata();
            has_header = false;
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
          headers_.push_back(current_header);
        }
      }
    }

    // If no headers found, create a default one
    if (headers_.empty())
    {
      headers_.push_back(PEFFDatabaseMetadata());
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

    outfile_.open(filename.c_str(), ofstream::out);

    if (!outfile_.good())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    // Write header
    outfile_ << formatHeader_(header);
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
    // Simple ProForma conversion - returns sequence only for now
    if (entry.modifications.empty())
    {
      return entry.sequence;
    }

    // TODO: Implement full ProForma conversion with modifications
    // This would require inserting modification annotations at the correct positions
    OPENMS_LOG_WARN << "toProForma: Modifications are not yet encoded in ProForma output" << std::endl;
    return entry.sequence;
  }

  void PEFFFile::parseHeaderLine_(const String& line, PEFFDatabaseMetadata& header, bool& new_db)
  {
    new_db = false;

    // Line includes the leading # character
    String trimmed(line);
    trimmed.trim();

    // Check for database separator
    if (trimmed == "#" || trimmed == "# //" || trimmed == "#//")
    {
      new_db = true;
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
    else if (key == "DbDate")
    {
      header.db_date = value;
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
      header.general_comment = value;
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
    else if (key == "OptionalTagDef")
    {
      header.optional_tag_defs.push_back(value);
    }
  }

  void PEFFFile::parseAnnotations_(const String& description, PEFFEntry& entry)
  {
    if (description.empty())
    {
      return;
    }

    // Split by backslash to find annotation blocks
    std::vector<String> parts;
    size_t pos = 0;
    size_t prev = 0;

    // First part before any backslash is not an annotation
    pos = description.find('\\');
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

    // Parse each annotation
    prev = pos;
    while (prev < description.size())
    {
      // Skip the backslash
      prev++;

      // Find the next backslash or end
      pos = description.find('\\', prev);
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
      else if (key == "ModRes" || key == "ModResPsi" || key == "ModResUnimod")
      {
        std::vector<String> mods = parseParenList_(value);
        for (const String& mod : mods)
        {
          entry.modifications.push_back(parseModification_(mod));
        }
      }
      else if (key == "VariantSimple")
      {
        std::vector<String> variants = parseParenList_(value);
        for (const String& var : variants)
        {
          entry.simple_variants.push_back(parseVariantSimple_(var));
        }
      }
      else if (key == "VariantComplex")
      {
        std::vector<String> variants = parseParenList_(value);
        for (const String& var : variants)
        {
          entry.complex_variants.push_back(parseVariantComplex_(var));
        }
      }
      else if (key == "Processed")
      {
        std::vector<String> regions = parseParenList_(value);
        for (const String& reg : regions)
        {
          entry.processed_regions.push_back(parseProcessedRegion_(reg));
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
    size_t pos = 0;
    while (pos < value.size())
    {
      if (value[pos] == '(')
      {
        // Find matching closing parenthesis
        int depth = 1;
        size_t start = pos + 1;
        size_t end = start;

        while (end < value.size() && depth > 0)
        {
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
    // Format: pos|accession|name|evidence  (evidence is optional)
    // Position can be ? for unknown position
    PEFFModification mod;

    std::vector<String> parts;
    tuple.split('|', parts);

    if (parts.size() >= 2)
    {
      // Position
      if (parts[0] == "?" || parts[0].empty())
      {
        mod.position = 0;
      }
      else
      {
        mod.position = parts[0].toInt();
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

      // Evidence
      if (parts.size() >= 4)
      {
        mod.evidence = parts[3];
      }
    }

    return mod;
  }

  PEFFVariantSimple PEFFFile::parseVariantSimple_(const String& tuple)
  {
    // Format: pos|aa|sources (sources optional)
    PEFFVariantSimple var;

    std::vector<String> parts;
    tuple.split('|', parts);

    if (parts.size() >= 2)
    {
      var.position = parts[0].toInt();
      var.variant_aa = parts[1].empty() ? '\0' : parts[1][0];

      if (parts.size() >= 3)
      {
        var.sources = parts[2];
      }
    }

    return var;
  }

  PEFFVariantComplex PEFFFile::parseVariantComplex_(const String& tuple)
  {
    // Format: start|end|replacement|sources (sources optional)
    PEFFVariantComplex var;

    std::vector<String> parts;
    tuple.split('|', parts);

    if (parts.size() >= 3)
    {
      var.start_position = parts[0].toInt();
      var.end_position = parts[1].toInt();
      var.replacement = parts[2];

      if (parts.size() >= 4)
      {
        var.sources = parts[3];
      }
    }

    return var;
  }

  PEFFProcessedRegion PEFFFile::parseProcessedRegion_(const String& tuple)
  {
    // Format: start|end|type|name|desc (name and desc optional)
    PEFFProcessedRegion reg;

    std::vector<String> parts;
    tuple.split('|', parts);

    if (parts.size() >= 3)
    {
      reg.start_position = parts[0].toInt();
      reg.end_position = parts[1].toInt();
      reg.type = parts[2];

      if (parts.size() >= 4)
      {
        reg.name = parts[3];
      }

      if (parts.size() >= 5)
      {
        reg.description = parts[4];
      }
    }

    return reg;
  }

  String PEFFFile::formatHeader_(const PEFFDatabaseMetadata& header) const
  {
    std::ostringstream out;

    out << "# PEFF " << header.version << "\n";

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

    if (!header.db_date.empty())
    {
      out << "# DbDate=" << header.db_date << "\n";
    }

    if (header.number_of_entries > 0)
    {
      out << "# NumberOfEntries=" << header.number_of_entries << "\n";
    }

    out << "# SequenceType=" << (header.sequence_type == PEFFDatabaseMetadata::SequenceType::AA ? "AA" : "NA") << "\n";

    if (!header.general_comment.empty())
    {
      out << "# GeneralComment=" << header.general_comment << "\n";
    }

    for (const auto& [key, desc] : header.specific_keys)
    {
      out << "# SpecificKey=" << key << ":" << desc << "\n";
    }

    for (const String& tag_def : header.optional_tag_defs)
    {
      out << "# OptionalTagDef=" << tag_def << "\n";
    }

    out << "# //\n";

    return out.str();
  }

  String PEFFFile::formatEntry_(const PEFFEntry& entry) const
  {
    std::ostringstream out;

    // Start with identifier
    out << ">" << entry.identifier;

    // Build description with annotations
    std::ostringstream desc;

    // Protein names
    if (!entry.protein_names.empty())
    {
      desc << " \\PName=";
      if (entry.protein_names.size() == 1)
      {
        desc << "(" << entry.protein_names[0] << ")";
      }
      else
      {
        for (const String& name : entry.protein_names)
        {
          desc << "(" << name << ")";
        }
      }
    }

    // Gene name
    if (!entry.gene_name.empty())
    {
      desc << " \\GName=" << entry.gene_name;
    }

    // NCBI Tax ID
    if (entry.ncbi_tax_id != 0)
    {
      desc << " \\NcbiTaxId=" << entry.ncbi_tax_id;
    }

    // Taxonomy name
    if (!entry.taxonomy_name.empty())
    {
      desc << " \\TaxName=" << entry.taxonomy_name;
    }

    // Length
    if (entry.sequence_length > 0)
    {
      desc << " \\Length=" << entry.sequence_length;
    }

    // Sequence version
    if (!entry.sequence_version.empty())
    {
      desc << " \\SV=" << entry.sequence_version;
    }

    // Entry version
    if (!entry.entry_version.empty())
    {
      desc << " \\EV=" << entry.entry_version;
    }

    // Protein existence
    if (entry.protein_existence > 0)
    {
      desc << " \\PE=" << entry.protein_existence;
    }

    // Modifications
    if (!entry.modifications.empty())
    {
      desc << " \\ModRes=";
      for (const PEFFModification& mod : entry.modifications)
      {
        desc << "(";
        if (mod.position == 0)
        {
          desc << "?";
        }
        else
        {
          desc << mod.position;
        }
        desc << "|" << mod.accession;
        if (!mod.name.empty())
        {
          desc << "|" << mod.name;
        }
        if (!mod.evidence.empty())
        {
          desc << "|" << mod.evidence;
        }
        desc << ")";
      }
    }

    // Simple variants
    if (!entry.simple_variants.empty())
    {
      desc << " \\VariantSimple=";
      for (const PEFFVariantSimple& var : entry.simple_variants)
      {
        desc << "(" << var.position << "|" << var.variant_aa;
        if (!var.sources.empty())
        {
          desc << "|" << var.sources;
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
        desc << "(" << var.start_position << "|" << var.end_position << "|" << var.replacement;
        if (!var.sources.empty())
        {
          desc << "|" << var.sources;
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
        desc << "(" << reg.start_position << "|" << reg.end_position << "|" << reg.type;
        if (!reg.name.empty())
        {
          desc << "|" << reg.name;
        }
        if (!reg.description.empty())
        {
          desc << "|" << reg.description;
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
        desc << "(" << pf << ")";
      }
    }

    // Custom annotations
    for (const auto& [key, value] : entry.custom_annotations)
    {
      desc << " \\" << key << "=" << value;
    }

    out << desc.str() << "\n";

    // Write sequence (80 characters per line)
    const String& seq = entry.sequence;
    int chunks = seq.size() / 80;
    Size chunk_pos = 0;

    for (int i = 0; i < chunks; ++i)
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

  String PEFFFile::formatModifications_(const std::vector<PEFFModification>& mods, const String& key) const
  {
    if (mods.empty())
    {
      return "";
    }

    std::ostringstream out;
    out << " \\" << key << "=";

    for (const PEFFModification& mod : mods)
    {
      out << "(";
      if (mod.position == 0)
      {
        out << "?";
      }
      else
      {
        out << mod.position;
      }
      out << "|" << mod.accession;
      if (!mod.name.empty())
      {
        out << "|" << mod.name;
      }
      if (!mod.evidence.empty())
      {
        out << "|" << mod.evidence;
      }
      out << ")";
    }

    return out.str();
  }

} // namespace OpenMS
