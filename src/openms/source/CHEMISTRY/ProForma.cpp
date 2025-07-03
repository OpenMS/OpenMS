// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace OpenMS
{
    void ProForma::validateCVModification(const String& modification)
    {
        size_t colon_pos = modification.find(':');
        if (colon_pos == String::npos)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
             "No CV prefix found in modification: " + modification);
        }

        String cv = modification.substr(0, colon_pos);
        String accession = modification.substr(colon_pos + 1);

        OPENMS_LOG_DEBUG << "Validating CV: " << cv << " with accession: " << accession << std::endl;

        if (supported_cvs_.find(cv) == supported_cvs_.end())
        {
            OPENMS_LOG_ERROR << "Unsupported CV detected: " << cv << std::endl;
            throw Exception::IllegalArgument(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Unsupported CV/ontology: " + cv);
        }

        if (accession.empty())
        {
            throw Exception::IllegalArgument(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Accession number cannot be empty in modification: " + modification);
        }
    }

    void ProForma::parseCVModificationNames(const String& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for CV modification.");
        }

        String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        OPENMS_LOG_DEBUG << "Parsing CV modification: " << modification << " at position " << residue_pos << std::endl;

        validateCVModification(modification);

        ModificationAttributes attributes;
        attributes.modification_name = modification;

        if (modString.size() > modEnd + 1 && modString[modEnd + 1] == '?')
        {
            attributes.ambiguous_start = true;
            modEnd++;
        }

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    void ProForma::parseStandardModification(const String& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for standard modification.");
        }

        String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        OPENMS_LOG_DEBUG << "Parsing standard modification: " << modification << " at position " << residue_pos << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    void ProForma::parseDeltaMassNotation(const String& modString, size_t& pos, size_t residue_pos)
    {
      size_t modStart = modString.find('[', pos);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == String::npos || modEnd == String::npos)
      {
        throwParseError("Invalid mass shift notation: Missing brackets.");
      }

      String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      OPENMS_LOG_DEBUG << "Parsing mass shift: " << modification << " at position " << residue_pos << std::endl;

      if (modification[0] != '+' && modification[0] != '-')
      {
        throwParseError("Invalid mass shift format: Missing +/- sign.");
      }

      try
      {
        double mass_shift = std::stod(modification);  // Convert to double
        ModificationAttributes attributes;
        attributes.mass_shift = mass_shift;  // Store the mass shift in attributes

        // Check if the modification is ambiguous
        if (modString.size() > modEnd + 1 && modString[modEnd + 1] == '?')
        {
          attributes.ambiguous_start = true;
          modEnd++;
        }

        modifications_[residue_pos] = attributes;  // Add to modifications map
      }
      catch (const std::invalid_argument&)
      {
        throwParseError("Invalid mass shift format: Could not convert to double.");
      }

      pos = modEnd + 1;
    }


    void ProForma::parseNTerminalModification(const String& modString, size_t& pos)
    {
        size_t modEnd = modString.find("]-", pos);
        if (modEnd == String::npos)
        {
            throwParseError("Invalid N-terminal modification format: Missing brackets and '-' indicator.");
        }

        String modification = modString.substr(1, modEnd - 1);
        OPENMS_LOG_DEBUG << "Parsing N-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[0] = attributes;

        pos = modEnd + 2;
    }

    void ProForma::parseCTerminalModification(const String& modString, size_t& pos)
    {
        size_t modStart = modString.find("-[", pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
            throwParseError("Invalid C-terminal modification format: Missing brackets and '-' indicator.");
        }

        String modification = modString.substr(modStart + 2, modEnd - modStart - 2);
        OPENMS_LOG_DEBUG << "Parsing C-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[sequence_.size() + 1] = attributes;

        pos = modEnd + 1;
    }

    void ProForma::parseRangeModification(const String& modString, size_t& pos)
    {
      size_t rangeStart = modString.find('(', pos);
      size_t rangeEnd = modString.find(')', rangeStart);
      if (rangeStart == String::npos || rangeEnd == String::npos)
      {
        throwParseError("Invalid range format: Missing parentheses.");
      }

      // Extract the range of residues (e.g., "ESFRMS" in "PRT(ESFRMS)[+19.0523]ISK")
      String range = modString.substr(rangeStart + 1, rangeEnd - rangeStart - 1);
      OPENMS_LOG_DEBUG << "Parsing range modification: " << range << std::endl;

      size_t modStart = modString.find('[', rangeEnd);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == String::npos || modEnd == String::npos)
      {
        throwParseError("Invalid modification format: Missing brackets after range.");
      }

      String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      OPENMS_LOG_DEBUG << "Applying modification: " << modification << " to range: " << range << std::endl;

      // The range starts after the sequence that we have already parsed
      size_t range_seq_start = sequence_.size() + 1; // Default to after the sequence end

      OPENMS_LOG_DEBUG << "Range starts at sequence position: " << range_seq_start << std::endl;

      sequence_ += AASequence::fromString(range); // Append the range to the sequence

      // Apply the modification to the entire range
      ModificationAttributes attributes;
      attributes.modification_name = modification;
      attributes.range = {range_seq_start, range_seq_start + range.size() - 1};  // Store the range as a pair (start, end)
      OPENMS_LOG_DEBUG << "Range of modification: " << attributes.range.first << " - " << attributes.range.second << std::endl;
      // Apply the modification from range start to range end
      for (size_t i = range_seq_start; i <= range_seq_start + range.size() - 1; ++i)
      {
        if (i > sequence_.size())
        {
          throw std::out_of_range("Modifying out of sequence bounds");
        }
        modifications_[i] = attributes; // Apply the same modification to all positions in the range
      }

      // Update the current position to after the closing bracket of the modification
      pos = modEnd + 1;
    }

    // UPDATED FUNCTION TO HANDLE RANGES
    void ProForma::fromProFormaString(const String& proforma_str)
    {
      //AASequence seq;
      size_t pos = 0;
      size_t residue_pos = 0;
      modifications_.clear(); // Clear previous modifications
      sequence_ = AASequence(); // Clear previous sequence

      while (pos < proforma_str.size())
      {
        // Handle amino acid residues
        if (std::isalpha(proforma_str[pos]))
        {
          sequence_ += ResidueDB::getInstance()->getResidue(proforma_str.substr(pos, 1));
          residue_pos = sequence_.size();
          pos++;
        }
        // Handle modifications (PSI-MOD, UNIMOD, RESID, etc.)
        else if (proforma_str[pos] == '[')
        {
          if (proforma_str.find(':', pos) != String::npos)
          {
            parseCVModificationNames(proforma_str, pos, residue_pos);
          }
          else if (proforma_str[pos + 1] == '+' || proforma_str[pos + 1] == '-')
          {
            parseDeltaMassNotation(proforma_str, pos, residue_pos);
          }
          else
          {
            parseStandardModification(proforma_str, pos, residue_pos);
          }
        }
        // Handle N-terminal and C-terminal modifications
        else if (proforma_str[pos] == '-' && proforma_str[pos + 1] == '[')
        {
          parseCTerminalModification(proforma_str, pos);
        }
        else if (proforma_str[pos] == '[' && proforma_str.find("]-") != String::npos)
        {
          parseNTerminalModification(proforma_str, pos);
        }
        // Handle range modifications
        else if (proforma_str[pos] == '(')
        {
          parseRangeModification(proforma_str, pos);  // Range modification handling
        }
        else
        {
          pos++;
        }
      }
    }


    String ProForma::toProFormaString() const
    {
      std::stringstream ss;

      // Handle N-terminal modification
      if (modifications_.find(0) != modifications_.end())
      {
        ss << "[" << modifications_.at(0).modification_name << "]-";
      }

      bool in_range = false;
      //size_t range_start = 0;
      size_t range_end = 0;
      std::stringstream mod_name;

      for (size_t i = 1; i <= sequence_.size(); ++i)
      {
        auto it = modifications_.find(i);

        // Open range if applicable
        if (it != modifications_.end() && it->second.range.first == i && it->second.range.second != i && !in_range)
        {
          //range_start = it->second.range.first;
          range_end = it->second.range.second;
          ss << "(";  // Open range
          in_range = true;
        }
        // Add the current residue
        ss << sequence_.getResidue(i - 1).getOneLetterCode();

        if (it != modifications_.end())
        {
          if (!it->second.modification_name.empty())
          {
            mod_name << "[" << it->second.modification_name << "]";
          }
          else if (it->second.mass_shift != 0.0)
          {
            mod_name << "[" << (it->second.mass_shift > 0 ? "+" : "") << it->second.mass_shift << "]";
          }

          if (it->second.range.first == it->second.range.second)
          {
            ss << mod_name.str();
            mod_name.clear();
          }

          if (it->second.ambiguous_start)
          {
            ss << "?";
          }
          mod_name.str("");
          mod_name.clear();
        }

        // Close range if applicable
        if (in_range && i == range_end)
        {
          ss << ")";  // Close range after the last residue in the range
          String mod_name_str = it->second.modification_name;
          if (!mod_name_str.empty())
          {
            ss << "[" << mod_name_str << "]";
          }
          in_range = false;  // End the range
        }
      }

      // Handle C-terminal modification
      if (modifications_.find(sequence_.size() + 1) != modifications_.end())
      {
        ss << "-[" << modifications_.at(sequence_.size() + 1).modification_name << "]";
      }

      return ss.str();
    }



/////////////
    void ProForma::throwParseError(const String& message) const
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error parsing ProForma string: " + message, "Please check the format and ensure all modifications are correctly specified.");
    }


    void ProForma::removeModification(size_t position)
    {
      OPENMS_LOG_DEBUG << "Attempting to remove modification at position: " << position << std::endl;
      if (position > sequence_.size())
      {
        throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      auto it = modifications_.find(position);
      if (it != modifications_.end())
      {
        OPENMS_LOG_DEBUG << "Removing modification: " << it->second.modification_name << " at position: " << position << std::endl;
        modifications_.erase(it);
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid position", "No modification found at the specified position: " + std::to_string(position));
      }
    }

    // N term mod : start_pos = 0. Modification after the first a.a. : start_pos = 1
    void ProForma::addModification(size_t start_pos, size_t end_pos, const String& mod_id, double mass_shift)
    {
      modifications_[start_pos] = {mass_shift, false, false, mod_id, std::pair<size_t, size_t>(start_pos, end_pos)};
    }
    }