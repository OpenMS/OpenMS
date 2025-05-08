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
    ProForma::ProForma(const AASequence& seq) : sequence_(seq) {}

    void ProForma::validateCVModification(const std::string& modification)
    {
        size_t colon_pos = modification.find(':');
        if (colon_pos == std::string::npos)
        {
            std::cout << "No CV prefix found in modification: " << modification << std::endl;
            return;
        }

        std::string cv = modification.substr(0, colon_pos);
        std::string accession = modification.substr(colon_pos + 1);

        std::cout << "Validating CV: " << cv << " with accession: " << accession << std::endl;

        if (supported_cvs_.find(cv) == supported_cvs_.end())
        {
            std::cout << "Unsupported CV detected: " << cv << std::endl;
            throw std::invalid_argument("Unsupported CV/ontology: " + cv);
        }

        if (accession.empty())
        {
            throw std::invalid_argument("Accession number cannot be empty in modification: " + modification);
        }
    }

    void ProForma::parseCVModificationNames(const std::string& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for CV modification.");
        }

        std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        std::cout << "Parsing CV modification: " << modification << " at position " << residue_pos << std::endl;

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

    void ProForma::parseStandardModification(const std::string& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid modification format: Missing brackets for standard modification.");
        }

        std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        std::cout << "Parsing standard modification: " << modification << " at position " << residue_pos << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    void ProForma::parseDeltaMassNotation(const std::string& modString, size_t& pos, size_t residue_pos)
    {
      size_t modStart = modString.find('[', pos);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == std::string::npos || modEnd == std::string::npos)
      {
        throwParseError("Invalid mass shift notation: Missing brackets.");
      }

      std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      std::cout << "Parsing mass shift: " << modification << " at position " << residue_pos << std::endl;

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


    void ProForma::parseNTerminalModification(const std::string& modString, size_t& pos)
    {
        size_t modEnd = modString.find("]-", pos);
        if (modEnd == std::string::npos)
        {
            throwParseError("Invalid N-terminal modification format: Missing brackets and '-' indicator.");
        }

        std::string modification = modString.substr(1, modEnd - 1);
        std::cout << "Parsing N-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[0] = attributes;

        pos = modEnd + 2;
    }

    void ProForma::parseCTerminalModification(const std::string& modString, size_t& pos)
    {
        size_t modStart = modString.find("-[", pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == std::string::npos || modEnd == std::string::npos)
        {
            throwParseError("Invalid C-terminal modification format: Missing brackets and '-' indicator.");
        }

        std::string modification = modString.substr(modStart + 2, modEnd - modStart - 2);
        std::cout << "Parsing C-terminal modification: " << modification << std::endl;

        ModificationAttributes attributes;
        attributes.modification_name = modification;
        modifications_[sequence_.size() + 1] = attributes;

        pos = modEnd + 1;
    }

    void ProForma::parseRangeModification(const std::string& modString, size_t& pos)
    {
      size_t rangeStart = modString.find('(', pos);
      size_t rangeEnd = modString.find(')', rangeStart);
      if (rangeStart == std::string::npos || rangeEnd == std::string::npos)
      {
        throwParseError("Invalid range format: Missing parentheses.");
      }

      // Extract the range of residues (e.g., "ESFRMS" in "PRT(ESFRMS)[+19.0523]ISK")
      std::string range = modString.substr(rangeStart + 1, rangeEnd - rangeStart - 1);
      std::cout << "Parsing range modification: " << range << std::endl;

      size_t modStart = modString.find('[', rangeEnd);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == std::string::npos || modEnd == std::string::npos)
      {
        throwParseError("Invalid modification format: Missing brackets after range.");
      }

      std::string modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      std::cout << "Applying modification: " << modification << " to range: " << range << std::endl;

      // Calculate the exact start position in the sequence for the range
      size_t range_seq_start = 0;
      for (size_t i = 0; i < sequence_.size(); ++i)
      {
        // Find the first occurrence of the range in the sequence (i.e., match "ESFRMS")
        std::string substring = sequence_.toString().substr(i, range.size());
        if (substring == range)
        {
          range_seq_start = i + 1; // +1 because the range starts from 1
          break;
        }
      }

      if (range_seq_start == 0)
      {
        throwParseError("Range not found in sequence.");
      }

      // Apply the modification to the entire range
      ModificationAttributes attributes;
      attributes.modification_name = modification;
      attributes.range = {range_seq_start, range_seq_start + range.size() - 1};  // Store the range as a pair (start, end)

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

      // If there's additional sequence after the modification (like "ISK"), we need to ensure it's processed properly.
      if (pos < modString.size())
      {
        std::string remaining_seq = modString.substr(pos);
        // Add the remaining residues (ISK in this case) to the sequence if it hasn't already been processed
        size_t expected_remaining_start = range_seq_start + range.size();

        if (remaining_seq.size() + expected_remaining_start - 1 <= sequence_.size())
        {
          // Verify that the next part of the string matches the unprocessed part of the sequence
          std::string sequence_part = sequence_.toString().substr(expected_remaining_start - 1);
          if (remaining_seq != sequence_part)
          {
            sequence_ = AASequence::fromString(sequence_.toString().substr(0, expected_remaining_start) + remaining_seq);
          }
        }

        pos = modString.size(); // Set position to the end of the string
      }
    }

    // UPDATED FUNCTION TO HANDLE RANGES
    AASequence ProForma::fromProFormaString(const std::string& proforma_str)
    {
      AASequence seq;
      size_t pos = 0;
      size_t residue_pos = 0;

      while (pos < proforma_str.size())
      {
        try
        {
          // Handle amino acid residues
          if (std::isalpha(proforma_str[pos]))
          {
            seq = seq + ResidueDB::getInstance()->getResidue(proforma_str.substr(pos, 1));
            residue_pos = seq.size();
            pos++;
          }
          // Handle modifications (PSI-MOD, UNIMOD, RESID, etc.)
          else if (proforma_str[pos] == '[')
          {
            if (proforma_str.find(':', pos) != std::string::npos)
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
          else if (proforma_str[pos] == '[' && proforma_str.find("]-") != std::string::npos)
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
        catch (const std::runtime_error& e)
        {
          throw std::runtime_error("Error parsing ProForma string: " + std::string(e.what()));
        }
      }
      return seq;
    }


    std::string ProForma::toProFormaString() const
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
        }

        // Close range if applicable
        if (in_range && i == range_end)
        {
          ss << ")";  // Close range after the last residue in the range
          String mod_name_str = mod_name.str();
          mod_name.clear();
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
    void ProForma::throwParseError(const std::string& message) const
    {
      throw std::runtime_error("ProForma parsing error: " + message);
    }


    void ProForma::removeModification(size_t position)
    {
      std::cout << "Attempting to remove modification at position: " << position << std::endl;

      auto it = modifications_.find(position);
      if (it != modifications_.end())
      {
        std::cout << "Removing modification: " << it->second.modification_name << " at position: " << position << std::endl;
        modifications_.erase(it);
      }
      else
      {
        std::cerr << "No modification found at position: " << position << std::endl;
      }
    }

    // N term mod : start_pos = 0. Modification after the first a.a. : start_pos = 1
    void ProForma::addModification(size_t start_pos, size_t end_pos, const std::string& mod_id, double mass_shift)
    {
      modifications_[start_pos] = {mass_shift, false, false, mod_id, std::pair<size_t, size_t>(start_pos, end_pos)};
    }
    }