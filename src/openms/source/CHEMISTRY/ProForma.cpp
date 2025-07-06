// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------

/**
 * @file ProForma.cpp
 * @brief Implementation of ProForma peptidoform notation parser
 *
 * This file implements the ProForma class for parsing and generating
 * ProForma notation strings according to the HUPO-PSI ProForma 2.0 specification.
 *
 * The implementation provides a state-machine based parser that processes
 * ProForma strings character by character, handling different syntax elements
 * including amino acids, modifications, and special notations.
 */

#include <OpenMS/CHEMISTRY/ProForma.h>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace OpenMS
{
    /**
     * @brief Validates controlled vocabulary modification syntax and support
     *
     * This method implements ProForma 2.0 section on "Controlled vocabularies"
     * by validating that:
     * 1. The modification follows CV:accession format
     * 2. The CV is among the supported vocabularies (UNIMOD, MOD, RESID, XLMOD, GNO)
     * 3. The accession is not empty
     *
     * @param modification The modification string without square brackets (e.g., "UNIMOD:35")
     *
     * @throws Exception::IllegalArgument If CV format is invalid, CV is unsupported,
     *                                    or accession is empty
     *
     * @note This validation ensures compliance with ProForma 2.0 specification
     *       which requires explicit CV prefixes for ontology-based modifications.
     */
    void ProForma::validateCVModification_(const String& modification)
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

        // Check against supported controlled vocabularies as defined in ProForma 2.0
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

    /**
     * @brief Parse controlled vocabulary (CV) based modifications
     *
     * Implements parsing for ProForma 2.0 "Controlled vocabularies" section.
     * Handles modifications specified using ontology accessions from supported
     * controlled vocabularies.
     *
     * Supported formats:
     * - UNIMOD: [UNIMOD:35] (oxidation)
     * - PSI-MOD: [MOD:00046] (phosphorylation)
     * - RESID: [RESID:AA0037] (phosphoserine)
     * - XLMOD: [XLMOD:02000] (cross-linking modifications)
     * - GNOme: [GNO:G00001] (glycan modifications)
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after modification
     * @param residue_pos 1-based position in sequence where modification applies
     *
     * @throws Exception::ParseError If bracket structure is malformed
     * @throws Exception::IllegalArgument If CV is unsupported or format invalid
     *
     * @note Also handles ambiguous modification markers (?) following the closing bracket
     */
    void ProForma::parseCVModificationNames_(const String& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid modification format: Missing brackets for CV modification.", "Missing brackets for CV modification.");
        }

        String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        OPENMS_LOG_DEBUG << "Parsing CV modification: " << modification << " at position " << residue_pos << std::endl;

        validateCVModification_(modification);

        ResidueModification attributes;

        if (modification.hasPrefix("UNIMOD:"))
        {
            // Handle UNIMOD modifications
            size_t colon_pos = modification.find(':');
            String unimod_id = modification.substr(colon_pos + 1);
            OPENMS_LOG_DEBUG << "UNIMOD ID: " << unimod_id << std::endl;
            attributes.setId(unimod_id);
            modifications_[residue_pos] = attributes;
        }
        else if (modification.hasPrefix("MOD:"))
        {
            // Handle MOD modifications
            size_t colon_pos = modification.find(':');
            String mod_id = modification.substr(colon_pos + 1);
            OPENMS_LOG_DEBUG << "MOD ID: " << mod_id << std::endl;
            attributes.setId(mod_id);
            modifications_[residue_pos] = attributes;
        }
        else if (modification.hasPrefix("RESID:"))
        {
            // Handle RESID modifications
            size_t colon_pos = modification.find(':');
            String resid_id = modification.substr(colon_pos + 1);
            OPENMS_LOG_DEBUG << "RESID ID: " << resid_id << std::endl;
            attributes.setRESIDAccession(resid_id);
            modifications_[residue_pos] = attributes;
        }
        else if (modification.hasPrefix("XLMOD:"))
        {
            // Handle XLMOD modifications
            size_t colon_pos = modification.find(':');
            String xlm_id = modification.substr(colon_pos + 1);
            OPENMS_LOG_DEBUG << "XLMOD ID: " << xlm_id << std::endl;
            attributes.setXLMODAccession(xlm_id);
            modifications_[residue_pos] = attributes;
        }
        else if (modification.hasPrefix("GNO:"))
        {
            // Handle GNO modifications
            size_t colon_pos = modification.find(':');
            String gno_id = modification.substr(colon_pos + 1);
            OPENMS_LOG_DEBUG << "GNO ID: " << gno_id << std::endl;
            attributes.setId(gno_id);
            modifications_[residue_pos] = attributes;
        }
        else
        {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unsupported CV/ontology in modification: " + modification, "Unsupported CV/ontology in modification.");
        }
        attributes.setName(modification);

        if (modString[modEnd + 1] == '?')
        {
            attributes.setAmbiguous(true);
            modEnd++;
        }

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    /**
     * @brief Parse named modifications without controlled vocabulary prefixes
     *
     * Handles user-defined modification names that don't use controlled vocabulary
     * notation. These are free-text modification names enclosed in square brackets.
     *
     * Examples:
     * - [Phospho] - phosphorylation
     * - [Acetyl] - acetylation
     * - [Oxidation] - oxidation
     * - [MyCustomMod] - user-defined modification
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after modification
     * @param residue_pos 1-based position in sequence where modification applies
     *
     * @throws Exception::ParseError If bracket structure is malformed
     *
     * @note This method handles the most basic form of ProForma modification notation
     *       and is commonly used for well-known modifications with standard names.
     */
    void ProForma::parseStandardModification_(const String& modString, size_t& pos, size_t residue_pos)
    {
        size_t modStart = modString.find('[', pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid modification format: Missing brackets for standard modification.", "Missing brackets for standard modification.");
        }

        String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
        OPENMS_LOG_DEBUG << "Parsing standard modification: " << modification << " at position " << residue_pos << std::endl;

        ResidueModification attributes;
        attributes.setName(modification);

        modifications_[residue_pos] = attributes;
        pos = modEnd + 1;
    }

    /**
     * @brief Parse delta mass notation modifications
     *
     * Implements ProForma 2.0 "Delta mass" notation for modifications specified
     * as mass shifts. This is useful when the exact modification is unknown but
     * the mass difference is measured.
     *
     * Format requirements (ProForma 2.0 compliant):
     * - Must start with explicit + or - sign
     * - Must be a valid decimal number
     * - Mass in Daltons (Da)
     * - Optional ? for ambiguous modifications
     *
     * Examples:
     * - [+15.9949] - oxidation mass shift
     * - [-17.0265] - loss of ammonia
     * - [+79.9663] - phosphorylation mass shift
     * - [+42.0106]? - ambiguous acetylation
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after modification
     * @param residue_pos 1-based position in sequence where modification applies
     *
     * @throws Exception::ParseError If brackets missing, sign missing, or invalid number format
     *
     * @note This parser enforces ProForma 2.0 requirement for explicit +/- signs
     *       to distinguish mass shifts from other modification types.
     */
    void ProForma::parseDeltaMassNotation_(const String& modString, size_t& pos, size_t residue_pos)
    {
      size_t modStart = modString.find('[', pos);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == String::npos || modEnd == String::npos)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid mass shift notation: Missing brackets.", "Missing brackets.");
      }

      String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      OPENMS_LOG_DEBUG << "Parsing mass shift: " << modification << " at position " << residue_pos << std::endl;

      if (modification[0] != '+' && modification[0] != '-')
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid mass shift format: Missing +/- sign.", "Invalid mass shift format.");
      }

      try
      {
        double mass_shift = modification.toDouble();  // Convert to double
        ResidueModification attributes;
        attributes.setDiffMonoMass(mass_shift);  // Store the mass shift in attributes

        // Check if the modification is ambiguous
        if (modString[modEnd + 1] == '?')
        {
          attributes.setAmbiguous(true);
          modEnd++;
        }

        modifications_[residue_pos] = attributes;  // Add to modifications map
      }
      catch (const Exception::ConversionError&)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid mass shift format: Could not convert to double.", "Invalid mass shift format.");
      }

      pos = modEnd + 1;
    }


    /**
     * @brief Parse N-terminal modifications
     *
     * Implements ProForma 2.0 "Terminal modifications" section for modifications
     * at the peptide N-terminus. These modifications affect the amino group of
     * the first amino acid or the peptide as a whole.
     *
     * Required format: [modification]-sequence
     * - Opening bracket immediately followed by modification name
     * - Closing bracket followed by dash
     * - Sequence follows the dash
     *
     * Examples:
     * - [Acetyl]-PEPTIDE - N-terminal acetylation
     * - [Carbamyl]-MSEQUENCE - N-terminal carbamylation
     * - [TMT6plex]-PEPTIDE - TMT labeling at N-terminus
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after "]-" sequence
     *
     * @throws Exception::ParseError If the expected "]-" pattern is not found
     *
     * @note N-terminal modifications are stored at position 0 in the modifications map
     *       to distinguish them from modifications on the first amino acid (position 1).
     */
    void ProForma::parseNTerminalModification_(const String& modString, size_t& pos)
    {
        size_t modEnd = modString.find("]-", pos);
        if (modEnd == String::npos)
        {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid N-terminal modification format:  " + modString, "Missing brackets and '-' indicator.");
        }

        String modification = modString.substr(1, modEnd - 1);  // Skip opening '[', extract to closing ']'
        OPENMS_LOG_DEBUG << "Parsing N-terminal modification: " << modification << std::endl;

        ResidueModification attributes;
        attributes.setName(modification);
        modifications_[0] = attributes;  // Position 0 = N-terminus

        pos = modEnd + 2;  // Move past the "]-" sequence
    }

    /**
     * @brief Parse C-terminal modifications
     *
     * Implements ProForma 2.0 "Terminal modifications" section for modifications
     * at the peptide C-terminus. These modifications affect the carboxyl group of
     * the last amino acid or the peptide as a whole.
     *
     * Required format: sequence-[modification]
     * - Sequence followed by dash
     * - Opening bracket with modification name
     * - Closing bracket
     *
     * Examples:
     * - PEPTIDE-[Amidation] - C-terminal amidation
     * - SEQUENCE-[Methylester] - C-terminal methylation
     * - PEPTIDE-[Deamidation] - C-terminal deamidation
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after closing bracket
     *
     * @throws Exception::ParseError If the expected "-[...]" pattern is not found
     *
     * @note C-terminal modifications are stored at position (sequence_length + 1)
     *       in the modifications map to distinguish them from modifications on
     *       the last amino acid.
     */
    void ProForma::parseCTerminalModification_(const String& modString, size_t& pos)
    {
        size_t modStart = modString.find("-[", pos);
        size_t modEnd = modString.find(']', modStart);
        if (modStart == String::npos || modEnd == String::npos)
        {
           throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid C-terminal modification format:  " + modString, "Missing brackets and '-' indicator.");
        }

        String modification = modString.substr(modStart + 2, modEnd - modStart - 2);  // Extract between "-[" and "]"
        OPENMS_LOG_DEBUG << "Parsing C-terminal modification: " << modification << std::endl;

        ResidueModification attributes;
        attributes.setName(modification);
        modifications_[sequence_.size() + 1] = attributes;  // Position after sequence end = C-terminus

        pos = modEnd + 1;  // Move past the closing bracket
    }

    /**
     * @brief Parse range modifications spanning multiple residues
     *
     * Implements parsing for modifications that span multiple consecutive amino acids.
     *
     * Required format: prefix(range)[modification]suffix
     * - Parentheses enclose the amino acid sequence affected by the modification
     * - Square brackets follow with the modification specification
     * - The range becomes part of the main sequence
     *
     * Examples:
     * - A(CDE)[+12.5]FGH - modification affecting residues CDE as a unit
     * - PEP(TIDE)[Oxidation] - oxidation affecting the TIDE region
     * - (MSEQ)[+15.99]UENCE - N-terminal region modification
     *
     * @param modString The complete ProForma string being parsed
     * @param pos Current parsing position, updated to character after modification
     *
     * @throws Exception::ParseError If parentheses or brackets are missing/malformed
     * @throws std::out_of_range If range extends beyond sequence boundaries
     *
     * @note Range modifications are applied to all positions within the specified range.
     *       The range sequence is appended to the main sequence during parsing.
     *       Each position in the range gets the same modification attributes.
     */
    void ProForma::parseRangeModification_(const String& modString, size_t& pos)
    {
      size_t rangeStart = modString.find('(', pos);
      size_t rangeEnd = modString.find(')', rangeStart);
      if (rangeStart == String::npos || rangeEnd == String::npos)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid range format:  " + modString, "Missing parentheses.");
      }

      // Extract the range of residues (e.g., "ESFRMS" in "PRT(ESFRMS)[+19.0523]ISK")
      String range = modString.substr(rangeStart + 1, rangeEnd - rangeStart - 1);
      OPENMS_LOG_DEBUG << "Parsing range modification: " << range << std::endl;

      size_t modStart = modString.find('[', rangeEnd);
      size_t modEnd = modString.find(']', modStart);
      if (modStart == String::npos || modEnd == String::npos)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid modification format:  " + modString, "Missing brackets after range.");
      }

      String modification = modString.substr(modStart + 1, modEnd - modStart - 1);
      OPENMS_LOG_DEBUG << "Applying modification: " << modification << " to range: " << range << std::endl;

      // The range starts after the sequence that we have already parsed
      size_t range_seq_start = sequence_.size() + 1; // Default to after the sequence end

      OPENMS_LOG_DEBUG << "Range starts at sequence position: " << range_seq_start << std::endl;

      sequence_ += AASequence::fromString(range); // Append the range to the sequence

      // Apply the modification to the entire range
      ResidueModification attributes;
      attributes.setName(modification);
      attributes.setPositionRange(range_seq_start, range_seq_start + range.size() - 1);  // Store the range as a pair (start, end)
      OPENMS_LOG_DEBUG << "Range of modification: " << attributes.getPositionRange().first << " - " << attributes.getPositionRange().second << std::endl;
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

    /**
     * @brief Main parser function for ProForma notation strings
     *
     * This is the primary entry point for parsing ProForma strings. It implements
     * a state-machine based parser that processes the input character by character,
     * identifying and dispatching to appropriate specialized parsing functions.
     *
     * Parsing algorithm:
     * 1. Clear any existing sequence and modifications
     * 2. Iterate through each character in the input string
     * 3. Identify syntax elements (amino acids, brackets, parentheses, etc.)
     * 4. Dispatch to specialized parsers for each element type
     * 5. Build the final sequence and modification map
     *
     * Supported ProForma 2.0 elements parsed:
     * - Amino acid residues (A-Z including unconventional and ambiguous)
     * - CV modifications ([CV:accession])
     * - Named modifications ([name])
     * - Delta mass modifications ([+/-mass])
     * - N-terminal modifications ([mod]-)
     * - C-terminal modifications (-[mod])
     * - Range modifications ((range)[mod])
     * - Ambiguous markers (?)
     *
     * @param proforma_str The ProForma notation string to parse
     *
     * @throws Exception::ParseError For malformed syntax elements
     * @throws Exception::IllegalArgument For unsupported CV ontologies
     * @throws std::out_of_range For invalid sequence positions
     *
     * @note This method clears any existing sequence and modification data
     *       before parsing the new input string.
     */
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
            parseCVModificationNames_(proforma_str, pos, residue_pos);
          }
          else if (proforma_str[pos + 1] == '+' || proforma_str[pos + 1] == '-')
          {
            parseDeltaMassNotation_(proforma_str, pos, residue_pos);
          }
          else
          {
            parseStandardModification_(proforma_str, pos, residue_pos);
          }
        }
        // Handle N-terminal and C-terminal modifications
        else if (proforma_str[pos] == '-' && proforma_str[pos + 1] == '[')
        {
          parseCTerminalModification_(proforma_str, pos);
        }
        else if (proforma_str[pos] == '[' && proforma_str.find("]-") != String::npos)
        {
          parseNTerminalModification_(proforma_str, pos);
        }
        // Handle range modifications
        else if (proforma_str[pos] == '(')
        {
          parseRangeModification_(proforma_str, pos);  // Range modification handling
        }
        else
        {
          pos++;
        }
      }
    }


    /**
     * @brief Generate ProForma notation string from internal representation
     *
     * Converts the current amino acid sequence and modifications back into
     * a valid ProForma notation string. The output follows ProForma 2.0
     * syntax conventions and can be parsed back by fromProFormaString().
     *
     * Generation algorithm:
     * 1. Add N-terminal modification if present ([mod]-)
     * 2. Iterate through sequence positions
     * 3. Handle range modifications with parentheses
     * 4. Add amino acid residue
     * 5. Add modifications for current position
     * 6. Add ambiguous markers if applicable
     * 7. Close ranges when appropriate
     * 8. Add C-terminal modification if present (-[mod])
     *
     * Output format elements:
     * - Named modifications: [ModName]
     * - Mass shift modifications: [+/-mass]
     * - CV modifications: [CV:accession]
     * - N-terminal: [mod]-sequence
     * - C-terminal: sequence-[mod]
     * - Range modifications: prefix(range)[mod]suffix
     * - Ambiguous: [mod]?
     *
     * @return ProForma notation string representing current peptidoform
     *
     */
    String ProForma::toProFormaString() const
    {
      std::stringstream ss;

      // Handle N-terminal modification
      if (modifications_.find(0) != modifications_.end())
      {
        ss << "[" << modifications_.at(0).getName() << "]-";
      }

      bool in_range = false;
      //size_t range_start = 0;
      size_t range_end = 0;
      std::stringstream mod_name;

      for (size_t i = 1; i <= sequence_.size(); ++i)
      {
        auto it = modifications_.find(i);

        // Open range if applicable
        if (it != modifications_.end() && it->second.getPositionRange().first == i && it->second.getPositionRange().second != i && !in_range)
        {
          //range_start = it->second.range.first;
          range_end = it->second.getPositionRange().second;
          ss << "(";  // Open range
          in_range = true;
        }
        // Add the current residue
        ss << sequence_.getResidue(i - 1).getOneLetterCode();

        if (it != modifications_.end())
        {
          if (!it->second.getName().empty())
          {
            mod_name << "[" << it->second.getName() << "]";
          }
          else if (it->second.getDiffMonoMass() != 0.0)
          {
            mod_name << "[" << (it->second.getDiffMonoMass() > 0 ? "+" : "") << it->second.getDiffMonoMass() << "]";
          }

          if (it->second.getPositionRange().first == it->second.getPositionRange().second)
          {
            ss << mod_name.str();
            mod_name.clear();
          }

          if (it->second.isAmbiguous())
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
          String mod_name_str = it->second.getName();
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
        ss << "-[" << modifications_.at(sequence_.size() + 1).getName() << "]";
        OPENMS_LOG_DEBUG << "Found C-terminal modification at position " << (sequence_.size() + 1) << ": " << modifications_.at(sequence_.size() + 1).getName() << std::endl;
      }
      else
      {
        OPENMS_LOG_DEBUG << "No C-terminal modification found at position " << (sequence_.size() + 1) << std::endl;
        // Check if there's a modification at the last amino acid position that should be treated as C-terminal
        if (modifications_.find(sequence_.size()) != modifications_.end())
        {
          OPENMS_LOG_DEBUG << "Found modification at last amino acid position " << sequence_.size() << ": " << modifications_.at(sequence_.size()).getName() << std::endl;
        }
      }

      return ss.str();
    }

    /**
     * @brief Remove modification at specified position
     *
     * Removes a modification from the internal modifications map at the given position.
     * This method allows for dynamic modification of the peptidoform after parsing.
     *
     * Position mapping:
     * - 0: N-terminal modification
     * - 1 to sequence_length: Amino acid positions (1-based)
     * - sequence_length + 1: C-terminal modification
     *
     * @param position The position from which to remove the modification
     *
     * @throws Exception::OutOfRange If position exceeds sequence length + 1
     * @throws Exception::InvalidValue If no modification exists at the specified position
     *
     * @note This method only removes the modification; it does not affect the
     *       underlying amino acid sequence.
     */
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
        OPENMS_LOG_DEBUG << "Removing modification: " << it->second.getName() << " at position: " << position << std::endl;
        modifications_.erase(it);
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid position", "No modification found at the specified position: " + std::to_string(position));
      }
    }

    /**
     * @brief Add modification to specified position or range
     *
     * Adds a modification to the internal modifications map at the specified position(s).
     * This method allows for programmatic construction of peptidoforms by adding
     * modifications after sequence initialization.
     *
     * Position mapping:
     * - 0: N-terminal modification
     * - 1 to sequence_length: Amino acid positions (1-based indexing)
     * - sequence_length + 1: C-terminal modification
     *
     * Modification types supported:
     * - Named modifications: Specify non-empty mod_id
     * - Mass shift modifications: Specify non-zero mass_shift
     * - Combined: Both mod_id and mass_shift can be provided
     * - Range modifications: start_pos != end_pos
     *
     * @param start_pos Starting position for the modification (0-based for terminals)
     * @param end_pos Ending position for range modifications (inclusive)
     * @param mod_id Modification identifier or name (empty string for mass-only modifications)
     * @param mass_shift Mass difference in Daltons (0.0 for name-only modifications)
     *
     * @throws Exception::OutOfRange If start_pos exceeds sequence length + 1
     *
     * @note For range modifications (start_pos != end_pos), the same modification
     *       attributes are applied to all positions within the range inclusive.
     *       Single position modifications have start_pos == end_pos.
     */
    void ProForma::addModification(size_t start_pos, size_t end_pos, const String& mod_id, double mass_shift)
    {
      OPENMS_LOG_DEBUG << "addModification called: mod_id=" << mod_id << ", start_pos=" << start_pos
                       << ", end_pos=" << end_pos << ", sequence_size=" << sequence_.size() << std::endl;
      
      if (start_pos > sequence_.size() + 1)
      {
        throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      ResidueModification mod;
      
      // Set modification properties
      if (!mod_id.empty())
      {
        mod.setName(mod_id);
      }
      
      if (mass_shift != 0.0)
      {
        mod.setDiffMonoMass(mass_shift);
      }
      
      // Set position range
      if (start_pos != end_pos)
      {
        // Range modification
        mod.setPositionRange(start_pos, end_pos);
        // Apply the modification to all positions in the range
        for (size_t i = start_pos; i <= end_pos; ++i)
        {
          modifications_[i] = mod;
        }
      }
      else
      {
        // Single position modification
        mod.setPositionRange(start_pos, start_pos);
        modifications_[start_pos] = mod;
      }
      
      OPENMS_LOG_DEBUG << "Added modification: " << mod_id << " with mass shift: " << mass_shift
                       << " at position(s): " << start_pos << " to " << end_pos << " (sequence size: " << sequence_.size() << ")" << std::endl;
    }
    }