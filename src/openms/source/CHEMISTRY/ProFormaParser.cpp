// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaParser.h>
#include <OpenMS/CHEMISTRY/ProFormaWriter.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <cstdlib>
#include <set>
#include <sstream>

namespace OpenMS
{

  // ============================================================================
  // Public static methods
  // ============================================================================

  Peptidoform ProFormaParser::parse(const String& input)
  {
    ProFormaParser parser(input);
    Peptidoform pf = parser.parsePeptidoform_();

    // Verify we consumed all input
    if (!parser.isAtEnd_())
    {
      parser.error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
                    "Unexpected characters after peptidoform");
    }

    return pf;
  }

  PeptidoformIon ProFormaParser::parseIon(const String& input)
  {
    ProFormaParser parser(input);
    return parser.parsePeptidoformIon_();
  }

  // ============================================================================
  // Constructor
  // ============================================================================

  ProFormaParser::ProFormaParser(std::string_view input)
    : tokenizer_(input),
      input_(input),
      current_token_{ProFormaTokenizer::TokenType::END, {}, 0},
      has_current_(false)
  {
  }

  // ============================================================================
  // High-level parsing methods
  // ============================================================================

  PeptidoformIon ProFormaParser::parsePeptidoformIon_()
  {
    PeptidoformIon ion;

    // Parse first chain (with optional per-chain charge for chimeric spectra)
    // First chain doesn't know if it's chimeric yet
    ion.chains.push_back(parsePeptidoformWithCharge_(false));

    // Parse additional chains separated by // (cross-linked) or + (chimeric)
    while (!isAtEnd_())
    {
      if (match_(ProFormaTokenizer::TokenType::SLASH))
      {
        if (!match_(ProFormaTokenizer::TokenType::SLASH))
        {
          // Single slash at top level - this is ion-level charge state
          ion.charge = parseChargeState_();
          break;
        }
        // Double slash - parse another cross-linked chain
        ion.chains.push_back(parsePeptidoformWithCharge_(false));
      }
      else if (match_(ProFormaTokenizer::TokenType::PLUS))
      {
        // Plus sign - parse another chimeric chain
        ion.is_chimeric = true;
        // Pass true to indicate chimeric context, so trailing charge is per-chain
        ion.chains.push_back(parsePeptidoformWithCharge_(true));
      }
      else
      {
        break;
      }
    }

    // Verify we consumed all input
    if (!isAtEnd_())
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
             "Unexpected characters after peptidoform ion");
    }

    return ion;
  }

  Peptidoform ProFormaParser::parsePeptidoformWithCharge_(bool is_chimeric_context)
  {
    Peptidoform pf = parsePeptidoform_();

    // Check for per-chain charge state (e.g., PEPTIDE/2 in chimeric context)
    // This is tricky: / can be per-chain charge OR ion-level charge OR chain separator (//)
    // We need to look ahead to distinguish:
    // - /2+ means per-chain charge 2, then chimeric separator
    // - /2 at end in chimeric context means per-chain charge
    // - /2 at end in non-chimeric context means ion-level charge (handled by caller)
    // - // means chain separator (handled by caller)
    if (check_(ProFormaTokenizer::TokenType::SLASH))
    {
      ProFormaTokenizer lookahead = createLookahead_();
      lookahead.next(); // consume /

      ProFormaTokenizer::Token next = lookahead.peek();
      if (next.type == ProFormaTokenizer::TokenType::SLASH)
      {
        // Double slash - chain separator, don't consume
        return pf;
      }

      // Check if this looks like a charge
      if (next.type == ProFormaTokenizer::TokenType::PLUS ||
          next.type == ProFormaTokenizer::TokenType::MINUS ||
          next.type == ProFormaTokenizer::TokenType::NUMBER ||
          next.type == ProFormaTokenizer::TokenType::LBRACKET)
      {
        // Could be charge state - look ahead to see what follows
        if (next.type == ProFormaTokenizer::TokenType::LBRACKET)
        {
          // Adduct list - scan past it
          int depth = 1;
          lookahead.next(); // consume [
          while (lookahead.hasMore() && depth > 0)
          {
            ProFormaTokenizer::Token tok = lookahead.next();
            if (tok.type == ProFormaTokenizer::TokenType::LBRACKET) depth++;
            else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET) depth--;
          }
        }
        else
        {
          // Simple charge - skip +/- and number
          if (next.type == ProFormaTokenizer::TokenType::PLUS ||
              next.type == ProFormaTokenizer::TokenType::MINUS)
          {
            lookahead.next();
          }
          if (lookahead.peek().type == ProFormaTokenizer::TokenType::NUMBER)
          {
            lookahead.next();
          }
        }

        // Check what follows the charge
        ProFormaTokenizer::Token after_charge = lookahead.peek();
        bool followed_by_chimeric = (after_charge.type == ProFormaTokenizer::TokenType::PLUS);
        bool at_end = (after_charge.type == ProFormaTokenizer::TokenType::END);

        // Parse as per-chain charge if:
        // - Followed by + (chimeric separator)
        // - At end AND we're in a chimeric context (this is the last chain in chimeric)
        if (followed_by_chimeric || (at_end && is_chimeric_context))
        {
          advance_(); // consume / in main parser
          pf.charge = parseChargeState_();
        }
        // Otherwise, leave it for the caller to handle as ion-level charge
      }
    }

    return pf;
  }

  Peptidoform ProFormaParser::parsePeptidoform_()
  {
    Peptidoform pf;

    // Check for gene/protein prefix: (>name) or (>>name) or (>>>name)
    // Multiple consecutive prefixes are allowed: (>>>Org)(>>Sub)(>Name)
    // This is a v2.1 extension for annotating the source gene/protein
    std::string combined_name;
    while (check_(ProFormaTokenizer::TokenType::LPAREN))
    {
      ProFormaTokenizer lookahead = createLookahead_();
      lookahead.next(); // consume (
      ProFormaTokenizer::Token next = lookahead.peek();

      // Check for (> or (>> or (>>> which indicates gene/protein prefix
      // Multiple > characters indicate hierarchy levels
      if (next.type != ProFormaTokenizer::TokenType::RANGLE)
      {
        break; // Not a gene prefix, might be ambiguous region or modified range
      }

      advance_(); // consume (

      // Consume all > characters (hierarchy level indicator)
      while (check_(ProFormaTokenizer::TokenType::RANGLE))
      {
        advance_(); // consume >
      }

      // Parse the name (can be multiple tokens, and may contain nested parentheses)
      std::string name;
      int paren_depth = 0;
      while (!isAtEnd_())
      {
        ProFormaTokenizer::Token tok = current_();

        if (tok.type == ProFormaTokenizer::TokenType::LPAREN)
        {
          paren_depth++;
          name += std::string(tok.text);
          advance_();
        }
        else if (tok.type == ProFormaTokenizer::TokenType::RPAREN)
        {
          if (paren_depth > 0)
          {
            paren_depth--;
            name += std::string(tok.text);
            advance_();
          }
          else
          {
            // This is the closing paren of the gene prefix
            break;
          }
        }
        else
        {
          name += std::string(tok.text);
          advance_();
        }
      }
      expect_(ProFormaTokenizer::TokenType::RPAREN, "')' to close gene/protein prefix");

      // Concatenate with previous prefixes
      if (!combined_name.empty())
      {
        combined_name += " / ";
      }
      combined_name += name;
    }

    if (!combined_name.empty())
    {
      pf.name = String(combined_name);
    }

    // Parse global modifications: <...>
    // There can be multiple consecutive blocks like <13C><15N>
    while (check_(ProFormaTokenizer::TokenType::LANGLE))
    {
      auto mods = parseGlobalMods_();
      for (auto& m : mods)
      {
        pf.global_mods.push_back(std::move(m));
      }
    }

    // Parse unlocalised modifications: [mod]? or [mod]^2?
    pf.unlocalised_mods = parseUnlocalisedMods_();

    // Parse labile modifications: {mod}
    pf.labile_mods = parseLabileModifications_();

    // Check for N-terminal modification: [mod]- or [mod1][mod2]-
    // Use lookahead to determine if brackets are N-terminal or part of sequence
    if (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      bool is_nterm = hasNTerminalModPattern_();

      if (is_nterm)
      {
        pf.n_term_mods = parseTerminalMods_();
        expect_(ProFormaTokenizer::TokenType::MINUS, "'-' after N-terminal modification");
      }
    }

    // Parse the sequence
    pf.sequence = parseSequence_();

    if (pf.sequence.empty())
    {
      error_(ProFormaErrorCode::EMPTY_SEQUENCE, "Empty sequence");
    }

    // Check for C-terminal modification: -[mod]
    if (match_(ProFormaTokenizer::TokenType::MINUS))
    {
      pf.c_term_mods = parseTerminalMods_();
      if (pf.c_term_mods.empty())
      {
        error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
               "Expected modification after '-' for C-terminal modification");
      }
    }

    return pf;
  }

  ProFormaTokenizer ProFormaParser::createLookahead_() const
  {
    // Create a tokenizer positioned at the current logical position
    // If we have a cached token, use its position; otherwise use the tokenizer's position
    size_t start_pos = has_current_ ? current_token_.position : tokenizer_.position();
    return ProFormaTokenizer(input_, start_pos);
  }

  bool ProFormaParser::hasNTerminalModPattern_()
  {
    // Create a lookahead tokenizer at the current position
    ProFormaTokenizer lookahead = createLookahead_();

    // Check for [...]- or [...][...]- pattern
    if (lookahead.peek().type != ProFormaTokenizer::TokenType::LBRACKET)
    {
      return false;
    }

    // Scan through all consecutive bracket groups
    while (lookahead.peek().type == ProFormaTokenizer::TokenType::LBRACKET)
    {
      int depth = 1;
      lookahead.next(); // consume [

      while (lookahead.hasMore() && depth > 0)
      {
        ProFormaTokenizer::Token tok = lookahead.next();
        if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
        {
          depth++;
        }
        else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET)
        {
          depth--;
        }
      }

      if (depth != 0)
      {
        // Unclosed bracket - not a valid N-term pattern
        return false;
      }
    }

    // Check if followed by MINUS
    return lookahead.peek().type == ProFormaTokenizer::TokenType::MINUS;
  }

  std::vector<GlobalModEntry> ProFormaParser::parseGlobalMods_()
  {
    std::vector<GlobalModEntry> entries;

    if (!match_(ProFormaTokenizer::TokenType::LANGLE))
    {
      return entries;
    }

    // Parse first entry
    entries.push_back(parseGlobalModEntry_());

    // Parse additional entries separated by commas
    while (match_(ProFormaTokenizer::TokenType::COMMA))
    {
      entries.push_back(parseGlobalModEntry_());
    }

    expect_(ProFormaTokenizer::TokenType::RANGLE, "'>' to close global modifications");

    return entries;
  }

  GlobalModEntry ProFormaParser::parseGlobalModEntry_()
  {
    // Check if this is an isotope replacement (just an identifier like 13C, 15N, D)
    // or a global modification ([mod]@locations)

    if (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      return parseGlobalModification_();
    }
    else
    {
      return parseIsotopeReplacement_();
    }
  }

  IsotopeReplacement ProFormaParser::parseIsotopeReplacement_()
  {
    IsotopeReplacement ir;

    // Isotope can be: 13C, 15N, D, 2H, etc.
    // It may start with a number token followed by identifier
    std::string isotope_str;

    if (check_(ProFormaTokenizer::TokenType::NUMBER))
    {
      ProFormaTokenizer::Token num = advance_();
      isotope_str = std::string(num.text);
    }

    if (check_(ProFormaTokenizer::TokenType::IDENTIFIER))
    {
      ProFormaTokenizer::Token id = advance_();
      isotope_str += std::string(id.text);
    }

    if (isotope_str.empty())
    {
      error_(ProFormaErrorCode::INVALID_FORMULA, "Expected isotope specification");
    }

    ir.isotope = isotope_str;
    return ir;
  }

  GlobalModification ProFormaParser::parseGlobalModification_()
  {
    GlobalModification gm;

    // Parse [modification]
    expect_(ProFormaTokenizer::TokenType::LBRACKET, "'['");
    gm.modification = parseModification_();
    expect_(ProFormaTokenizer::TokenType::RBRACKET, "']'");

    // Labels are not allowed on global modifications
    for (const auto& alt : gm.modification.alternatives)
    {
      if (alt.second.has_value())
      {
        error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
               "Labels are not allowed on global modifications");
      }
    }

    // Parse @locations
    expect_(ProFormaTokenizer::TokenType::AT, "'@' for global modification locations");

    // Parse locations: K, N-term, C-term, C-term:K, etc.
    // Locations are comma-separated within the global mod context
    std::string location;

    // First location
    if (check_(ProFormaTokenizer::TokenType::IDENTIFIER))
    {
      ProFormaTokenizer::Token id = advance_();
      location = std::string(id.text);

      // Check for -term suffix
      if (match_(ProFormaTokenizer::TokenType::MINUS))
      {
        ProFormaTokenizer::Token term = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "term identifier");
        location += "-" + std::string(term.text);
      }

      // Check for :X suffix (e.g., C-term:K)
      if (match_(ProFormaTokenizer::TokenType::COLON))
      {
        ProFormaTokenizer::Token suffix = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "location suffix");
        location += ":" + std::string(suffix.text);
      }

      gm.locations.push_back(location);
    }

    // Additional locations separated by commas (but within <>)
    while (check_(ProFormaTokenizer::TokenType::COMMA))
    {
      // Peek ahead to see if this comma separates locations or global mod entries
      // If the token after comma is [ or a number (isotope), it's a new global mod entry
      // Note: tokenizer_ is already past the comma (current_() consumed it), so we just peek
      ProFormaTokenizer lookahead = tokenizer_;
      ProFormaTokenizer::Token after_comma = lookahead.peek();

      if (after_comma.type == ProFormaTokenizer::TokenType::LBRACKET ||
          after_comma.type == ProFormaTokenizer::TokenType::NUMBER)
      {
        // New global mod entry, stop parsing locations
        break;
      }

      advance_(); // consume comma
      location.clear();

      if (check_(ProFormaTokenizer::TokenType::IDENTIFIER))
      {
        ProFormaTokenizer::Token id = advance_();
        location = std::string(id.text);

        if (match_(ProFormaTokenizer::TokenType::MINUS))
        {
          ProFormaTokenizer::Token term = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "term identifier");
          location += "-" + std::string(term.text);
        }

        if (match_(ProFormaTokenizer::TokenType::COLON))
        {
          ProFormaTokenizer::Token suffix = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "location suffix");
          location += ":" + std::string(suffix.text);
        }

        gm.locations.push_back(location);
      }
    }

    return gm;
  }

  std::vector<UnlocalisedMod> ProFormaParser::parseUnlocalisedMods_()
  {
    std::vector<UnlocalisedMod> mods;

    // Unlocalised mods look like: [mod]? or [mod]^2? or [mod1][mod2]?
    // Multiple consecutive brackets can share a single ?
    // They appear before the sequence and end with ?

    while (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      // Look ahead to see if this group ends with ?
      // A group can be [mod]? or [mod1][mod2]? etc.
      ProFormaTokenizer lookahead = createLookahead_();

      // Count how many consecutive [...] groups there are
      int group_count = 0;
      bool found_question = false;

      while (true)
      {
        ProFormaTokenizer::Token tok = lookahead.peek();
        if (tok.type != ProFormaTokenizer::TokenType::LBRACKET)
        {
          break;
        }

        // Scan past this bracket content
        int depth = 0;
        tok = lookahead.next(); // [
        depth = 1;
        group_count++;

        while (lookahead.hasMore() && depth > 0)
        {
          tok = lookahead.next();
          if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
          {
            depth++;
          }
          else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET)
          {
            depth--;
          }
        }

        // After ], check what's next
        tok = lookahead.peek();

        // Check for ^N
        if (tok.type == ProFormaTokenizer::TokenType::CARET)
        {
          lookahead.next(); // consume ^
          tok = lookahead.peek();
          if (tok.type == ProFormaTokenizer::TokenType::NUMBER)
          {
            lookahead.next(); // consume number
            tok = lookahead.peek();
          }
        }

        // Check for ?
        if (tok.type == ProFormaTokenizer::TokenType::QUESTION)
        {
          found_question = true;
          break;
        }

        // If next is another [, continue the group
        // Otherwise, this is not an unlocalised mod group
        if (tok.type != ProFormaTokenizer::TokenType::LBRACKET)
        {
          break;
        }
      }

      if (!found_question)
      {
        // Not an unlocalised mod
        break;
      }

      // Now actually parse the group
      // Each [mod] or [mod]^N is a separate UnlocalisedMod entry
      while (check_(ProFormaTokenizer::TokenType::LBRACKET))
      {
        UnlocalisedMod um;

        expect_(ProFormaTokenizer::TokenType::LBRACKET, "'['");
        um.modifications.push_back(parseModification_());
        expect_(ProFormaTokenizer::TokenType::RBRACKET, "']'");

        // Check for occurrence specifier ^N
        if (match_(ProFormaTokenizer::TokenType::CARET))
        {
          ProFormaTokenizer::Token num = expect_(ProFormaTokenizer::TokenType::NUMBER, "occurrence count");
          try
          {
            um.occurrence = std::stoi(std::string(num.text));
          }
          catch (const std::exception&)
          {
            errorAt_(ProFormaErrorCode::INVALID_MASS_VALUE, num.position, "Invalid occurrence count");
          }
        }

        mods.push_back(std::move(um));

        // Continue if next is another [ (more mods in this group before ?)
        // Break if next is ? (end of unlocalised group)
        if (check_(ProFormaTokenizer::TokenType::QUESTION))
        {
          break;
        }
      }

      expect_(ProFormaTokenizer::TokenType::QUESTION, "'?' for unlocalised modification");
    }

    return mods;
  }

  std::vector<LabileModification> ProFormaParser::parseLabileModifications_()
  {
    std::vector<LabileModification> mods;

    // Labile mods look like: {mod}
    while (match_(ProFormaTokenizer::TokenType::LBRACE))
    {
      LabileModification lm;
      lm.modification = parseModification_();

      // Labels are not allowed on labile modifications
      for (const auto& alt : lm.modification.alternatives)
      {
        if (alt.second.has_value())
        {
          error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
                 "Labels are not allowed on labile modifications");
        }
      }

      expect_(ProFormaTokenizer::TokenType::RBRACE, "'}'");
      mods.push_back(std::move(lm));
    }

    return mods;
  }

  std::vector<SequenceSection> ProFormaParser::parseSequence_()
  {
    std::vector<SequenceSection> sections;

    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      // Check for end of sequence markers
      if (tok.type == ProFormaTokenizer::TokenType::MINUS ||
          tok.type == ProFormaTokenizer::TokenType::SLASH ||
          tok.type == ProFormaTokenizer::TokenType::END)
      {
        break;
      }

      // Check for ambiguous region: (?XY)
      if (tok.type == ProFormaTokenizer::TokenType::LPAREN)
      {
        // Look ahead to see if this is (?...) or (...)
        ProFormaTokenizer lookahead = createLookahead_();
        lookahead.next(); // consume (
        ProFormaTokenizer::Token after_paren = lookahead.peek();

        if (after_paren.type == ProFormaTokenizer::TokenType::QUESTION)
        {
          sections.push_back(parseAmbiguousRegion_());
        }
        else
        {
          sections.push_back(parseModifiedRange_());
        }
        continue;
      }

      // Regular amino acid
      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        // Each letter is a separate amino acid
        std::string_view text = tok.text;
        advance_(); // consume identifier

        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c))
          {
            errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position,
                     "Invalid amino acid character");
          }

          SequenceElement elem;
          elem.amino_acid = c;

          // Only the last character can have modifications attached
          // E.g., "EM[mod]" means M has the modification, not E
          if (i == text.size() - 1 && check_(ProFormaTokenizer::TokenType::LBRACKET))
          {
            elem.modifications = parseModificationList_();
          }

          sections.push_back(elem);
        }
        continue;
      }

      // If we get here, it's an unexpected token
      break;
    }

    return sections;
  }

  SequenceElement ProFormaParser::parseSequenceElement_()
  {
    SequenceElement elem;

    ProFormaTokenizer::Token tok = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "amino acid");

    if (tok.text.size() != 1 || !isAminoAcid_(tok.text[0]))
    {
      errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position,
               "Expected single amino acid letter");
    }

    elem.amino_acid = tok.text[0];

    // Check for modifications
    if (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      elem.modifications = parseModificationList_();
    }

    return elem;
  }

  AmbiguousRegion ProFormaParser::parseAmbiguousRegion_()
  {
    AmbiguousRegion region;

    expect_(ProFormaTokenizer::TokenType::LPAREN, "'('");
    expect_(ProFormaTokenizer::TokenType::QUESTION, "'?' for ambiguous region");

    // Parse amino acids until )
    // The tokenizer combines consecutive letters into one identifier,
    // so we need to iterate through multi-letter identifiers (e.g., "DQ" -> D, Q)
    while (!check_(ProFormaTokenizer::TokenType::RPAREN) && !isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();

        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c))
          {
            errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position,
                     "Invalid amino acid in ambiguous region");
          }

          SequenceElement elem;
          elem.amino_acid = c;

          // Only the last character can have modifications attached
          if (i == text.size() - 1 && check_(ProFormaTokenizer::TokenType::LBRACKET))
          {
            elem.modifications = parseModificationList_();
          }

          region.elements.push_back(elem);
        }
      }
      else
      {
        break;
      }
    }

    expect_(ProFormaTokenizer::TokenType::RPAREN, "')'");

    return region;
  }

  ModifiedRange ProFormaParser::parseModifiedRange_()
  {
    ModifiedRange range;

    expect_(ProFormaTokenizer::TokenType::LPAREN, "'('");

    // Parse sequence elements until )
    // Elements can have modifications: (EOC[Carbamidomethyl]FORMS)
    while (!check_(ProFormaTokenizer::TokenType::RPAREN) && !isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();

        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c))
          {
            errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position + i,
                     "Invalid amino acid in range");
          }

          SequenceElement elem;
          elem.amino_acid = c;

          // Check if this is the last character and there's a modification following
          if (i == text.size() - 1 && check_(ProFormaTokenizer::TokenType::LBRACKET))
          {
            elem.modifications = parseModificationList_();
          }

          range.elements.push_back(elem);
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
      {
        // Modification on the last element (shouldn't happen if IDENTIFIER handled it)
        if (!range.elements.empty())
        {
          auto mods = parseModificationList_();
          for (auto& m : mods)
          {
            range.elements.back().modifications.push_back(std::move(m));
          }
        }
        else
        {
          break;
        }
      }
      else
      {
        break;
      }
    }

    expect_(ProFormaTokenizer::TokenType::RPAREN, "')'");

    // Empty ranges are not allowed
    if (range.elements.empty())
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Empty range is not allowed");
    }

    // Parse modifications for the range (applied to entire range)
    if (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      range.modifications = parseModificationList_();
    }

    return range;
  }

  std::vector<Modification> ProFormaParser::parseTerminalMods_()
  {
    std::vector<Modification> mods;

    while (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      expect_(ProFormaTokenizer::TokenType::LBRACKET, "'['");
      mods.push_back(parseModification_());
      expect_(ProFormaTokenizer::TokenType::RBRACKET, "']'");
    }

    return mods;
  }

  // ============================================================================
  // Modification parsing
  // ============================================================================

  std::vector<Modification> ProFormaParser::parseModificationList_()
  {
    std::vector<Modification> mods;

    while (match_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      mods.push_back(parseModification_());
      expect_(ProFormaTokenizer::TokenType::RBRACKET, "']'");
    }

    return mods;
  }

  Modification ProFormaParser::parseModification_()
  {
    Modification mod;

    // Parse first alternative
    auto [tag, label] = parseModificationTagWithLabel_();
    mod.alternatives.emplace_back(std::move(tag), std::move(label));

    // Parse additional alternatives separated by |
    while (match_(ProFormaTokenizer::TokenType::PIPE))
    {
      auto [alt_tag, alt_label] = parseModificationTagWithLabel_();
      mod.alternatives.emplace_back(std::move(alt_tag), std::move(alt_label));
    }

    return mod;
  }

  std::pair<ModificationTag, std::optional<Label>> ProFormaParser::parseModificationTagWithLabel_()
  {
    // ProForma allows label-only modifications like [#XL1] that reference
    // earlier defined cross-links. Check for this case first.
    if (check_(ProFormaTokenizer::TokenType::HASH))
    {
      // Label-only modification - no tag, just a reference
      Label label = parseLabel_();
      // Use empty InfoTag as placeholder for the tag
      InfoTag empty_tag;
      empty_tag.text = "";
      return {std::move(empty_tag), std::move(label)};
    }

    ModificationTag tag = parseModificationTag_();
    std::optional<Label> label;

    // Check for label: #identifier
    if (check_(ProFormaTokenizer::TokenType::HASH))
    {
      label = parseLabel_();
    }

    return {std::move(tag), std::move(label)};
  }

  ModificationTag ProFormaParser::parseModificationTag_()
  {
    ProFormaTokenizer::Token tok = current_();

    // Mass delta: starts with + or - or number
    if (tok.type == ProFormaTokenizer::TokenType::PLUS ||
        tok.type == ProFormaTokenizer::TokenType::MINUS ||
        tok.type == ProFormaTokenizer::TokenType::NUMBER)
    {
      return parseMassDelta_();
    }

    // Must be an identifier-based tag
    if (tok.type != ProFormaTokenizer::TokenType::IDENTIFIER)
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Expected modification");
    }

    std::string_view id = tok.text;

    // Check for specific prefixes
    if (id == "Formula")
    {
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after Formula");
      return parseFormulaTag_();
    }
    else if (id == "Glycan")
    {
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after Glycan");
      return parseGlycanComposition_();
    }
    else if (id == "INFO" || id == "info" || id == "Info")
    {
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after INFO");
      return parseInfoTag_();
    }
    else if (id == "Position" || id == "position" || id == "POSITION")
    {
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after Position");
      return parsePositionConstraint_();
    }
    else if (id == "Cation")
    {
      // Cation tag like Cation:Mg[II]
      // Parse as a NamedMod but include any bracketed charge notation
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after Cation");

      NamedMod nm;
      std::string cation_str;

      // Parse element name
      if (check_(ProFormaTokenizer::TokenType::IDENTIFIER))
      {
        ProFormaTokenizer::Token elem = advance_();
        cation_str = std::string(elem.text);
      }

      // Parse optional charge notation like [II], [III]
      if (match_(ProFormaTokenizer::TokenType::LBRACKET))
      {
        cation_str += "[";
        while (!check_(ProFormaTokenizer::TokenType::RBRACKET) && !isAtEnd_())
        {
          ProFormaTokenizer::Token tok = advance_();
          cation_str += std::string(tok.text);
        }
        expect_(ProFormaTokenizer::TokenType::RBRACKET, "']' for cation charge");
        cation_str += "]";
      }

      nm.name = "Cation:" + cation_str;
      return nm;
    }
    else if (id == "Obs")
    {
      // Mass delta with Obs: prefix
      advance_();
      expect_(ProFormaTokenizer::TokenType::COLON, "':' after Obs");
      MassDelta md = parseMassDelta_();
      md.source = MassDelta::Source::OBS;
      return md;
    }
    else if (id == "UNIMOD" || id == "MOD" || id == "RESID" || id == "XLMOD" || id == "GNO" ||
             id == "unimod" || id == "mod" || id == "resid" || id == "xlmod" || id == "gno" ||
             id == "xlink")  // xlink is an alias for xlmod
    {
      return parseCvAccession_();
    }
    else if (id.size() == 1)
    {
      // Could be CV hint prefix: U:name, M:name, etc.
      char prefix = id[0];
      if (prefix == 'U' || prefix == 'M' || prefix == 'R' || prefix == 'X' || prefix == 'G')
      {
        // Look ahead for colon
        ProFormaTokenizer lookahead = tokenizer_;
        ProFormaTokenizer::Token next_tok = lookahead.peek();
        if (next_tok.type == ProFormaTokenizer::TokenType::COLON)
        {
          advance_(); // consume prefix
          advance_(); // consume colon
          tok = current_();

          // Check if followed by number (mass delta) or identifier (named mod)
          if (tok.type == ProFormaTokenizer::TokenType::PLUS ||
              tok.type == ProFormaTokenizer::TokenType::MINUS ||
              tok.type == ProFormaTokenizer::TokenType::NUMBER)
          {
            MassDelta md = parseMassDelta_();
            switch (prefix)
            {
              case 'U': md.source = MassDelta::Source::U; break;
              case 'M': md.source = MassDelta::Source::M; break;
              case 'R': md.source = MassDelta::Source::R; break;
              case 'X': md.source = MassDelta::Source::X; break;
              case 'G': md.source = MassDelta::Source::G; break;
            }
            return md;
          }
          else
          {
            // Named mod with CV hint
            return parseNamedMod_(prefix);
          }
        }
      }
    }

    // Default: named modification
    return parseNamedMod_();
  }

  NamedMod ProFormaParser::parseNamedMod_()
  {
    NamedMod nm;
    nm.cv_hint = std::nullopt;

    // Named mods can contain letters, numbers, and certain characters
    // E.g., TMT6plex, iTRAQ4plex, half cystine, L-methionine sulfoxide
    // Also: L-cystine (cross-link), Gln->pyro-Glu, dss[138]
    // Tokenizer splits at letter/digit boundaries, so we concatenate tokens
    std::string name;
    int paren_depth = 0;   // Track parentheses depth for names like "(cross-link)"
    int bracket_depth = 0; // Track bracket depth for names like "dss[138]"

    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        name += std::string(tok.text);
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::NUMBER)
      {
        name += std::string(tok.text);
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::MINUS)
      {
        // Hyphens can be part of names like "L-methionine" or "a-type-ion"
        // Also check for arrow: -> (minus followed by RANGLE)
        ProFormaTokenizer lookahead = tokenizer_;
        ProFormaTokenizer::Token next = lookahead.peek();
        if (next.type == ProFormaTokenizer::TokenType::RANGLE)
        {
          // This is an arrow: -> (e.g., Gln->pyro-Glu)
          name += "->";
          advance_(); // consume minus
          advance_(); // consume >
        }
        else if (next.type == ProFormaTokenizer::TokenType::IDENTIFIER ||
                 next.type == ProFormaTokenizer::TokenType::NUMBER ||
                 next.type == ProFormaTokenizer::TokenType::LPAREN ||
                 next.type == ProFormaTokenizer::TokenType::LBRACKET)
        {
          name += "-";
          advance_();
        }
        else if (paren_depth > 0 || bracket_depth > 0)
        {
          // Inside parentheses or brackets, allow hyphen
          name += "-";
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::LPAREN)
      {
        // Parentheses can be part of names like "L-cystine (cross-link)"
        // Only if we already have some name content (not at start)
        if (!name.empty() || paren_depth > 0 || bracket_depth > 0)
        {
          name += "(";
          paren_depth++;
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::RPAREN)
      {
        if (paren_depth > 0)
        {
          name += ")";
          paren_depth--;
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
      {
        // Brackets can be part of names like "dss[138]" (XLMOD accession formats)
        // Only if we already have some name content (not at start)
        if (!name.empty() || paren_depth > 0 || bracket_depth > 0)
        {
          name += "[";
          bracket_depth++;
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET)
      {
        if (bracket_depth > 0)
        {
          name += "]";
          bracket_depth--;
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::RANGLE)
      {
        // Check if this is part of an arrow: -> which is handled above with MINUS
        // Standalone > after name ends the name
        break;
      }
      else
      {
        break;
      }
    }

    if (name.empty())
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Expected modification name");
    }

    nm.name = String(name);
    return nm;
  }

  // Helper overload for when we already know the CV hint
  NamedMod ProFormaParser::parseNamedMod_(char cv_hint)
  {
    // Parse the name using the same logic as the no-hint version
    NamedMod nm = parseNamedMod_();

    // Set the CV hint
    switch (cv_hint)
    {
      case 'U': nm.cv_hint = CvDatabase::UNIMOD; break;
      case 'M': nm.cv_hint = CvDatabase::MOD; break;
      case 'R': nm.cv_hint = CvDatabase::RESID; break;
      case 'X': nm.cv_hint = CvDatabase::XLMOD; break;
      case 'G': nm.cv_hint = CvDatabase::GNO; break;
    }

    return nm;
  }

  CvAccession ProFormaParser::parseCvAccession_()
  {
    CvAccession cv;

    ProFormaTokenizer::Token prefix_tok = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "CV prefix");
    std::string_view prefix = prefix_tok.text;

    // Handle uppercase and lowercase CV prefixes
    if (prefix == "UNIMOD" || prefix == "unimod")
    {
      cv.database = CvDatabase::UNIMOD;
    }
    else if (prefix == "MOD" || prefix == "mod")
    {
      cv.database = CvDatabase::MOD;
    }
    else if (prefix == "RESID" || prefix == "resid")
    {
      cv.database = CvDatabase::RESID;
    }
    else if (prefix == "XLMOD" || prefix == "xlmod" || prefix == "xlink")
    {
      cv.database = CvDatabase::XLMOD;
    }
    else if (prefix == "GNO" || prefix == "gno")
    {
      cv.database = CvDatabase::GNO;
    }
    else
    {
      errorAt_(ProFormaErrorCode::INVALID_CV_PREFIX, prefix_tok.position,
               "Invalid CV prefix");
    }

    expect_(ProFormaTokenizer::TokenType::COLON, "':' after CV prefix");

    // GNO, RESID, and XLMOD accessions can be alphanumeric (e.g., G59626AS, AA0581, dss[138])
    // UNIMOD and MOD accessions are numeric
    if (cv.database == CvDatabase::GNO || cv.database == CvDatabase::RESID || cv.database == CvDatabase::XLMOD)
    {
      // Parse alphanumeric accession (letters, numbers, and brackets for XLMOD like dss[138])
      std::string accession;
      int bracket_depth = 0;
      while (!isAtEnd_())
      {
        ProFormaTokenizer::Token tok = current_();
        if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER ||
            tok.type == ProFormaTokenizer::TokenType::NUMBER)
        {
          accession += std::string(tok.text);
          advance_();
        }
        else if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
        {
          // Brackets in accessions like dss[138]
          accession += "[";
          bracket_depth++;
          advance_();
        }
        else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET)
        {
          if (bracket_depth > 0)
          {
            accession += "]";
            bracket_depth--;
            advance_();
          }
          else
          {
            break; // End of accession
          }
        }
        else
        {
          break;
        }
      }
      if (accession.empty())
      {
        error_(ProFormaErrorCode::INVALID_CV_ACCESSION, "Expected accession");
      }
      cv.accession = accession;
    }
    else
    {
      ProFormaTokenizer::Token num = expect_(ProFormaTokenizer::TokenType::NUMBER, "accession number");
      cv.accession = String(num.text);
    }

    return cv;
  }

  MassDelta ProFormaParser::parseMassDelta_()
  {
    MassDelta md;
    md.source = MassDelta::Source::NONE;

    std::string original;
    ProFormaTokenizer::Token tok = current_();

    // Handle sign
    if (tok.type == ProFormaTokenizer::TokenType::PLUS)
    {
      original = "+";
      advance_();
      tok = current_();
    }
    else if (tok.type == ProFormaTokenizer::TokenType::MINUS)
    {
      original = "-";
      advance_();
      tok = current_();
    }

    // The number token may already include the sign
    if (tok.type == ProFormaTokenizer::TokenType::NUMBER)
    {
      original += std::string(tok.text);
      advance_();
    }
    else
    {
      error_(ProFormaErrorCode::INVALID_MASS_VALUE, "Expected mass value");
    }

    md.original_text = original;

    // Parse the mass value
    try
    {
      md.mass = std::stod(original);
    }
    catch (const std::exception&)
    {
      error_(ProFormaErrorCode::INVALID_MASS_VALUE, "Invalid mass value format");
    }

    return md;
  }

  FormulaTag ProFormaParser::parseFormulaTag_()
  {
    FormulaTag ft;

    // Parse formula string: can contain letters, numbers, parentheses, and
    // square brackets for isotope notation like [13C2][12C-2]H2N
    // Also handles minus sign inside brackets for negative counts
    std::string formula;
    int bracket_depth = 0;  // Track square bracket depth for isotope notation

    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        formula += std::string(tok.text);
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::NUMBER)
      {
        formula += std::string(tok.text);
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::LPAREN)
      {
        formula += "(";
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::RPAREN)
      {
        formula += ")";
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
      {
        // Square bracket for isotope notation like [13C2]
        formula += "[";
        bracket_depth++;
        advance_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::RBRACKET)
      {
        if (bracket_depth > 0)
        {
          // Closing isotope bracket
          formula += "]";
          bracket_depth--;
          advance_();
        }
        else
        {
          // This is the closing bracket of the modification, not part of formula
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::MINUS)
      {
        // Minus can be:
        // 1. Inside isotope bracket for negative counts like [12C-2]
        // 2. Outside brackets for negative element counts like H-1C-1O-2
        // Peek ahead to see if followed by a number
        ProFormaTokenizer lookahead = tokenizer_;
        ProFormaTokenizer::Token next = lookahead.peek();
        if (bracket_depth > 0 || next.type == ProFormaTokenizer::TokenType::NUMBER)
        {
          // Minus is part of the formula
          formula += "-";
          advance_();
        }
        else
        {
          // Minus not part of formula (e.g., C-terminal delimiter)
          break;
        }
      }
      else if (tok.type == ProFormaTokenizer::TokenType::COLON)
      {
        // Check for :z+N charge specifier
        ProFormaTokenizer lookahead = tokenizer_;
        ProFormaTokenizer::Token next = lookahead.peek();
        if (next.type == ProFormaTokenizer::TokenType::IDENTIFIER &&
            next.text == "z")
        {
          advance_(); // consume :
          advance_(); // consume z

          // Parse charge: +N or -N
          int sign = 1;
          if (match_(ProFormaTokenizer::TokenType::PLUS))
          {
            sign = 1;
          }
          else if (match_(ProFormaTokenizer::TokenType::MINUS))
          {
            sign = -1;
          }

          ProFormaTokenizer::Token num = expect_(ProFormaTokenizer::TokenType::NUMBER, "charge value");
          try
          {
            ft.charge = sign * std::stoi(std::string(num.text));
          }
          catch (const std::exception&)
          {
            errorAt_(ProFormaErrorCode::INVALID_CHARGE, num.position, "Invalid charge value");
          }
          break;
        }
        else
        {
          break;
        }
      }
      else
      {
        break;
      }
    }

    if (formula.empty())
    {
      error_(ProFormaErrorCode::INVALID_FORMULA, "Empty formula");
    }

    ft.formula_string = formula;
    return ft;
  }

  GlycanComposition ProFormaParser::parseGlycanComposition_()
  {
    GlycanComposition gc;

    // Parse monosaccharide components: Name1Count1Name2Count2...
    // e.g., HexNAc2Hex3
    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type != ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        break;
      }

      std::string mono_name(tok.text);
      advance_();

      // Check for count
      int count = 1;
      if (check_(ProFormaTokenizer::TokenType::NUMBER))
      {
        ProFormaTokenizer::Token num = advance_();
        try
        {
          count = std::stoi(std::string(num.text));
        }
        catch (const std::exception&)
        {
          errorAt_(ProFormaErrorCode::INVALID_MASS_VALUE, num.position, "Invalid monosaccharide count");
        }
      }

      gc.components.emplace_back(String(mono_name), count);
    }

    if (gc.components.empty())
    {
      error_(ProFormaErrorCode::UNKNOWN_MONOSACCHARIDE, "Empty glycan composition");
    }

    return gc;
  }

  InfoTag ProFormaParser::parseInfoTag_()
  {
    InfoTag it;

    // Collect all remaining text until the next structural delimiter
    std::string text;

    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      // Stop at structural delimiters
      if (tok.type == ProFormaTokenizer::TokenType::RBRACKET ||
          tok.type == ProFormaTokenizer::TokenType::PIPE ||
          tok.type == ProFormaTokenizer::TokenType::HASH ||
          tok.type == ProFormaTokenizer::TokenType::COMMA)
      {
        break;
      }

      text += std::string(tok.text);
      advance_();
    }

    it.text = text;
    return it;
  }

  PositionConstraint ProFormaParser::parsePositionConstraint_()
  {
    PositionConstraint pc;

    // Collect amino acid characters or terminal positions until the next structural delimiter
    while (!isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      // Stop at structural delimiters
      if (tok.type == ProFormaTokenizer::TokenType::RBRACKET ||
          tok.type == ProFormaTokenizer::TokenType::PIPE ||
          tok.type == ProFormaTokenizer::TokenType::HASH ||
          tok.type == ProFormaTokenizer::TokenType::RBRACE)
      {
        break;
      }

      // Check for N-term or C-term (identifier followed by -term)
      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER &&
          (tok.text == "N" || tok.text == "C"))
      {
        // Look ahead for -term pattern
        ProFormaTokenizer lookahead = createLookahead_();
        lookahead.next(); // consume N or C
        ProFormaTokenizer::Token next1 = lookahead.peek();
        if (next1.type == ProFormaTokenizer::TokenType::MINUS)
        {
          lookahead.next(); // consume -
          ProFormaTokenizer::Token next2 = lookahead.peek();
          if (next2.type == ProFormaTokenizer::TokenType::IDENTIFIER && next2.text == "term")
          {
            // This is N-term or C-term
            if (tok.text == "N")
            {
              pc.n_term = true;
            }
            else
            {
              pc.c_term = true;
            }
            advance_(); // consume N or C
            advance_(); // consume -
            advance_(); // consume term

            // Check for comma separator to continue
            if (check_(ProFormaTokenizer::TokenType::COMMA))
            {
              advance_(); // consume comma
            }
            continue;
          }
        }
      }

      // Check for comma (separator between positions)
      if (tok.type == ProFormaTokenizer::TokenType::COMMA)
      {
        advance_();
        continue;
      }

      // Accept identifiers (which contain the amino acid letters)
      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        for (char c : tok.text)
        {
          if (isAminoAcid_(c))
          {
            pc.residues.push_back(c);
          }
          else
          {
            errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position,
                     "Invalid amino acid in position constraint");
          }
        }
        advance_();
      }
      else
      {
        error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
               "Expected amino acid residues or terminal position after Position:");
      }
    }

    if (pc.residues.empty() && !pc.n_term && !pc.c_term)
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
             "Position constraint requires at least one residue or terminal position");
    }

    return pc;
  }

  Label ProFormaParser::parseLabel_()
  {
    Label label;

    expect_(ProFormaTokenizer::TokenType::HASH, "'#'");

    // Labels can be alphanumeric (XL1, g1, s1) or just identifiers (BRANCH)
    // Tokenizer splits at letter/digit boundaries, so we need to combine tokens
    // E.g., "XL1" becomes IDENTIFIER("XL") + NUMBER("1")
    // Also handle pure numeric labels like "#1" for ambiguous modifications

    std::string label_str;

    if (check_(ProFormaTokenizer::TokenType::IDENTIFIER))
    {
      ProFormaTokenizer::Token id = advance_();
      label_str = std::string(id.text);

      // Check if followed by a number (e.g., XL1, g1)
      if (check_(ProFormaTokenizer::TokenType::NUMBER))
      {
        ProFormaTokenizer::Token num = advance_();
        label_str += std::string(num.text);
      }
    }
    else if (check_(ProFormaTokenizer::TokenType::NUMBER))
    {
      // Pure numeric label (e.g., #1, #2)
      ProFormaTokenizer::Token num = advance_();
      label_str = std::string(num.text);
    }
    else
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Expected label identifier");
    }

    label.identifier = String(label_str);

    // Determine label type based on identifier
    if (label.identifier == "BRANCH")
    {
      label.type = Label::Type::BRANCH;
    }
    else if (label.identifier.hasPrefix("XL"))
    {
      label.type = Label::Type::CROSSLINK;
    }
    else
    {
      label.type = Label::Type::AMBIGUOUS;
    }

    // Check for optional score: (0.90)
    if (match_(ProFormaTokenizer::TokenType::LPAREN))
    {
      ProFormaTokenizer::Token score_tok = expect_(ProFormaTokenizer::TokenType::NUMBER, "score value");
      try
      {
        label.score = std::stod(std::string(score_tok.text));
      }
      catch (const std::exception&)
      {
        errorAt_(ProFormaErrorCode::INVALID_MASS_VALUE, score_tok.position, "Invalid score value");
      }
      expect_(ProFormaTokenizer::TokenType::RPAREN, "')'");
    }

    return label;
  }

  // ============================================================================
  // Charge state parsing
  // ============================================================================

  std::optional<ChargeState> ProFormaParser::parseChargeState_()
  {
    // At this point, the first / has been consumed
    // Check if this is a simple charge or adduct list

    ProFormaTokenizer::Token tok = current_();

    if (tok.type == ProFormaTokenizer::TokenType::LBRACKET)
    {
      // Adduct list: [Na:z+1, H:z+1]
      return parseAdductIons_();
    }
    else
    {
      // Simple charge: +2, -1, 2
      int sign = 1;
      if (tok.type == ProFormaTokenizer::TokenType::PLUS)
      {
        sign = 1;
        advance_();
        tok = current_();
      }
      else if (tok.type == ProFormaTokenizer::TokenType::MINUS)
      {
        sign = -1;
        advance_();
        tok = current_();
      }

      if (tok.type == ProFormaTokenizer::TokenType::NUMBER)
      {
        int charge;
        try
        {
          charge = sign * std::stoi(std::string(tok.text));
        }
        catch (const std::exception&)
        {
          errorAt_(ProFormaErrorCode::INVALID_CHARGE, tok.position, "Invalid charge value");
        }
        advance_();
        return charge;
      }
      else
      {
        error_(ProFormaErrorCode::INVALID_CHARGE, "Expected charge value");
      }
    }
  }

  std::vector<AdductIon> ProFormaParser::parseAdductIons_()
  {
    std::vector<AdductIon> adducts;

    expect_(ProFormaTokenizer::TokenType::LBRACKET, "'['");

    // Parse first adduct
    adducts.push_back(parseAdductIon_());

    // Parse additional adducts
    while (match_(ProFormaTokenizer::TokenType::COMMA))
    {
      adducts.push_back(parseAdductIon_());
    }

    expect_(ProFormaTokenizer::TokenType::RBRACKET, "']'");

    return adducts;
  }

  AdductIon ProFormaParser::parseAdductIon_()
  {
    AdductIon adduct;

    // Parse formula (e.g., Na, H, K, Al H-3)
    // Formula can contain spaces and special chars, continues until :z
    std::string formula;
    while (!isAtEnd_())
    {
      // Look ahead for :z pattern which ends the formula
      if (check_(ProFormaTokenizer::TokenType::COLON))
      {
        ProFormaTokenizer lookahead = createLookahead_();
        lookahead.next(); // consume :
        ProFormaTokenizer::Token next = lookahead.peek();
        if (next.type == ProFormaTokenizer::TokenType::IDENTIFIER && next.text == "z")
        {
          break; // End of formula
        }
      }

      // Collect this token as part of formula
      ProFormaTokenizer::Token tok = current_();
      if (tok.type == ProFormaTokenizer::TokenType::RBRACKET ||
          tok.type == ProFormaTokenizer::TokenType::COMMA)
      {
        error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Unexpected token in adduct formula, expected ':z'");
      }
      formula += std::string(tok.text);
      advance_();
    }

    if (formula.empty())
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER, "Expected adduct formula");
    }
    adduct.formula = String(formula);

    // Parse :z+N or :z-N
    expect_(ProFormaTokenizer::TokenType::COLON, "':'");
    expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "'z'"); // Should be 'z'

    int sign = 1;
    if (match_(ProFormaTokenizer::TokenType::PLUS))
    {
      sign = 1;
    }
    else if (match_(ProFormaTokenizer::TokenType::MINUS))
    {
      sign = -1;
    }

    ProFormaTokenizer::Token charge_tok = expect_(ProFormaTokenizer::TokenType::NUMBER, "charge value");
    try
    {
      adduct.charge = sign * std::stoi(std::string(charge_tok.text));
    }
    catch (const std::exception&)
    {
      errorAt_(ProFormaErrorCode::INVALID_CHARGE, charge_tok.position, "Invalid adduct charge value");
    }

    // Check for occurrence ^N
    if (match_(ProFormaTokenizer::TokenType::CARET))
    {
      ProFormaTokenizer::Token occ = expect_(ProFormaTokenizer::TokenType::NUMBER, "occurrence count");
      try
      {
        adduct.occurrence = std::stoi(std::string(occ.text));
      }
      catch (const std::exception&)
      {
        errorAt_(ProFormaErrorCode::INVALID_MASS_VALUE, occ.position, "Invalid adduct occurrence count");
      }
    }

    return adduct;
  }

  // ============================================================================
  // Helper methods
  // ============================================================================

  ProFormaTokenizer::Token ProFormaParser::current_()
  {
    if (!has_current_)
    {
      current_token_ = tokenizer_.next();
      has_current_ = true;
    }
    return current_token_;
  }

  ProFormaTokenizer::Token ProFormaParser::peek_()
  {
    return tokenizer_.peek();
  }

  ProFormaTokenizer::Token ProFormaParser::advance_()
  {
    ProFormaTokenizer::Token tok = current_();
    has_current_ = false;
    return tok;
  }

  bool ProFormaParser::check_(ProFormaTokenizer::TokenType type)
  {
    return current_().type == type;
  }

  bool ProFormaParser::match_(ProFormaTokenizer::TokenType type)
  {
    if (check_(type))
    {
      advance_();
      return true;
    }
    return false;
  }

  ProFormaTokenizer::Token ProFormaParser::expect_(ProFormaTokenizer::TokenType type, const char* expected_desc)
  {
    ProFormaTokenizer::Token tok = current_();
    if (tok.type != type)
    {
      std::string msg = std::string("Expected ") + expected_desc;
      errorAt_(ProFormaErrorCode::UNEXPECTED_CHARACTER, tok.position, msg.c_str());
    }
    return advance_();
  }

  bool ProFormaParser::isAtEnd_()
  {
    return current_().type == ProFormaTokenizer::TokenType::END;
  }

  void ProFormaParser::error_(ProFormaErrorCode code, const char* message)
  {
    errorAt_(code, current_().position, message);
  }

  void ProFormaParser::errorAt_(ProFormaErrorCode code, size_t pos, const char* message)
  {
    throw ProFormaParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                             code, pos, input_, message);
  }

  std::optional<CvDatabase> ProFormaParser::parseCvDatabasePrefix_(const std::string_view& id)
  {
    if (id == "UNIMOD") return CvDatabase::UNIMOD;
    if (id == "MOD") return CvDatabase::MOD;
    if (id == "RESID") return CvDatabase::RESID;
    if (id == "XLMOD") return CvDatabase::XLMOD;
    if (id == "GNO") return CvDatabase::GNO;
    return std::nullopt;
  }

  bool ProFormaParser::isAminoAcid_(char c)
  {
    // Standard amino acid one-letter codes (upper and lower case)
    // A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    // Plus: U (selenocysteine), O (pyrrolysine), X (unknown), B (Asx), Z (Glx), J (Xle)
    // ProForma allows lowercase amino acids
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
  }

  bool ProFormaParser::looksLikeModificationTagContent_()
  {
    ProFormaTokenizer::Token tok = current_();
    return tok.type == ProFormaTokenizer::TokenType::IDENTIFIER ||
           tok.type == ProFormaTokenizer::TokenType::NUMBER ||
           tok.type == ProFormaTokenizer::TokenType::PLUS ||
           tok.type == ProFormaTokenizer::TokenType::MINUS;
  }

  // ============================================================================
  // ToString methods (delegate to ProFormaWriter)
  // ============================================================================

  String ProFormaParser::toString(const Peptidoform& pf, ProFormaWriteMode mode)
  {
    return ProFormaWriter::toString(pf, mode);
  }

  String ProFormaParser::toString(const PeptidoformIon& pfi, ProFormaWriteMode mode)
  {
    return ProFormaWriter::toString(pfi, mode);
  }

  // ============================================================================
  // AASequence Conversion Methods
  // ============================================================================

  namespace
  {
    // Helper to resolve a single modification tag to a ResidueModification
    const ResidueModification* resolveModificationTag_(
      const ModificationTag& tag,
      char residue = '\0',
      ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY)
    {
      // Need non-const pointer because getBestModificationByDiffMonoMass is non-const
      ModificationsDB* mod_db = ModificationsDB::getInstance();

      return std::visit([&](auto&& arg) -> const ResidueModification* {
        using T = std::decay_t<decltype(arg)>;

        if constexpr (std::is_same_v<T, CvAccession>)
        {
          // Look up by CV accession
          String full_accession;
          switch (arg.database)
          {
            case CvDatabase::UNIMOD:
              full_accession = "UNIMOD:" + arg.accession;
              break;
            case CvDatabase::MOD:
              full_accession = "MOD:" + arg.accession;
              break;
            case CvDatabase::RESID:
              full_accession = "RESID:" + arg.accession;
              break;
            case CvDatabase::XLMOD:
              full_accession = "XLMOD:" + arg.accession;
              break;
            case CvDatabase::GNO:
              full_accession = "GNO:" + arg.accession;
              break;
          }

          try
          {
            // Pass residue to get correct modification variant (e.g., Oxidation on M vs D)
            String residue_str = (residue != '\0') ? String(1, residue) : "";
            return mod_db->getModification(full_accession, residue_str, term_spec);
          }
          catch (const Exception::ElementNotFound&)
          {
            return nullptr;
          }
        }
        else if constexpr (std::is_same_v<T, NamedMod>)
        {
          // Look up by name using searchModificationsFast
          bool multiple_matches = false;
          String residue_str = (residue != '\0') ? String(1, residue) : "";
          return mod_db->searchModificationsFast(arg.name, multiple_matches, residue_str, term_spec);
        }
        else if constexpr (std::is_same_v<T, MassDelta>)
        {
          // Look up by mass delta
          String residue_str = (residue != '\0') ? String(1, residue) : "";
          return mod_db->getBestModificationByDiffMonoMass(
            arg.mass, 0.01, residue_str, term_spec);
        }
        else if constexpr (std::is_same_v<T, FormulaTag>)
        {
          // Calculate mass from formula and look up
          // For now, we don't resolve formulas - they're used as-is
          return nullptr;
        }
        else if constexpr (std::is_same_v<T, GlycanComposition>)
        {
          // Glycans are complex and typically not in ModificationsDB
          return nullptr;
        }
        else if constexpr (std::is_same_v<T, InfoTag>)
        {
          // Info tags have no associated modification
          return nullptr;
        }
        else
        {
          return nullptr;
        }
      }, tag);
    }

    // Helper to resolve a modification (taking first alternative)
    void resolveModification_(Modification& mod, char residue, ResidueModification::TermSpecificity term_spec)
    {
      if (mod.alternatives.empty())
      {
        return;
      }

      // Try to resolve the first (primary) alternative
      const auto& [tag, label] = mod.alternatives[0];
      mod.resolved_mod = resolveModificationTag_(tag, residue, term_spec);
    }
  }

  void ProFormaParser::resolveModifications(Peptidoform& pf)
  {
    // Resolve sequence modifications
    for (auto& section : pf.sequence)
    {
      if (auto* elem = std::get_if<SequenceElement>(&section))
      {
        for (auto& mod : elem->modifications)
        {
          resolveModification_(mod, elem->amino_acid, ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
        }
      }
      else if (auto* region = std::get_if<AmbiguousRegion>(&section))
      {
        for (auto& elem : region->elements)
        {
          for (auto& mod : elem.modifications)
          {
            resolveModification_(mod, elem.amino_acid, ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
          }
        }
      }
      else if (auto* range = std::get_if<ModifiedRange>(&section))
      {
        for (auto& mod : range->modifications)
        {
          // Range modifications don't have a specific residue
          resolveModification_(mod, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
        }
      }
    }

    // Resolve N-terminal modifications
    for (auto& mod : pf.n_term_mods)
    {
      resolveModification_(mod, '\0', ResidueModification::N_TERM);
    }

    // Resolve C-terminal modifications
    for (auto& mod : pf.c_term_mods)
    {
      resolveModification_(mod, '\0', ResidueModification::C_TERM);
    }

    // Resolve unlocalised modifications
    for (auto& um : pf.unlocalised_mods)
    {
      for (auto& mod : um.modifications)
      {
        resolveModification_(mod, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
      }
    }

    // Resolve labile modifications
    for (auto& lm : pf.labile_mods)
    {
      resolveModification_(lm.modification, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    }

    // Resolve global modifications
    for (auto& entry : pf.global_mods)
    {
      if (auto* gm = std::get_if<GlobalModification>(&entry))
      {
        resolveModification_(gm->modification, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
      }
    }
  }

  std::vector<ConversionIssue> ProFormaParser::getAASequenceConversionIssues(const Peptidoform& pf)
  {
    std::vector<ConversionIssue> issues;

    // Make a copy to resolve modifications (don't modify original)
    Peptidoform pf_copy = pf;
    resolveModifications(pf_copy);

    // Check for unlocalised modifications
    if (!pf_copy.unlocalised_mods.empty())
    {
      issues.push_back({
        ConversionIssueType::UNLOCALISED_MOD,
        "Peptidoform contains unlocalised modifications that cannot be represented in AASequence",
        SIZE_MAX
      });
    }

    // Check for labile modifications
    if (!pf_copy.labile_mods.empty())
    {
      issues.push_back({
        ConversionIssueType::LABILE_MOD,
        "Peptidoform contains labile modifications that cannot be represented in AASequence",
        SIZE_MAX
      });
    }

    // Check for global modifications
    for (const auto& entry : pf_copy.global_mods)
    {
      if (std::holds_alternative<GlobalModification>(entry))
      {
        issues.push_back({
          ConversionIssueType::GLOBAL_MOD,
          "Peptidoform contains global modifications that cannot be represented in AASequence",
          SIZE_MAX
        });
        break;
      }
    }

    // Check sequence sections for issues
    size_t position = 0;
    for (const auto& section : pf_copy.sequence)
    {
      if (std::holds_alternative<AmbiguousRegion>(section))
      {
        issues.push_back({
          ConversionIssueType::AMBIGUOUS_REGION,
          "Peptidoform contains ambiguous amino acid region at position " + std::to_string(position),
          position
        });
      }
      else if (std::holds_alternative<ModifiedRange>(section))
      {
        issues.push_back({
          ConversionIssueType::MODIFIED_RANGE,
          "Peptidoform contains modified range at position " + std::to_string(position),
          position
        });
      }
      else if (const auto* elem = std::get_if<SequenceElement>(&section))
      {
        for (const auto& mod : elem->modifications)
        {
          // Check for unresolved modifications
          if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
          {
            // Check if the first alternative is a meaningful tag (not empty InfoTag)
            const auto& [tag, label] = mod.alternatives[0];
            bool is_empty_info = std::holds_alternative<InfoTag>(tag) &&
                                 std::get<InfoTag>(tag).text.empty();
            if (!is_empty_info)
            {
              issues.push_back({
                ConversionIssueType::UNRESOLVED_MOD,
                "Modification at position " + std::to_string(position) + " could not be resolved",
                position
              });
            }
          }

          // Check for alternative modifications
          if (mod.alternatives.size() > 1)
          {
            issues.push_back({
              ConversionIssueType::ALTERNATIVE_MODS,
              "Modification at position " + std::to_string(position) + " has multiple alternatives",
              position
            });
          }

          // Check for cross-links
          for (const auto& [tag, label] : mod.alternatives)
          {
            if (label.has_value() && label->type == Label::Type::CROSSLINK)
            {
              issues.push_back({
                ConversionIssueType::CROSS_LINK,
                "Modification at position " + std::to_string(position) + " is part of a cross-link",
                position
              });
              break;
            }
          }
        }
      }

      position++;
    }

    // Check N-terminal modifications
    for (const auto& mod : pf_copy.n_term_mods)
    {
      if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
      {
        const auto& [tag, label] = mod.alternatives[0];
        bool is_empty_info = std::holds_alternative<InfoTag>(tag) &&
                             std::get<InfoTag>(tag).text.empty();
        if (!is_empty_info)
        {
          issues.push_back({
            ConversionIssueType::UNRESOLVED_MOD,
            "N-terminal modification could not be resolved",
            SIZE_MAX
          });
        }
      }
      if (mod.alternatives.size() > 1)
      {
        issues.push_back({
          ConversionIssueType::ALTERNATIVE_MODS,
          "N-terminal modification has multiple alternatives",
          SIZE_MAX
        });
      }
    }

    // Check C-terminal modifications
    for (const auto& mod : pf_copy.c_term_mods)
    {
      if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
      {
        const auto& [tag, label] = mod.alternatives[0];
        bool is_empty_info = std::holds_alternative<InfoTag>(tag) &&
                             std::get<InfoTag>(tag).text.empty();
        if (!is_empty_info)
        {
          issues.push_back({
            ConversionIssueType::UNRESOLVED_MOD,
            "C-terminal modification could not be resolved",
            SIZE_MAX
          });
        }
      }
      if (mod.alternatives.size() > 1)
      {
        issues.push_back({
          ConversionIssueType::ALTERNATIVE_MODS,
          "C-terminal modification has multiple alternatives",
          SIZE_MAX
        });
      }
    }

    return issues;
  }

  bool ProFormaParser::isRepresentableAsAASequence(const Peptidoform& pf)
  {
    return getAASequenceConversionIssues(pf).empty();
  }

  AASequence ProFormaParser::toAASequence(
    const Peptidoform& pf,
    AASequenceConversionPolicy policy)
  {
    // Get conversion issues
    std::vector<ConversionIssue> issues = getAASequenceConversionIssues(pf);

    // Check policy
    if (policy == AASequenceConversionPolicy::STRICT_MODE && !issues.empty())
    {
      std::string error_msg = "Cannot convert Peptidoform to AASequence: ";
      for (const auto& issue : issues)
      {
        error_msg += issue.description + "; ";
      }
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
    }

    // Make a copy and resolve modifications
    Peptidoform pf_copy = pf;
    resolveModifications(pf_copy);

    // Build unmodified sequence string
    std::string unmod_seq;
    for (const auto& section : pf_copy.sequence)
    {
      if (const auto* elem = std::get_if<SequenceElement>(&section))
      {
        unmod_seq += elem->amino_acid;
      }
      else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
      {
        // For ambiguous regions, use the first amino acid
        if (!region->elements.empty())
        {
          unmod_seq += region->elements[0].amino_acid;
        }
      }
      else if (const auto* range = std::get_if<ModifiedRange>(&section))
      {
        // Add all amino acids from the range
        for (const auto& elem : range->elements)
        {
          unmod_seq += elem.amino_acid;
        }
      }
    }

    // Create AASequence from unmodified string
    AASequence seq = AASequence::fromString(unmod_seq);

    // Apply modifications
    size_t seq_pos = 0;
    for (const auto& section : pf_copy.sequence)
    {
      if (const auto* elem = std::get_if<SequenceElement>(&section))
      {
        for (const auto& mod : elem->modifications)
        {
          if (mod.resolved_mod != nullptr)
          {
            seq.setModification(seq_pos, mod.resolved_mod);
          }
          else if (policy == AASequenceConversionPolicy::STRICT_MODE)
          {
            // Already checked above, but double-check
            throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Unresolved modification at position " + std::to_string(seq_pos));
          }
          // For BEST_EFFORT and DROP_UNLOCALISED, skip unresolved
        }
        seq_pos++;
      }
      else if (std::holds_alternative<AmbiguousRegion>(section))
      {
        seq_pos++; // Counted as one position above
      }
      else if (const auto* range = std::get_if<ModifiedRange>(&section))
      {
        // Skip modifications for ranges in non-STRICT mode
        seq_pos += range->elements.size();
      }
    }

    // Apply N-terminal modifications (use the first one if multiple)
    if (!pf_copy.n_term_mods.empty() && pf_copy.n_term_mods[0].resolved_mod != nullptr)
    {
      seq.setNTerminalModification(pf_copy.n_term_mods[0].resolved_mod);
    }

    // Apply C-terminal modifications (use the first one if multiple)
    if (!pf_copy.c_term_mods.empty() && pf_copy.c_term_mods[0].resolved_mod != nullptr)
    {
      seq.setCTerminalModification(pf_copy.c_term_mods[0].resolved_mod);
    }

    return seq;
  }

  Peptidoform ProFormaParser::fromAASequence(const AASequence& seq)
  {
    Peptidoform pf;

    // Build sequence
    for (Size i = 0; i < seq.size(); ++i)
    {
      SequenceElement elem;
      elem.amino_acid = seq[i].getOneLetterCode()[0];

      // Check for modifications on this residue
      if (seq[i].isModified())
      {
        const ResidueModification* mod = seq[i].getModification();
        if (mod != nullptr)
        {
          Modification proforma_mod;

          // Prefer UNIMOD accession if available
          String unimod_acc = mod->getUniModAccession();
          if (!unimod_acc.empty() && unimod_acc.hasPrefix("UniMod:"))
          {
            CvAccession cv;
            cv.database = CvDatabase::UNIMOD;
            cv.accession = unimod_acc.substr(7); // Remove "UniMod:" prefix
            proforma_mod.alternatives.emplace_back(std::move(cv), std::nullopt);
          }
          else
          {
            // Use modification name
            NamedMod nm;
            nm.name = mod->getId();
            nm.cv_hint = std::nullopt;
            proforma_mod.alternatives.emplace_back(std::move(nm), std::nullopt);
          }

          proforma_mod.resolved_mod = mod;
          elem.modifications.push_back(std::move(proforma_mod));
        }
      }

      pf.sequence.push_back(std::move(elem));
    }

    // Check for N-terminal modification
    if (seq.hasNTerminalModification())
    {
      const ResidueModification* mod = seq.getNTerminalModification();
      if (mod != nullptr)
      {
        Modification proforma_mod;

        String unimod_acc = mod->getUniModAccession();
        if (!unimod_acc.empty() && unimod_acc.hasPrefix("UniMod:"))
        {
          CvAccession cv;
          cv.database = CvDatabase::UNIMOD;
          cv.accession = unimod_acc.substr(7);
          proforma_mod.alternatives.emplace_back(std::move(cv), std::nullopt);
        }
        else
        {
          NamedMod nm;
          nm.name = mod->getId();
          nm.cv_hint = std::nullopt;
          proforma_mod.alternatives.emplace_back(std::move(nm), std::nullopt);
        }

        proforma_mod.resolved_mod = mod;
        pf.n_term_mods.push_back(std::move(proforma_mod));
      }
    }

    // Check for C-terminal modification
    if (seq.hasCTerminalModification())
    {
      const ResidueModification* mod = seq.getCTerminalModification();
      if (mod != nullptr)
      {
        Modification proforma_mod;

        String unimod_acc = mod->getUniModAccession();
        if (!unimod_acc.empty() && unimod_acc.hasPrefix("UniMod:"))
        {
          CvAccession cv;
          cv.database = CvDatabase::UNIMOD;
          cv.accession = unimod_acc.substr(7);
          proforma_mod.alternatives.emplace_back(std::move(cv), std::nullopt);
        }
        else
        {
          NamedMod nm;
          nm.name = mod->getId();
          nm.cv_hint = std::nullopt;
          proforma_mod.alternatives.emplace_back(std::move(nm), std::nullopt);
        }

        proforma_mod.resolved_mod = mod;
        pf.c_term_mods.push_back(std::move(proforma_mod));
      }
    }

    return pf;
  }

  // ============================================================================
  // Mass Calculation Methods
  // ============================================================================

  namespace
  {
    /// Helper to get modification mass from a Modification struct
    /// Returns pair<bool, double> where first indicates if mass is available
    std::pair<bool, double> getModificationMass_(const Modification& mod)
    {
      // If resolved, use the resolved modification's mass
      if (mod.resolved_mod != nullptr)
      {
        return {true, mod.resolved_mod->getDiffMonoMass()};
      }

      // Check first alternative for mass delta or formula
      if (mod.alternatives.empty())
      {
        return {false, 0.0};
      }

      const auto& tag = mod.alternatives[0].first;

      // MassDelta has explicit mass
      if (const auto* md = std::get_if<MassDelta>(&tag))
      {
        return {true, md->mass};
      }

      // FormulaTag can be converted to mass
      if (const auto* ft = std::get_if<FormulaTag>(&tag))
      {
        try
        {
          EmpiricalFormula formula(ft->formula_string);
          return {true, formula.getMonoWeight()};
        }
        catch (...)
        {
          return {false, 0.0};
        }
      }

      // InfoTag and PositionConstraint don't contribute mass
      if (std::holds_alternative<InfoTag>(tag) ||
          std::holds_alternative<PositionConstraint>(tag))
      {
        return {true, 0.0};
      }

      // Other types (CvAccession, NamedMod, GlycanComposition) need resolution
      return {false, 0.0};
    }

    /// Helper to check a single modification and collect issues
    void checkModificationForMass_(const Modification& mod, size_t position,
                                   std::vector<ConversionIssue>& issues)
    {
      auto [has_mass, mass] = getModificationMass_(mod);
      (void)mass; // unused in check

      if (!has_mass)
      {
        issues.push_back({
          ConversionIssueType::UNRESOLVED_MOD,
          "Modification at position " + std::to_string(position) + " has no resolvable mass",
          position
        });
      }
    }

    /// Helper to calculate mass for a single chain, sharing cross-link tracking across chains
    /// @param pf_resolved A Peptidoform with modifications already resolved
    /// @param counted_crosslinks Shared set to track which cross-link labels have been counted
    double calculateChainMass_(const Peptidoform& pf_resolved, std::set<String>& counted_crosslinks)
    {
      double mass = 0.0;

      // Helper lambda to add modification mass, handling cross-links
      auto addModMass = [&](const Modification& mod) {
        // Check for cross-link label - only count mass once per label
        if (!mod.alternatives.empty() && mod.alternatives[0].second.has_value())
        {
          const auto& label = mod.alternatives[0].second.value();
          if (label.type == Label::Type::CROSSLINK)
          {
            if (counted_crosslinks.count(label.identifier) > 0)
            {
              return; // Already counted
            }
            counted_crosslinks.insert(label.identifier);
          }
        }

        auto [has_mass, mod_mass] = getModificationMass_(mod);
        if (has_mass)
        {
          mass += mod_mass;
        }
      };

      // Calculate sequence mass
      for (const auto& section : pf_resolved.sequence)
      {
        if (const auto* elem = std::get_if<SequenceElement>(&section))
        {
          const Residue* res = ResidueDB::getInstance()->getResidue(elem->amino_acid);
          mass += res->getMonoWeight(Residue::Internal);

          for (const auto& mod : elem->modifications)
          {
            addModMass(mod);
          }
        }
        else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
        {
          if (!region->elements.empty())
          {
            const Residue* res = ResidueDB::getInstance()->getResidue(region->elements[0].amino_acid);
            mass += res->getMonoWeight(Residue::Internal);

            for (const auto& mod : region->elements[0].modifications)
            {
              addModMass(mod);
            }
          }
        }
        else if (const auto* range = std::get_if<ModifiedRange>(&section))
        {
          for (const auto& elem : range->elements)
          {
            const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
            mass += res->getMonoWeight(Residue::Internal);
          }
          for (const auto& mod : range->modifications)
          {
            addModMass(mod);
          }
        }
      }

      // Add terminal masses (H at N-term, OH at C-term)
      mass += EmpiricalFormula("H2O").getMonoWeight();

      // Add terminal modifications
      for (const auto& mod : pf_resolved.n_term_mods)
      {
        addModMass(mod);
      }
      for (const auto& mod : pf_resolved.c_term_mods)
      {
        addModMass(mod);
      }

      // Add unlocalised modifications (with optional occurrence count)
      for (const auto& um : pf_resolved.unlocalised_mods)
      {
        int occurrence = um.occurrence.value_or(1);
        for (const auto& mod : um.modifications)
        {
          auto [has_mass, mod_mass] = getModificationMass_(mod);
          if (has_mass)
          {
            mass += mod_mass * occurrence;
          }
        }
      }

      // Add labile modifications
      for (const auto& lm : pf_resolved.labile_mods)
      {
        auto [has_mass, mod_mass] = getModificationMass_(lm.modification);
        if (has_mass)
        {
          mass += mod_mass;
        }
      }

      // Handle global modifications
      for (const auto& entry : pf_resolved.global_mods)
      {
        if (const auto* gm = std::get_if<GlobalModification>(&entry))
        {
          auto [has_mass, mod_mass] = getModificationMass_(gm->modification);
          if (has_mass)
          {
            // Count matching residues
            int count = 0;
            for (const auto& section : pf_resolved.sequence)
            {
              if (const auto* elem = std::get_if<SequenceElement>(&section))
              {
                for (const String& loc : gm->locations)
                {
                  if (loc.size() == 1 && elem->amino_acid == loc[0])
                  {
                    ++count;
                    break;
                  }
                }
              }
            }
            mass += mod_mass * count;
          }
        }
      }

      return mass;
    }
  } // anonymous namespace

  std::vector<ConversionIssue> ProFormaParser::getMassCalculationIssues(const Peptidoform& pf)
  {
    std::vector<ConversionIssue> issues;

    // Make a copy and resolve modifications
    Peptidoform pf_copy = pf;
    resolveModifications(pf_copy);

    // Check sequence elements
    size_t position = 0;
    for (const auto& section : pf_copy.sequence)
    {
      if (const auto* elem = std::get_if<SequenceElement>(&section))
      {
        // Check if amino acid is valid
        const Residue* res = ResidueDB::getInstance()->getResidue(elem->amino_acid);
        if (res == nullptr)
        {
          issues.push_back({
            ConversionIssueType::UNSUPPORTED_FEATURE,
            String("Unknown amino acid '") + elem->amino_acid + "' at position " + String(position),
            position
          });
        }

        // Check modifications
        for (const auto& mod : elem->modifications)
        {
          checkModificationForMass_(mod, position, issues);
        }
        ++position;
      }
      else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
      {
        // Check if all amino acids in ambiguous region have the same mass
        std::set<double> masses;
        for (const auto& elem : region->elements)
        {
          const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
          if (res != nullptr)
          {
            masses.insert(res->getMonoWeight(Residue::Internal));
          }
          else
          {
            issues.push_back({
              ConversionIssueType::UNSUPPORTED_FEATURE,
              String("Unknown amino acid '") + elem.amino_acid + "' in ambiguous region at position " + String(position),
              position
            });
          }
        }
        if (masses.size() > 1)
        {
          issues.push_back({
            ConversionIssueType::AMBIGUOUS_REGION,
            "Ambiguous region at position " + std::to_string(position) +
            " contains amino acids with different masses",
            position
          });
        }
        ++position;
      }
      else if (const auto* range = std::get_if<ModifiedRange>(&section))
      {
        for (const auto& elem : range->elements)
        {
          const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
          if (res == nullptr)
          {
            issues.push_back({
              ConversionIssueType::UNSUPPORTED_FEATURE,
              String("Unknown amino acid '") + elem.amino_acid + "' in range at position " + String(position),
              position
            });
          }
          ++position;
        }
        for (const auto& mod : range->modifications)
        {
          checkModificationForMass_(mod, position - range->elements.size(), issues);
        }
      }
    }

    // Check terminal modifications
    for (const auto& mod : pf_copy.n_term_mods)
    {
      checkModificationForMass_(mod, 0, issues);
    }
    for (const auto& mod : pf_copy.c_term_mods)
    {
      checkModificationForMass_(mod, position > 0 ? position - 1 : 0, issues);
    }

    // Check unlocalised modifications
    for (const auto& um : pf_copy.unlocalised_mods)
    {
      for (const auto& mod : um.modifications)
      {
        checkModificationForMass_(mod, SIZE_MAX, issues);
      }
    }

    // Check labile modifications
    for (const auto& lm : pf_copy.labile_mods)
    {
      checkModificationForMass_(lm.modification, SIZE_MAX, issues);
    }

    // Check global modifications
    for (const auto& entry : pf_copy.global_mods)
    {
      if (const auto* gm = std::get_if<GlobalModification>(&entry))
      {
        checkModificationForMass_(gm->modification, SIZE_MAX, issues);
      }
      // IsotopeReplacement doesn't prevent mass calculation but changes mass
      // We can handle isotope replacement by adjusting elemental masses
    }

    return issues;
  }

  std::vector<ConversionIssue> ProFormaParser::getMassCalculationIssues(const PeptidoformIon& pfi)
  {
    std::vector<ConversionIssue> issues;

    for (size_t i = 0; i < pfi.chains.size(); ++i)
    {
      auto chain_issues = getMassCalculationIssues(pfi.chains[i]);
      for (auto& issue : chain_issues)
      {
        issue.description = "Chain " + std::to_string(i) + ": " + issue.description;
        issues.push_back(std::move(issue));
      }
    }

    return issues;
  }

  bool ProFormaParser::canCalculateMass(const Peptidoform& pf)
  {
    return getMassCalculationIssues(pf).empty();
  }

  bool ProFormaParser::canCalculateMass(const PeptidoformIon& pfi)
  {
    return getMassCalculationIssues(pfi).empty();
  }

  double ProFormaParser::getMonoWeight(const Peptidoform& pf)
  {
    auto issues = getMassCalculationIssues(pf);
    if (!issues.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate mass: " + issues[0].description, "");
    }

    Peptidoform pf_copy = pf;
    resolveModifications(pf_copy);

    std::set<String> counted_crosslinks;
    return calculateChainMass_(pf_copy, counted_crosslinks);
  }

  double ProFormaParser::getMonoWeight(const PeptidoformIon& pfi)
  {
    auto issues = getMassCalculationIssues(pfi);
    if (!issues.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate mass: " + issues[0].description, "");
    }

    if (pfi.chains.empty())
    {
      return 0.0;
    }

    // For chimeric spectra, mass is ambiguous - each chain is an independent peptide
    if (pfi.is_chimeric)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate single mass for chimeric spectra. "
        "Use getMonoWeight() on individual chains.", "");
    }

    // For cross-linked peptides, share cross-link tracking so mass is counted once
    std::set<String> counted_crosslinks;
    double total = 0.0;
    for (const auto& chain : pfi.chains)
    {
      Peptidoform chain_copy = chain;
      resolveModifications(chain_copy);
      total += calculateChainMass_(chain_copy, counted_crosslinks);
    }
    return total;
  }

  double ProFormaParser::getMZ(const PeptidoformIon& pfi)
  {
    if (!pfi.charge.has_value())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate m/z: no charge state specified", "");
    }

    // ChargeState is a variant: either int or vector<AdductIon>
    int charge = 0;
    if (const int* simple_charge = std::get_if<int>(&pfi.charge.value()))
    {
      charge = *simple_charge;
    }
    else if (const auto* adducts = std::get_if<std::vector<AdductIon>>(&pfi.charge.value()))
    {
      // Sum the charges from all adducts
      for (const auto& adduct : *adducts)
      {
        int occurrence = adduct.occurrence.value_or(1);
        charge += adduct.charge * occurrence;
      }
    }

    if (charge == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate m/z: charge state is zero", "");
    }

    double mass = getMonoWeight(pfi);
    return (mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
  }

  double ProFormaParser::getMZ(const Peptidoform& pf, int charge)
  {
    if (charge == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot calculate m/z: charge state is zero", "");
    }

    double mass = getMonoWeight(pf);
    return (mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
  }

  // ============================================================================
  // Non-throwing variants (single-pass, efficient)
  // ============================================================================

  std::optional<double> ProFormaParser::tryGetMonoWeight(const Peptidoform& pf)
  {
    std::vector<ConversionIssue> issues;
    return tryGetMonoWeight(pf, issues);
  }

  std::optional<double> ProFormaParser::tryGetMonoWeight(const Peptidoform& pf,
                                                          std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    // Single pass: resolve and collect issues
    Peptidoform pf_copy = pf;
    resolveModifications(pf_copy);

    // Collect issues (reuse existing logic but on resolved copy)
    issues_out = getMassCalculationIssues(pf_copy);

    if (!issues_out.empty())
    {
      return std::nullopt;
    }

    // Calculate mass using resolved copy
    std::set<String> counted_crosslinks;
    return calculateChainMass_(pf_copy, counted_crosslinks);
  }

  std::optional<double> ProFormaParser::tryGetMonoWeight(const PeptidoformIon& pfi)
  {
    std::vector<ConversionIssue> issues;
    return tryGetMonoWeight(pfi, issues);
  }

  std::optional<double> ProFormaParser::tryGetMonoWeight(const PeptidoformIon& pfi,
                                                          std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    if (pfi.chains.empty())
    {
      return 0.0;
    }

    // Chimeric spectra have ambiguous mass
    if (pfi.is_chimeric)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Cannot calculate single mass for chimeric spectra. Use tryGetMonoWeight() on individual chains.",
        0
      });
      return std::nullopt;
    }

    // Check all chains for issues first
    for (size_t i = 0; i < pfi.chains.size(); ++i)
    {
      Peptidoform chain_copy = pfi.chains[i];
      resolveModifications(chain_copy);

      auto chain_issues = getMassCalculationIssues(chain_copy);
      for (auto& issue : chain_issues)
      {
        issue.description = "Chain " + std::to_string(i) + ": " + issue.description;
        issues_out.push_back(std::move(issue));
      }
    }

    if (!issues_out.empty())
    {
      return std::nullopt;
    }

    // Calculate mass for cross-linked peptides with shared tracking
    std::set<String> counted_crosslinks;
    double total = 0.0;
    for (const auto& chain : pfi.chains)
    {
      Peptidoform chain_copy = chain;
      resolveModifications(chain_copy);
      total += calculateChainMass_(chain_copy, counted_crosslinks);
    }
    return total;
  }

  std::optional<double> ProFormaParser::tryGetMZ(const Peptidoform& pf, int charge)
  {
    std::vector<ConversionIssue> issues;
    return tryGetMZ(pf, charge, issues);
  }

  std::optional<double> ProFormaParser::tryGetMZ(const Peptidoform& pf, int charge,
                                                  std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    if (charge == 0)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Charge state is zero",
        0
      });
      return std::nullopt;
    }

    auto mass = tryGetMonoWeight(pf, issues_out);
    if (!mass.has_value())
    {
      return std::nullopt;
    }

    return (*mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
  }

  std::optional<double> ProFormaParser::tryGetMZ(const PeptidoformIon& pfi)
  {
    std::vector<ConversionIssue> issues;
    return tryGetMZ(pfi, issues);
  }

  std::optional<double> ProFormaParser::tryGetMZ(const PeptidoformIon& pfi,
                                                  std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    if (!pfi.charge.has_value())
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "No charge state specified",
        0
      });
      return std::nullopt;
    }

    // Extract charge from variant
    int charge = 0;
    if (const int* simple_charge = std::get_if<int>(&pfi.charge.value()))
    {
      charge = *simple_charge;
    }
    else if (const auto* adducts = std::get_if<std::vector<AdductIon>>(&pfi.charge.value()))
    {
      for (const auto& adduct : *adducts)
      {
        int occurrence = adduct.occurrence.value_or(1);
        charge += adduct.charge * occurrence;
      }
    }

    if (charge == 0)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Charge state is zero",
        0
      });
      return std::nullopt;
    }

    auto mass = tryGetMonoWeight(pfi, issues_out);
    if (!mass.has_value())
    {
      return std::nullopt;
    }

    return (*mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
  }

  // ============================================================================
  // Theoretical Spectrum Generation
  // ============================================================================

  std::optional<MSSpectrum> ProFormaParser::tryGenerateSpectrum(
    const Peptidoform& pf,
    int min_charge,
    int max_charge,
    bool add_losses,
    bool add_metainfo)
  {
    std::vector<ConversionIssue> issues;
    return tryGenerateSpectrum(pf, min_charge, max_charge, add_losses, add_metainfo, issues);
  }

  std::optional<MSSpectrum> ProFormaParser::tryGenerateSpectrum(
    const Peptidoform& pf,
    int min_charge,
    int max_charge,
    bool add_losses,
    bool add_metainfo,
    std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    // Check if convertible to AASequence
    if (!isRepresentableAsAASequence(pf))
    {
      auto conv_issues = getAASequenceConversionIssues(pf);
      for (auto& issue : conv_issues)
      {
        issues_out.push_back(std::move(issue));
      }
      return std::nullopt;
    }

    // Convert to AASequence
    AASequence seq;
    try
    {
      seq = toAASequence(pf, AASequenceConversionPolicy::STRICT_MODE);
    }
    catch (const Exception::ConversionError& e)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        String("AASequence conversion failed: ") + e.what(),
        0
      });
      return std::nullopt;
    }

    // Configure the generator
    TheoreticalSpectrumGenerator generator;
    Param param = generator.getParameters();
    param.setValue("add_b_ions", "true");
    param.setValue("add_y_ions", "true");
    param.setValue("add_losses", add_losses ? "true" : "false");
    param.setValue("add_metainfo", add_metainfo ? "true" : "false");
    generator.setParameters(param);

    // Generate spectrum
    MSSpectrum spectrum;
    try
    {
      generator.getSpectrum(spectrum, seq, min_charge, max_charge);
    }
    catch (const Exception::BaseException& e)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        String("Spectrum generation failed: ") + e.what(),
        0
      });
      return std::nullopt;
    }

    return spectrum;
  }

  std::optional<MSSpectrum> ProFormaParser::tryGenerateSpectrum(
    const PeptidoformIon& pfi,
    int min_charge,
    int max_charge,
    bool add_losses,
    bool add_metainfo)
  {
    std::vector<ConversionIssue> issues;
    return tryGenerateSpectrum(pfi, min_charge, max_charge, add_losses, add_metainfo, issues);
  }

  std::optional<MSSpectrum> ProFormaParser::tryGenerateSpectrum(
    const PeptidoformIon& pfi,
    int min_charge,
    int max_charge,
    bool add_losses,
    bool add_metainfo,
    std::vector<ConversionIssue>& issues_out)
  {
    issues_out.clear();

    if (pfi.chains.empty())
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "No peptide chains to fragment",
        0
      });
      return std::nullopt;
    }

    // Chimeric spectra not supported
    if (pfi.is_chimeric)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Theoretical spectrum generation not supported for chimeric spectra. Generate spectra for individual chains.",
        0
      });
      return std::nullopt;
    }

    // Single chain - delegate to Peptidoform version
    if (pfi.chains.size() == 1)
    {
      return tryGenerateSpectrum(pfi.chains[0], min_charge, max_charge, add_losses, add_metainfo, issues_out);
    }

    // Multiple chains = cross-linked peptides
    // Need to extract cross-link information
    if (pfi.chains.size() != 2)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Only two-chain cross-links are currently supported for spectrum generation",
        0
      });
      return std::nullopt;
    }

    // Find cross-link positions and mass in both chains
    const Peptidoform& alpha_pf = pfi.chains[0];
    const Peptidoform& beta_pf = pfi.chains[1];

    // Helper to find cross-link position and mass in a chain
    auto findCrossLink = [](const Peptidoform& chain) -> std::tuple<bool, size_t, double, String> {
      size_t position = 0;
      for (const auto& section : chain.sequence)
      {
        if (const auto* elem = std::get_if<SequenceElement>(&section))
        {
          for (const auto& mod : elem->modifications)
          {
            if (!mod.alternatives.empty() && mod.alternatives[0].second.has_value())
            {
              const auto& label = mod.alternatives[0].second.value();
              if (label.type == Label::Type::CROSSLINK)
              {
                // Get the mass
                double mass = 0.0;
                if (const auto* md = std::get_if<MassDelta>(&mod.alternatives[0].first))
                {
                  mass = md->mass;
                }
                else if (mod.resolved_mod != nullptr)
                {
                  mass = mod.resolved_mod->getDiffMonoMass();
                }
                return {true, position, mass, label.identifier};
              }
            }
          }
          ++position;
        }
      }
      return {false, 0, 0.0, ""};
    };

    auto [alpha_found, alpha_pos, alpha_mass, alpha_label] = findCrossLink(alpha_pf);
    auto [beta_found, beta_pos, beta_mass, beta_label] = findCrossLink(beta_pf);

    if (!alpha_found || !beta_found)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Cross-link label not found in both chains",
        0
      });
      return std::nullopt;
    }

    if (alpha_label != beta_label)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        "Cross-link labels don't match between chains",
        0
      });
      return std::nullopt;
    }

    // Convert both chains to AASequence using BEST_EFFORT to skip cross-link mods
    // (the cross-linker mass is handled separately via OPXLDataStructs)
    AASequence alpha_seq, beta_seq;
    try
    {
      alpha_seq = toAASequence(alpha_pf, AASequenceConversionPolicy::BEST_EFFORT);
      beta_seq = toAASequence(beta_pf, AASequenceConversionPolicy::BEST_EFFORT);
    }
    catch (const Exception::ConversionError& e)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        String("AASequence conversion failed: ") + e.what(),
        0
      });
      return std::nullopt;
    }

    // Determine cross-linker mass (use the one with mass, typically first chain)
    double linker_mass = (alpha_mass > 0.001) ? alpha_mass : beta_mass;

    // Build cross-link structure
    OPXLDataStructs::ProteinProteinCrossLink crosslink;
    crosslink.alpha = &alpha_seq;
    crosslink.beta = &beta_seq;
    crosslink.cross_link_position = {static_cast<SignedSize>(alpha_pos), static_cast<SignedSize>(beta_pos)};
    crosslink.cross_linker_mass = linker_mass;
    crosslink.cross_linker_name = alpha_label;

    // Configure the XL generator
    TheoreticalSpectrumGeneratorXLMS generator;
    Param param = generator.getParameters();
    param.setValue("add_b_ions", "true");
    param.setValue("add_y_ions", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_losses", add_losses ? "true" : "false");
    param.setValue("add_metainfo", add_metainfo ? "true" : "false");
    param.setValue("add_precursor_peaks", "true");
    generator.setParameters(param);

    // Generate spectrum
    MSSpectrum spectrum;
    try
    {
      // Generate fragments from alpha chain
      generator.getXLinkIonSpectrum(spectrum, crosslink, true, min_charge, max_charge);
      // Generate fragments from beta chain
      generator.getXLinkIonSpectrum(spectrum, crosslink, false, min_charge, max_charge);
      // Sort by m/z
      spectrum.sortByPosition();
    }
    catch (const Exception::BaseException& e)
    {
      issues_out.push_back({
        ConversionIssueType::UNSUPPORTED_FEATURE,
        String("XL spectrum generation failed: ") + e.what(),
        0
      });
      return std::nullopt;
    }

    return spectrum;
  }

} // namespace OpenMS
