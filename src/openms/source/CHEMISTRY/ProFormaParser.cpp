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
#include <OpenMS/CONCEPT/Exception.h>

#include <cstdlib>
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

    // Parse first chain
    ion.chains.push_back(parsePeptidoform_());

    // Parse additional chains separated by //
    while (match_(ProFormaTokenizer::TokenType::SLASH))
    {
      if (!match_(ProFormaTokenizer::TokenType::SLASH))
      {
        // Single slash - this is charge state
        ion.charge = parseChargeState_();
        break;
      }
      // Double slash - parse another chain
      ion.chains.push_back(parsePeptidoform_());
    }

    // Verify we consumed all input
    if (!isAtEnd_())
    {
      error_(ProFormaErrorCode::UNEXPECTED_CHARACTER,
             "Unexpected characters after peptidoform ion");
    }

    return ion;
  }

  Peptidoform ProFormaParser::parsePeptidoform_()
  {
    Peptidoform pf;

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
    }

    return pf;
  }

  bool ProFormaParser::hasNTerminalModPattern_()
  {
    // Create a lookahead tokenizer at the current position
    ProFormaTokenizer lookahead(input_);

    // Skip to current position
    size_t target = has_current_ ? current_token_.position : tokenizer_.position();
    while (lookahead.position() < target && lookahead.hasMore())
    {
      lookahead.next();
    }

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
      ProFormaTokenizer lookahead = tokenizer_;
      lookahead.next(); // consume comma
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

    // Unlocalised mods look like: [mod]? or [mod]^2?
    // They appear before the sequence and end with ?

    while (check_(ProFormaTokenizer::TokenType::LBRACKET))
    {
      // Look ahead to see if this ends with ?
      ProFormaTokenizer lookahead(input_);
      size_t target = has_current_ ? current_token_.position : tokenizer_.position();
      while (lookahead.position() < target && lookahead.hasMore())
      {
        lookahead.next();
      }

      // Scan past the bracket content
      int depth = 0;
      bool found_question = false;
      ProFormaTokenizer::Token tok = lookahead.next(); // [
      depth = 1;

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

      // Check for ^N and/or ?
      tok = lookahead.peek();
      std::optional<int> occurrence;

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

      if (tok.type == ProFormaTokenizer::TokenType::QUESTION)
      {
        found_question = true;
      }

      if (!found_question)
      {
        // Not an unlocalised mod
        break;
      }

      // Now actually parse it
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

      expect_(ProFormaTokenizer::TokenType::QUESTION, "'?' for unlocalised modification");

      mods.push_back(std::move(um));
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
        ProFormaTokenizer lookahead(input_);
        size_t target = tok.position;
        while (lookahead.position() < target && lookahead.hasMore())
        {
          lookahead.next();
        }
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
    while (!check_(ProFormaTokenizer::TokenType::RPAREN) && !isAtEnd_())
    {
      ProFormaTokenizer::Token tok = current_();

      if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();

        for (char c : text)
        {
          if (!isAminoAcid_(c))
          {
            errorAt_(ProFormaErrorCode::INVALID_AMINO_ACID, tok.position,
                     "Invalid amino acid in range");
          }

          SequenceElement elem;
          elem.amino_acid = c;
          range.elements.push_back(elem);
        }
      }
      else
      {
        break;
      }
    }

    expect_(ProFormaTokenizer::TokenType::RPAREN, "')'");

    // Parse modifications for the range
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
    else if (id == "UNIMOD" || id == "MOD" || id == "RESID" || id == "XLMOD" || id == "GNO")
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
    // Tokenizer splits at letter/digit boundaries, so we concatenate tokens
    std::string name;

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
        // But we need to be careful - hyphens at end might be structural
        // Peek ahead to see if next token is part of the name
        ProFormaTokenizer lookahead = tokenizer_;
        ProFormaTokenizer::Token next = lookahead.peek();
        if (next.type == ProFormaTokenizer::TokenType::IDENTIFIER ||
            next.type == ProFormaTokenizer::TokenType::NUMBER)
        {
          name += "-";
          advance_();
        }
        else
        {
          break;
        }
      }
      else if (tok.text == " ")
      {
        // Space can be part of names like "half cystine"
        // but we can't see spaces in the tokenizer - they're not tokenized
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

    if (prefix == "UNIMOD")
    {
      cv.database = CvDatabase::UNIMOD;
    }
    else if (prefix == "MOD")
    {
      cv.database = CvDatabase::MOD;
    }
    else if (prefix == "RESID")
    {
      cv.database = CvDatabase::RESID;
    }
    else if (prefix == "XLMOD")
    {
      cv.database = CvDatabase::XLMOD;
    }
    else if (prefix == "GNO")
    {
      cv.database = CvDatabase::GNO;
    }
    else
    {
      errorAt_(ProFormaErrorCode::INVALID_CV_PREFIX, prefix_tok.position,
               "Invalid CV prefix");
    }

    expect_(ProFormaTokenizer::TokenType::COLON, "':' after CV prefix");

    // GNO and RESID accessions can be alphanumeric (e.g., G59626AS, AA0581)
    // UNIMOD, MOD, XLMOD accessions are numeric
    if (cv.database == CvDatabase::GNO || cv.database == CvDatabase::RESID)
    {
      // Parse alphanumeric accession (letters and numbers)
      std::string accession;
      while (!isAtEnd_())
      {
        ProFormaTokenizer::Token tok = current_();
        if (tok.type == ProFormaTokenizer::TokenType::IDENTIFIER ||
            tok.type == ProFormaTokenizer::TokenType::NUMBER)
        {
          accession += std::string(tok.text);
          advance_();
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

    // Parse formula (e.g., Na, H, K)
    ProFormaTokenizer::Token formula = expect_(ProFormaTokenizer::TokenType::IDENTIFIER, "adduct formula");
    adduct.formula = String(formula.text);

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
    // Standard amino acid one-letter codes
    // A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    // Plus: U (selenocysteine), O (pyrrolysine), X (unknown), B (Asx), Z (Glx), J (Xle)
    return (c >= 'A' && c <= 'Z');
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
            return mod_db->getModification(full_accession);
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
    if (policy == AASequenceConversionPolicy::STRICT && !issues.empty())
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
          else if (policy == AASequenceConversionPolicy::STRICT)
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

} // namespace OpenMS
