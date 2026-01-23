// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <ostream>
#include <sstream>

namespace OpenMS
{

  //--------------------------------------------------------------------------
  // Error code to string conversion
  //--------------------------------------------------------------------------

  const char* mzPAFErrorCodeToString(MzPAFErrorCode code)
  {
    switch (code)
    {
      case MzPAFErrorCode::UNEXPECTED_CHARACTER: return "Unexpected character";
      case MzPAFErrorCode::UNCLOSED_BRACKET: return "Unclosed bracket";
      case MzPAFErrorCode::INVALID_ION_SERIES: return "Invalid ion series";
      case MzPAFErrorCode::INVALID_NUMBER: return "Invalid number";
      case MzPAFErrorCode::INVALID_FORMULA: return "Invalid formula";
      case MzPAFErrorCode::INVALID_CHARGE: return "Invalid charge";
      case MzPAFErrorCode::INVALID_DELTA: return "Invalid mass delta";
      case MzPAFErrorCode::INVALID_CONFIDENCE: return "Invalid confidence";
      case MzPAFErrorCode::EMPTY_INPUT: return "Empty input";
      case MzPAFErrorCode::UNEXPECTED_END_OF_INPUT: return "Unexpected end of input";
      case MzPAFErrorCode::INTERNAL_ERROR: return "Internal error";
      default: return "Unknown error";
    }
  }

  //--------------------------------------------------------------------------
  // MzPAFParseError implementation
  //--------------------------------------------------------------------------

  MzPAFParseError::MzPAFParseError(
    const char* file,
    int line,
    const char* function,
    MzPAFErrorCode error_code,
    size_t error_position,
    const String& input,
    const String& message
  ) noexcept :
    Exception::ParseError(file, line, function, input, message),
    code_(error_code),
    position_(error_position)
  {
    extractContext_(input, error_position);
  }

  void MzPAFParseError::extractContext_(const String& input, size_t pos)
  {
    if (pos > 20)
    {
      context_before_ = input.substr(pos - 20, 20);
    }
    else
    {
      context_before_ = input.substr(0, pos);
    }

    if (pos < input.size())
    {
      size_t remaining = input.size() - pos;
      context_after_ = input.substr(pos, std::min(remaining, size_t(20)));
    }
  }

  String MzPAFParseError::getFormattedMessage() const
  {
    std::ostringstream oss;
    oss << "mzPAF parse error at position " << position_ << ": "
        << mzPAFErrorCodeToString(code_) << "\n";
    oss << "Context: " << context_before_ << ">>>" << context_after_ << "<<<";
    return String(oss.str());
  }

  //--------------------------------------------------------------------------
  // MzPAFNeutralLoss implementation
  //--------------------------------------------------------------------------

  bool MzPAFNeutralLoss::operator==(const MzPAFNeutralLoss& other) const
  {
    // Compare only formula for semantic equality (original_text is for lossless roundtrip)
    return formula == other.formula;
  }

  //--------------------------------------------------------------------------
  // MzPAFMassDelta implementation
  //--------------------------------------------------------------------------

  bool MzPAFMassDelta::operator==(const MzPAFMassDelta& other) const
  {
    return std::abs(value - other.value) < 1e-9 && unit == other.unit;
  }

  //--------------------------------------------------------------------------
  // MzPAFAnnotation implementation
  //--------------------------------------------------------------------------

  bool MzPAFAnnotation::isValid() const
  {
    // Must have a known ion series (not UNKNOWN) to be valid
    if (ion_series == MzPAFIonSeries::UNKNOWN)
    {
      return false;
    }

    // Standard fragment ions need ordinal
    if (ion_series == MzPAFIonSeries::A || ion_series == MzPAFIonSeries::B ||
        ion_series == MzPAFIonSeries::C || ion_series == MzPAFIonSeries::X ||
        ion_series == MzPAFIonSeries::Y || ion_series == MzPAFIonSeries::Z)
    {
      if (!ordinal.has_value())
      {
        return false;
      }
    }

    // Immonium needs residue
    if (ion_series == MzPAFIonSeries::IMMONIUM)
    {
      if (!immonium_residue.has_value())
      {
        return false;
      }
    }

    // Internal needs range
    if (ion_series == MzPAFIonSeries::INTERNAL)
    {
      if (!internal_range.has_value())
      {
        return false;
      }
    }

    // Reporter needs name
    if (ion_series == MzPAFIonSeries::REPORTER)
    {
      if (!reporter_name.has_value())
      {
        return false;
      }
    }

    // Formula needs formula
    if (ion_series == MzPAFIonSeries::FORMULA)
    {
      if (!formula.has_value())
      {
        return false;
      }
    }

    // Named needs name
    if (ion_series == MzPAFIonSeries::NAMED)
    {
      if (!named_compound.has_value())
      {
        return false;
      }
    }

    return true;
  }

  bool MzPAFAnnotation::operator==(const MzPAFAnnotation& other) const
  {
    return analyte_index == other.analyte_index &&
           ion_series == other.ion_series &&
           ordinal == other.ordinal &&
           immonium_residue == other.immonium_residue &&
           internal_range == other.internal_range &&
           reporter_name == other.reporter_name &&
           formula == other.formula &&
           named_compound == other.named_compound &&
           neutral_losses == other.neutral_losses &&
           isotope_offset == other.isotope_offset &&
           adduct == other.adduct &&
           charge == other.charge &&
           mass_delta == other.mass_delta &&
           confidence == other.confidence &&
           embedded_sequence == other.embedded_sequence;
  }

  std::ostream& operator<<(std::ostream& os, const MzPAFAnnotation& ann)
  {
    os << MzPAF::toString(ann, MzPAFWriteMode::CANONICAL);
    return os;
  }

  //--------------------------------------------------------------------------
  // MzPAFPeakAnnotations implementation
  //--------------------------------------------------------------------------

  bool MzPAFPeakAnnotations::operator==(const MzPAFPeakAnnotations& other) const
  {
    return annotations == other.annotations;
  }

  //--------------------------------------------------------------------------
  // Internal tokenizer
  //--------------------------------------------------------------------------

  namespace
  {
    enum class TokenType
    {
      NUMBER,
      IDENTIFIER,
      LBRACKET,   // [
      RBRACKET,   // ]
      LBRACE,     // {
      RBRACE,     // }
      PLUS,       // +
      MINUS,      // -
      CARET,      // ^
      SLASH,      // /
      ASTERISK,   // *
      COLON,      // :
      COMMA,      // ,
      AT,         // @
      END
    };

    struct Token
    {
      TokenType type;
      std::string_view text;
      size_t position;
    };

    class Tokenizer
    {
    public:
      explicit Tokenizer(std::string_view input) : input_(input), pos_(0) {}

      Token next()
      {
        skipWhitespace_();

        if (isAtEnd_())
        {
          return {TokenType::END, {}, pos_};
        }

        size_t start = pos_;
        char c = advance_();

        switch (c)
        {
          case '[': return {TokenType::LBRACKET, input_.substr(start, 1), start};
          case ']': return {TokenType::RBRACKET, input_.substr(start, 1), start};
          case '{': return {TokenType::LBRACE, input_.substr(start, 1), start};
          case '}': return {TokenType::RBRACE, input_.substr(start, 1), start};
          case '+': return scanNumberOrPlus_(start);
          case '-': return scanNumberOrMinus_(start);
          case '^': return {TokenType::CARET, input_.substr(start, 1), start};
          case '/': return {TokenType::SLASH, input_.substr(start, 1), start};
          case '*': return {TokenType::ASTERISK, input_.substr(start, 1), start};
          case ':': return {TokenType::COLON, input_.substr(start, 1), start};
          case ',': return {TokenType::COMMA, input_.substr(start, 1), start};
          case '@': return {TokenType::AT, input_.substr(start, 1), start};
          default:
            if (std::isdigit(c))
            {
              return scanNumber_(start);
            }
            if (std::isalpha(c) || c == '_')
            {
              return scanIdentifier_(start);
            }
            // Unknown character - return as identifier
            return {TokenType::IDENTIFIER, input_.substr(start, 1), start};
        }
      }

      Token peek()
      {
        size_t saved_pos = pos_;
        Token tok = next();
        pos_ = saved_pos;
        return tok;
      }

      bool hasMore() const { return pos_ < input_.size(); }
      size_t position() const { return pos_; }

    private:
      void skipWhitespace_()
      {
        while (!isAtEnd_() && std::isspace(static_cast<unsigned char>(input_[pos_])))
        {
          pos_++;
        }
      }

      bool isAtEnd_() const { return pos_ >= input_.size(); }

      char current_() const
      {
        return isAtEnd_() ? '\0' : input_[pos_];
      }

      char advance_()
      {
        return input_[pos_++];
      }

      Token scanNumber_(size_t start)
      {
        // Already consumed first digit
        while (!isAtEnd_() && std::isdigit(current_()))
        {
          advance_();
        }
        // Check for decimal
        if (!isAtEnd_() && current_() == '.')
        {
          advance_();
          while (!isAtEnd_() && std::isdigit(current_()))
          {
            advance_();
          }
        }
        return {TokenType::NUMBER, input_.substr(start, pos_ - start), start};
      }

      Token scanNumberOrPlus_(size_t start)
      {
        // In mzPAF, '+' is always a separate token (for isotope offset, adducts)
        // Don't combine with following digits
        return {TokenType::PLUS, input_.substr(start, 1), start};
      }

      Token scanNumberOrMinus_(size_t start)
      {
        // In mzPAF, '-' is always a separate token (for neutral losses)
        // Don't combine with following digits
        return {TokenType::MINUS, input_.substr(start, 1), start};
      }

      Token scanIdentifier_(size_t start)
      {
        // Already consumed first char
        while (!isAtEnd_() && (std::isalnum(current_()) || current_() == '_'))
        {
          advance_();
        }
        return {TokenType::IDENTIFIER, input_.substr(start, pos_ - start), start};
      }

      std::string_view input_;
      size_t pos_;
    };

    // Parser helper class
    class Parser
    {
    public:
      explicit Parser(const String& input) :
        input_(input), tokenizer_(input)
      {
        advance_();  // Load first token
      }

      MzPAFPeakAnnotations parseAll()
      {
        MzPAFPeakAnnotations result;

        if (current_.type == TokenType::END)
        {
          return result;  // Empty input
        }

        result.annotations.push_back(parseAnnotation_());

        // Parse comma-separated additional annotations
        while (match_(TokenType::COMMA))
        {
          result.annotations.push_back(parseAnnotation_());
        }

        return result;
      }

    private:
      MzPAFAnnotation parseAnnotation_()
      {
        MzPAFAnnotation ann;

        // Check for analyte index: N@
        if (current_.type == TokenType::NUMBER)
        {
          Token num = current_;
          advance_();
          if (current_.type == TokenType::AT)
          {
            advance_();
            try
            {
              ann.analyte_index = String(num.text).toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid analyte index");
            }
          }
          else
          {
            // Not analyte index, must be something else - error
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected ion series after number");
          }
        }

        // Parse ion series
        parseIonSeries_(ann);

        // Parse optional modifiers: neutral losses, isotope, charge, delta, confidence
        parseModifiers_(ann);

        return ann;
      }

      void parseIonSeries_(MzPAFAnnotation& ann)
      {
        if (current_.type != TokenType::IDENTIFIER)
        {
          error_(MzPAFErrorCode::INVALID_ION_SERIES, "Expected ion series identifier");
        }

        std::string_view text = current_.text;

        if (text.empty())
        {
          error_(MzPAFErrorCode::INVALID_ION_SERIES, "Empty identifier");
        }

        char first = text[0];

        // Standard fragment ions: a, b, c, x, y, z
        if (first == 'a' || first == 'b' || first == 'c' ||
            first == 'x' || first == 'y' || first == 'z')
        {
          MzPAF::charToIonSeries(first, ann.ion_series);

          // Check if ordinal is embedded in the identifier (e.g., "y4" as single token)
          if (text.size() > 1 && std::isdigit(text[1]))
          {
            // Parse ordinal from the identifier itself
            try
            {
              ann.ordinal = String(text.substr(1)).toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid ordinal number");
            }
            advance_();
          }
          else
          {
            advance_();
            // Parse ordinal from next token
            if (current_.type == TokenType::NUMBER)
            {
              try
              {
                ann.ordinal = String(current_.text).toInt();
              }
              catch (const Exception::ConversionError&)
              {
                error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid ordinal number");
              }
              advance_();
            }
            else
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Expected ordinal number after ion series");
            }
          }

          // Check for embedded sequence: {sequence}
          if (current_.type == TokenType::LBRACE)
          {
            parseEmbeddedSequence_(ann);
          }
        }
        // Precursor: p or pN
        else if (first == 'p')
        {
          ann.ion_series = MzPAFIonSeries::PRECURSOR;
          advance_();

          // Optional ordinal for precursor
          if (current_.type == TokenType::NUMBER)
          {
            try
            {
              ann.ordinal = String(current_.text).toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid precursor ordinal");
            }
            advance_();
          }
        }
        // Immonium: I followed by single letter
        else if (first == 'I' && text.size() >= 2)
        {
          ann.ion_series = MzPAFIonSeries::IMMONIUM;
          ann.immonium_residue = text[1];
          advance_();
        }
        // Internal fragment: mN:M
        else if (first == 'm')
        {
          ann.ion_series = MzPAFIonSeries::INTERNAL;

          int start_pos;
          // Check if start position is embedded in the identifier (e.g., "m3" as single token)
          if (text.size() > 1 && std::isdigit(text[1]))
          {
            try
            {
              start_pos = String(text.substr(1)).toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment start position");
            }
            advance_();
          }
          else
          {
            advance_();
            // Parse start position from next token
            if (current_.type != TokenType::NUMBER)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Expected start position for internal fragment");
            }
            try
            {
              start_pos = String(current_.text).toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment start position");
            }
            advance_();
          }

          // Expect colon
          if (current_.type != TokenType::COLON)
          {
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected ':' in internal fragment");
          }
          advance_();

          // Parse end position
          if (current_.type != TokenType::NUMBER)
          {
            error_(MzPAFErrorCode::INVALID_NUMBER, "Expected end position for internal fragment");
          }
          int end_pos;
          try
          {
            end_pos = String(current_.text).toInt();
          }
          catch (const Exception::ConversionError&)
          {
            error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment end position");
          }
          advance_();

          ann.internal_range = std::make_pair(start_pos, end_pos);
        }
        // Reporter ion: r[name]
        else if (first == 'r')
        {
          ann.ion_series = MzPAFIonSeries::REPORTER;
          advance_();

          if (current_.type != TokenType::LBRACKET)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Expected '[' after 'r' for reporter ion");
          }
          advance_();

          // Parse reporter name
          String name;
          while (current_.type != TokenType::RBRACKET && current_.type != TokenType::END)
          {
            name += String(current_.text);
            advance_();
          }

          if (current_.type != TokenType::RBRACKET)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Unclosed bracket in reporter ion");
          }
          advance_();

          ann.reporter_name = name;
        }
        // Formula ion: f{formula}
        else if (first == 'f')
        {
          ann.ion_series = MzPAFIonSeries::FORMULA;
          advance_();

          if (current_.type != TokenType::LBRACE)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Expected '{' after 'f' for formula ion");
          }
          advance_();

          // Parse formula
          String formula_str;
          while (current_.type != TokenType::RBRACE && current_.type != TokenType::END)
          {
            formula_str += String(current_.text);
            advance_();
          }

          if (current_.type != TokenType::RBRACE)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Unclosed brace in formula ion");
          }
          advance_();

          try
          {
            ann.formula = EmpiricalFormula(formula_str);
          }
          catch (...)
          {
            error_(MzPAFErrorCode::INVALID_FORMULA, "Invalid chemical formula");
          }
        }
        // Named compound: _[name]
        else if (first == '_')
        {
          ann.ion_series = MzPAFIonSeries::NAMED;
          advance_();

          if (current_.type != TokenType::LBRACKET)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Expected '[' after '_' for named compound");
          }
          advance_();

          // Parse compound name
          String name;
          while (current_.type != TokenType::RBRACKET && current_.type != TokenType::END)
          {
            name += String(current_.text);
            advance_();
          }

          if (current_.type != TokenType::RBRACKET)
          {
            error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Unclosed bracket in named compound");
          }
          advance_();

          ann.named_compound = name;
        }
        else
        {
          error_(MzPAFErrorCode::INVALID_ION_SERIES, "Unrecognized ion series");
        }
      }

      void parseEmbeddedSequence_(MzPAFAnnotation& ann)
      {
        // Current token is LBRACE
        advance_();

        String seq_str;
        int brace_depth = 1;

        while (brace_depth > 0 && current_.type != TokenType::END)
        {
          if (current_.type == TokenType::LBRACE)
          {
            brace_depth++;
          }
          else if (current_.type == TokenType::RBRACE)
          {
            brace_depth--;
            if (brace_depth == 0)
            {
              break;
            }
          }
          seq_str += String(current_.text);
          advance_();
        }

        if (current_.type != TokenType::RBRACE)
        {
          error_(MzPAFErrorCode::UNCLOSED_BRACKET, "Unclosed brace in embedded sequence");
        }
        advance_();

        ann.embedded_sequence = seq_str;
      }

      void parseModifiers_(MzPAFAnnotation& ann)
      {
        // Parse in a loop since modifiers can appear in various orders
        while (true)
        {
          // Neutral loss: -formula
          if (current_.type == TokenType::MINUS)
          {
            parseNeutralLoss_(ann);
          }
          // Isotope offset: +Ni
          else if (current_.type == TokenType::PLUS)
          {
            parseIsotopeOrAdduct_(ann);
          }
          // Charge: ^N
          else if (current_.type == TokenType::CARET)
          {
            parseCharge_(ann);
          }
          // Mass delta: /value[ppm]
          else if (current_.type == TokenType::SLASH)
          {
            parseMassDelta_(ann);
          }
          // Confidence: *value
          else if (current_.type == TokenType::ASTERISK)
          {
            parseConfidence_(ann);
          }
          else
          {
            break;  // No more modifiers
          }
        }
      }

      void parseNeutralLoss_(MzPAFAnnotation& ann)
      {
        // Current token is MINUS
        advance_();

        // Parse the loss formula (identifier)
        if (current_.type != TokenType::IDENTIFIER)
        {
          error_(MzPAFErrorCode::INVALID_FORMULA, "Expected formula after '-'");
        }

        MzPAFNeutralLoss loss;
        loss.original_text = String(current_.text);

        try
        {
          loss.formula = EmpiricalFormula(loss.original_text);
        }
        catch (...)
        {
          error_(MzPAFErrorCode::INVALID_FORMULA, "Invalid neutral loss formula");
        }

        advance_();
        ann.neutral_losses.push_back(loss);
      }

      void parseIsotopeOrAdduct_(MzPAFAnnotation& ann)
      {
        // Current token is PLUS
        advance_();

        if (current_.type == TokenType::NUMBER)
        {
          // Could be isotope offset: +Ni
          String num_str(current_.text);
          advance_();

          // Check for 'i' suffix
          if (current_.type == TokenType::IDENTIFIER &&
              current_.text.size() >= 1 && current_.text[0] == 'i')
          {
            try
            {
              ann.isotope_offset = num_str.toInt();
            }
            catch (const Exception::ConversionError&)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Invalid isotope offset");
            }
            // Skip past 'i' - it might be just 'i' or 'i' followed by more
            if (current_.text.size() == 1)
            {
              advance_();
            }
            // If it's like "i" followed by something else, we only consumed the identifier
          }
          else
          {
            // Not isotope, might be something else - treat as error for now
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected 'i' after number for isotope offset");
          }
        }
        else if (current_.type == TokenType::IDENTIFIER)
        {
          // Adduct: +Na, +K, etc.
          try
          {
            ann.adduct = EmpiricalFormula(String(current_.text));
          }
          catch (...)
          {
            error_(MzPAFErrorCode::INVALID_FORMULA, "Invalid adduct formula");
          }
          advance_();
        }
        else
        {
          error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected number or formula after '+'");
        }
      }

      void parseCharge_(MzPAFAnnotation& ann)
      {
        // Current token is CARET
        advance_();

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_CHARGE, "Expected charge number after '^'");
        }

        try
        {
          ann.charge = String(current_.text).toInt();
        }
        catch (const Exception::ConversionError&)
        {
          error_(MzPAFErrorCode::INVALID_CHARGE, "Invalid charge number");
        }
        advance_();
      }

      void parseMassDelta_(MzPAFAnnotation& ann)
      {
        // Current token is SLASH
        advance_();

        MzPAFMassDelta delta;
        delta.unit = MzPAFDeltaUnit::DALTON;

        // Handle optional sign
        double sign = 1.0;
        String sign_str;
        if (current_.type == TokenType::MINUS)
        {
          sign = -1.0;
          sign_str = "-";
          advance_();
        }
        else if (current_.type == TokenType::PLUS)
        {
          sign_str = "+";
          advance_();
        }

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_DELTA, "Expected number after '/'");
        }

        String num_str(current_.text);
        try
        {
          delta.value = sign * num_str.toDouble();
        }
        catch (const Exception::ConversionError&)
        {
          error_(MzPAFErrorCode::INVALID_DELTA, "Invalid mass delta value");
        }
        delta.original_text = sign_str + num_str;
        advance_();

        // Check for 'ppm' suffix
        if (current_.type == TokenType::IDENTIFIER)
        {
          String suffix(current_.text);
          if (suffix == "ppm")
          {
            delta.unit = MzPAFDeltaUnit::PPM;
            delta.original_text += "ppm";
            advance_();
          }
        }

        ann.mass_delta = delta;
      }

      void parseConfidence_(MzPAFAnnotation& ann)
      {
        // Current token is ASTERISK
        advance_();

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_CONFIDENCE, "Expected confidence value after '*'");
        }

        try
        {
          ann.confidence = String(current_.text).toDouble();
        }
        catch (const Exception::ConversionError&)
        {
          error_(MzPAFErrorCode::INVALID_CONFIDENCE, "Invalid confidence value");
        }
        advance_();
      }

      void advance_()
      {
        current_ = tokenizer_.next();
      }

      bool match_(TokenType type)
      {
        if (current_.type == type)
        {
          advance_();
          return true;
        }
        return false;
      }

      [[noreturn]] void error_(MzPAFErrorCode code, const char* message)
      {
        throw MzPAFParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                              code, current_.position, input_, message);
      }

      String input_;
      Tokenizer tokenizer_;
      Token current_;
    };

  } // anonymous namespace

  //--------------------------------------------------------------------------
  // MzPAF implementation
  //--------------------------------------------------------------------------

  MzPAFAnnotation MzPAF::parse(const String& input)
  {
    if (input.empty())
    {
      throw MzPAFParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           MzPAFErrorCode::EMPTY_INPUT, 0, input, "Empty input");
    }

    Parser parser(input);
    MzPAFPeakAnnotations result = parser.parseAll();

    if (result.empty())
    {
      throw MzPAFParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           MzPAFErrorCode::EMPTY_INPUT, 0, input, "No annotations parsed");
    }

    return result.annotations[0];
  }

  MzPAFPeakAnnotations MzPAF::parseMultiple(const String& input)
  {
    if (input.empty())
    {
      throw MzPAFParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           MzPAFErrorCode::EMPTY_INPUT, 0, input, "Empty input");
    }

    Parser parser(input);
    return parser.parseAll();
  }

  std::optional<MzPAFAnnotation> MzPAF::tryParse(const String& input)
  {
    try
    {
      return parse(input);
    }
    catch (...)
    {
      return std::nullopt;
    }
  }

  std::optional<MzPAFPeakAnnotations> MzPAF::tryParseMultiple(const String& input)
  {
    try
    {
      return parseMultiple(input);
    }
    catch (...)
    {
      return std::nullopt;
    }
  }

  //--------------------------------------------------------------------------
  // Writer implementation
  //--------------------------------------------------------------------------

  String MzPAF::toString(const MzPAFAnnotation& ann, MzPAFWriteMode mode)
  {
    std::ostringstream oss;

    // Analyte index
    if (ann.analyte_index.has_value())
    {
      oss << ann.analyte_index.value() << "@";
    }

    // Ion series and specifics
    switch (ann.ion_series)
    {
      case MzPAFIonSeries::A:
      case MzPAFIonSeries::B:
      case MzPAFIonSeries::C:
      case MzPAFIonSeries::X:
      case MzPAFIonSeries::Y:
      case MzPAFIonSeries::Z:
        oss << ionSeriesToChar(ann.ion_series);
        if (ann.ordinal.has_value())
        {
          oss << ann.ordinal.value();
        }
        if (ann.embedded_sequence.has_value())
        {
          oss << "{" << ann.embedded_sequence.value() << "}";
        }
        break;

      case MzPAFIonSeries::PRECURSOR:
        oss << "p";
        if (ann.ordinal.has_value())
        {
          oss << ann.ordinal.value();
        }
        break;

      case MzPAFIonSeries::IMMONIUM:
        oss << "I";
        if (ann.immonium_residue.has_value())
        {
          oss << ann.immonium_residue.value();
        }
        break;

      case MzPAFIonSeries::INTERNAL:
        oss << "m";
        if (ann.internal_range.has_value())
        {
          oss << ann.internal_range.value().first << ":"
              << ann.internal_range.value().second;
        }
        break;

      case MzPAFIonSeries::REPORTER:
        oss << "r[";
        if (ann.reporter_name.has_value())
        {
          oss << ann.reporter_name.value();
        }
        oss << "]";
        break;

      case MzPAFIonSeries::FORMULA:
        oss << "f{";
        if (ann.formula.has_value())
        {
          oss << ann.formula.value().toString();
        }
        oss << "}";
        break;

      case MzPAFIonSeries::NAMED:
        oss << "_[";
        if (ann.named_compound.has_value())
        {
          oss << ann.named_compound.value();
        }
        oss << "]";
        break;

      case MzPAFIonSeries::UNKNOWN:
        oss << "?";
        break;
    }

    // Neutral losses
    for (const auto& loss : ann.neutral_losses)
    {
      oss << "-";
      if (mode == MzPAFWriteMode::LOSSLESS && !loss.original_text.empty())
      {
        oss << loss.original_text;
      }
      else
      {
        oss << loss.formula.toString();
      }
    }

    // Isotope offset
    if (ann.isotope_offset.has_value())
    {
      oss << "+" << ann.isotope_offset.value() << "i";
    }

    // Adduct
    if (ann.adduct.has_value())
    {
      oss << "+" << ann.adduct.value().toString();
    }

    // Charge
    if (ann.charge.has_value())
    {
      oss << "^" << ann.charge.value();
    }

    // Mass delta
    if (ann.mass_delta.has_value())
    {
      oss << "/";
      if (mode == MzPAFWriteMode::LOSSLESS && !ann.mass_delta.value().original_text.empty())
      {
        oss << ann.mass_delta.value().original_text;
      }
      else
      {
        oss << ann.mass_delta.value().value;
        if (ann.mass_delta.value().unit == MzPAFDeltaUnit::PPM)
        {
          oss << "ppm";
        }
      }
    }

    // Confidence
    if (ann.confidence.has_value())
    {
      oss << "*" << ann.confidence.value();
    }

    return String(oss.str());
  }

  String MzPAF::toString(const MzPAFPeakAnnotations& anns, MzPAFWriteMode mode)
  {
    if (anns.empty())
    {
      return "";
    }

    std::ostringstream oss;
    for (size_t i = 0; i < anns.annotations.size(); ++i)
    {
      if (i > 0)
      {
        oss << ",";
      }
      oss << toString(anns.annotations[i], mode);
    }

    return String(oss.str());
  }

  //--------------------------------------------------------------------------
  // PeakAnnotation integration
  //--------------------------------------------------------------------------

  PeptideHit::PeakAnnotation MzPAF::toPeakAnnotation(
    const MzPAFAnnotation& mzpaf, double mz, double intensity)
  {
    PeptideHit::PeakAnnotation pa;
    pa.annotation = toString(mzpaf, MzPAFWriteMode::CANONICAL);
    pa.charge = mzpaf.charge.value_or(1);
    pa.mz = mz;
    pa.intensity = intensity;
    return pa;
  }

  MzPAFPeakAnnotations MzPAF::fromPeakAnnotation(
    const PeptideHit::PeakAnnotation& peak_annotation)
  {
    auto result = tryParseMultiple(peak_annotation.annotation);
    return result.value_or(MzPAFPeakAnnotations{});
  }

  //--------------------------------------------------------------------------
  // Utilities
  //--------------------------------------------------------------------------

  bool MzPAF::isMzPAFFormat(const String& annotation)
  {
    if (annotation.empty())
    {
      return false;
    }

    // Quick heuristic: should start with ion series letter or digit (for analyte)
    char first = annotation[0];

    // Valid starting characters for mzPAF
    if (first == 'a' || first == 'b' || first == 'c' ||
        first == 'x' || first == 'y' || first == 'z' ||
        first == 'p' || first == 'I' || first == 'm' ||
        first == 'r' || first == 'f' || first == '_' ||
        std::isdigit(first))
    {
      // Try parsing
      return tryParse(annotation).has_value();
    }

    return false;
  }

  std::optional<double> MzPAF::calculateTheoreticalMZ(
    const MzPAFAnnotation& ann, const AASequence& sequence)
  {
    // Only handle standard fragment ions for now
    if (ann.ion_series != MzPAFIonSeries::A &&
        ann.ion_series != MzPAFIonSeries::B &&
        ann.ion_series != MzPAFIonSeries::C &&
        ann.ion_series != MzPAFIonSeries::X &&
        ann.ion_series != MzPAFIonSeries::Y &&
        ann.ion_series != MzPAFIonSeries::Z)
    {
      return std::nullopt;
    }

    if (!ann.ordinal.has_value())
    {
      return std::nullopt;
    }

    int pos = ann.ordinal.value();
    int charge = ann.charge.value_or(1);

    if (pos < 1 || pos >= static_cast<int>(sequence.size()))
    {
      return std::nullopt;
    }

    double mass = 0.0;

    // Calculate fragment mass based on ion type
    switch (ann.ion_series)
    {
      case MzPAFIonSeries::A:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::Internal) - 27.9949;  // -CO
        break;
      case MzPAFIonSeries::B:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::Internal);
        break;
      case MzPAFIonSeries::C:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::Internal) + 17.0265;  // +NH3
        break;
      case MzPAFIonSeries::X:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::Internal) + 25.9793;  // +CO-H
        break;
      case MzPAFIonSeries::Y:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::Internal) + 18.0106;  // +H2O
        break;
      case MzPAFIonSeries::Z:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::Internal) + 1.9919;   // +H2O-NH3
        break;
      default:
        return std::nullopt;
    }

    // Apply neutral losses
    for (const auto& loss : ann.neutral_losses)
    {
      mass -= loss.formula.getMonoWeight();
    }

    // Apply isotope offset
    if (ann.isotope_offset.has_value())
    {
      mass += ann.isotope_offset.value() * 1.00335;  // Approximate mass difference per isotope
    }

    // Calculate m/z
    double mz = (mass + charge * Constants::PROTON_MASS_U) / charge;

    return mz;
  }

  char MzPAF::ionSeriesToChar(MzPAFIonSeries series)
  {
    switch (series)
    {
      case MzPAFIonSeries::A: return 'a';
      case MzPAFIonSeries::B: return 'b';
      case MzPAFIonSeries::C: return 'c';
      case MzPAFIonSeries::X: return 'x';
      case MzPAFIonSeries::Y: return 'y';
      case MzPAFIonSeries::Z: return 'z';
      case MzPAFIonSeries::PRECURSOR: return 'p';
      case MzPAFIonSeries::IMMONIUM: return 'I';
      case MzPAFIonSeries::INTERNAL: return 'm';
      case MzPAFIonSeries::REPORTER: return 'r';
      case MzPAFIonSeries::FORMULA: return 'f';
      case MzPAFIonSeries::NAMED: return '_';
      case MzPAFIonSeries::UNKNOWN: return '?';
      default: return '?';
    }
  }

  bool MzPAF::charToIonSeries(char c, MzPAFIonSeries& series)
  {
    switch (c)
    {
      case 'a': series = MzPAFIonSeries::A; return true;
      case 'b': series = MzPAFIonSeries::B; return true;
      case 'c': series = MzPAFIonSeries::C; return true;
      case 'x': series = MzPAFIonSeries::X; return true;
      case 'y': series = MzPAFIonSeries::Y; return true;
      case 'z': series = MzPAFIonSeries::Z; return true;
      case 'p': series = MzPAFIonSeries::PRECURSOR; return true;
      case 'I': series = MzPAFIonSeries::IMMONIUM; return true;
      case 'm': series = MzPAFIonSeries::INTERNAL; return true;
      case 'r': series = MzPAFIonSeries::REPORTER; return true;
      case 'f': series = MzPAFIonSeries::FORMULA; return true;
      case '_': series = MzPAFIonSeries::NAMED; return true;
      default: return false;
    }
  }

} // namespace OpenMS
