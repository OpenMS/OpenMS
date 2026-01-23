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
    if (ion_series == MzPAFIonSeries::UNKNOWN)
    {
      return false;
    }

    // Standard fragment ions need ordinal
    if (MzPAF::isStandardFragmentIon(ion_series) && !ordinal.has_value())
    {
      return false;
    }

    // Immonium needs residue
    if (ion_series == MzPAFIonSeries::IMMONIUM && !immonium_residue.has_value())
    {
      return false;
    }

    // Internal needs range
    if (ion_series == MzPAFIonSeries::INTERNAL && !internal_range.has_value())
    {
      return false;
    }

    // Reporter needs name
    if (ion_series == MzPAFIonSeries::REPORTER && !reporter_name.has_value())
    {
      return false;
    }

    // Formula needs formula
    if (ion_series == MzPAFIonSeries::FORMULA && !formula.has_value())
    {
      return false;
    }

    // Named needs name
    if (ion_series == MzPAFIonSeries::NAMED && !named_compound.has_value())
    {
      return false;
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
    os << MzPAF::toString(ann);
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
          case '+': return {TokenType::PLUS, input_.substr(start, 1), start};
          case '-': return {TokenType::MINUS, input_.substr(start, 1), start};
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
            return {TokenType::IDENTIFIER, input_.substr(start, 1), start};
        }
      }

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
        while (!isAtEnd_() && std::isdigit(current_()))
        {
          advance_();
        }
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

      Token scanIdentifier_(size_t start)
      {
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
        advance_();
      }

      MzPAFPeakAnnotations parseAll()
      {
        MzPAFPeakAnnotations result;

        if (current_.type == TokenType::END)
        {
          return result;
        }

        result.annotations.push_back(parseAnnotation_());

        while (match_(TokenType::COMMA))
        {
          result.annotations.push_back(parseAnnotation_());
        }

        return result;
      }

    private:
      //------------------------------------------------------------------------
      // Helper methods (DRY)
      //------------------------------------------------------------------------

      int parseInt_(std::string_view text, MzPAFErrorCode code, const char* msg)
      {
        try
        {
          return String(text).toInt();
        }
        catch (const Exception::ConversionError&)
        {
          error_(code, msg);
        }
      }

      double parseDouble_(std::string_view text, MzPAFErrorCode code, const char* msg)
      {
        try
        {
          return String(text).toDouble();
        }
        catch (const Exception::ConversionError&)
        {
          error_(code, msg);
        }
      }

      String parseBracketedContent_(TokenType open, TokenType close, const char* open_err, const char* close_err)
      {
        if (current_.type != open)
        {
          error_(MzPAFErrorCode::UNCLOSED_BRACKET, open_err);
        }
        advance_();

        String content;
        while (current_.type != close && current_.type != TokenType::END)
        {
          content += String(current_.text);
          advance_();
        }

        if (current_.type != close)
        {
          error_(MzPAFErrorCode::UNCLOSED_BRACKET, close_err);
        }
        advance_();

        return content;
      }

      //------------------------------------------------------------------------
      // Parsing methods
      //------------------------------------------------------------------------

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
            ann.analyte_index = parseInt_(num.text, MzPAFErrorCode::INVALID_NUMBER, "Invalid analyte index");
          }
          else
          {
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected ion series after number");
          }
        }

        parseIonSeries_(ann);
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

          if (text.size() > 1 && std::isdigit(text[1]))
          {
            ann.ordinal = parseInt_(text.substr(1), MzPAFErrorCode::INVALID_NUMBER, "Invalid ordinal number");
            advance_();
          }
          else
          {
            advance_();
            if (current_.type == TokenType::NUMBER)
            {
              ann.ordinal = parseInt_(current_.text, MzPAFErrorCode::INVALID_NUMBER, "Invalid ordinal number");
              advance_();
            }
            else
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Expected ordinal number after ion series");
            }
          }

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

          if (current_.type == TokenType::NUMBER)
          {
            ann.ordinal = parseInt_(current_.text, MzPAFErrorCode::INVALID_NUMBER, "Invalid precursor ordinal");
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
          if (text.size() > 1 && std::isdigit(text[1]))
          {
            start_pos = parseInt_(text.substr(1), MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment start");
            advance_();
          }
          else
          {
            advance_();
            if (current_.type != TokenType::NUMBER)
            {
              error_(MzPAFErrorCode::INVALID_NUMBER, "Expected start position for internal fragment");
            }
            start_pos = parseInt_(current_.text, MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment start");
            advance_();
          }

          if (current_.type != TokenType::COLON)
          {
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected ':' in internal fragment");
          }
          advance_();

          if (current_.type != TokenType::NUMBER)
          {
            error_(MzPAFErrorCode::INVALID_NUMBER, "Expected end position for internal fragment");
          }
          int end_pos = parseInt_(current_.text, MzPAFErrorCode::INVALID_NUMBER, "Invalid internal fragment end");
          advance_();

          ann.internal_range = std::make_pair(start_pos, end_pos);
        }
        // Reporter ion: r[name]
        else if (first == 'r')
        {
          ann.ion_series = MzPAFIonSeries::REPORTER;
          advance_();
          ann.reporter_name = parseBracketedContent_(TokenType::LBRACKET, TokenType::RBRACKET,
                                                     "Expected '[' after 'r'", "Unclosed bracket in reporter ion");
        }
        // Formula ion: f{formula}
        else if (first == 'f')
        {
          ann.ion_series = MzPAFIonSeries::FORMULA;
          advance_();
          String formula_str = parseBracketedContent_(TokenType::LBRACE, TokenType::RBRACE,
                                                      "Expected '{' after 'f'", "Unclosed brace in formula ion");
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
          ann.named_compound = parseBracketedContent_(TokenType::LBRACKET, TokenType::RBRACKET,
                                                      "Expected '[' after '_'", "Unclosed bracket in named compound");
        }
        else
        {
          error_(MzPAFErrorCode::INVALID_ION_SERIES, "Unrecognized ion series");
        }
      }

      void parseEmbeddedSequence_(MzPAFAnnotation& ann)
      {
        advance_(); // consume LBRACE

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
        while (true)
        {
          if (current_.type == TokenType::MINUS)
          {
            parseNeutralLoss_(ann);
          }
          else if (current_.type == TokenType::PLUS)
          {
            parseIsotopeOrAdduct_(ann);
          }
          else if (current_.type == TokenType::CARET)
          {
            parseCharge_(ann);
          }
          else if (current_.type == TokenType::SLASH)
          {
            parseMassDelta_(ann);
          }
          else if (current_.type == TokenType::ASTERISK)
          {
            parseConfidence_(ann);
          }
          else
          {
            break;
          }
        }
      }

      void parseNeutralLoss_(MzPAFAnnotation& ann)
      {
        advance_();

        if (current_.type != TokenType::IDENTIFIER)
        {
          error_(MzPAFErrorCode::INVALID_FORMULA, "Expected formula after '-'");
        }

        MzPAFNeutralLoss loss;
        try
        {
          loss.formula = EmpiricalFormula(String(current_.text));
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
        advance_();

        if (current_.type == TokenType::NUMBER)
        {
          String num_str(current_.text);
          advance_();

          if (current_.type == TokenType::IDENTIFIER &&
              !current_.text.empty() && current_.text[0] == 'i')
          {
            ann.isotope_offset = parseInt_(num_str, MzPAFErrorCode::INVALID_NUMBER, "Invalid isotope offset");
            if (current_.text.size() == 1)
            {
              advance_();
            }
          }
          else
          {
            error_(MzPAFErrorCode::UNEXPECTED_CHARACTER, "Expected 'i' after number for isotope offset");
          }
        }
        else if (current_.type == TokenType::IDENTIFIER)
        {
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
        advance_();

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_CHARGE, "Expected charge number after '^'");
        }

        ann.charge = parseInt_(current_.text, MzPAFErrorCode::INVALID_CHARGE, "Invalid charge number");
        advance_();
      }

      void parseMassDelta_(MzPAFAnnotation& ann)
      {
        advance_();

        MzPAFMassDelta delta;
        delta.unit = MzPAFDeltaUnit::DALTON;

        double sign = 1.0;
        if (current_.type == TokenType::MINUS)
        {
          sign = -1.0;
          advance_();
        }
        else if (current_.type == TokenType::PLUS)
        {
          advance_();
        }

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_DELTA, "Expected number after '/'");
        }

        delta.value = sign * parseDouble_(current_.text, MzPAFErrorCode::INVALID_DELTA, "Invalid mass delta value");
        advance_();

        if (current_.type == TokenType::IDENTIFIER)
        {
          String suffix(current_.text);
          if (suffix == "ppm")
          {
            delta.unit = MzPAFDeltaUnit::PPM;
            advance_();
          }
        }

        ann.mass_delta = delta;
      }

      void parseConfidence_(MzPAFAnnotation& ann)
      {
        advance_();

        if (current_.type != TokenType::NUMBER)
        {
          error_(MzPAFErrorCode::INVALID_CONFIDENCE, "Expected confidence value after '*'");
        }

        ann.confidence = parseDouble_(current_.text, MzPAFErrorCode::INVALID_CONFIDENCE, "Invalid confidence value");
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

  String MzPAF::toString(const MzPAFAnnotation& ann)
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
      oss << "-" << loss.formula.toString();
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
      oss << "/" << ann.mass_delta.value().value;
      if (ann.mass_delta.value().unit == MzPAFDeltaUnit::PPM)
      {
        oss << "ppm";
      }
    }

    // Confidence
    if (ann.confidence.has_value())
    {
      oss << "*" << ann.confidence.value();
    }

    return String(oss.str());
  }

  String MzPAF::toString(const MzPAFPeakAnnotations& anns)
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
      oss << toString(anns.annotations[i]);
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
    pa.annotation = toString(mzpaf);
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

  bool MzPAF::isStandardFragmentIon(MzPAFIonSeries series)
  {
    return series == MzPAFIonSeries::A || series == MzPAFIonSeries::B ||
           series == MzPAFIonSeries::C || series == MzPAFIonSeries::X ||
           series == MzPAFIonSeries::Y || series == MzPAFIonSeries::Z;
  }

  bool MzPAF::isMzPAFFormat(const String& annotation)
  {
    if (annotation.empty())
    {
      return false;
    }

    char first = annotation[0];

    if (first == 'a' || first == 'b' || first == 'c' ||
        first == 'x' || first == 'y' || first == 'z' ||
        first == 'p' || first == 'I' || first == 'm' ||
        first == 'r' || first == 'f' || first == '_' ||
        std::isdigit(first))
    {
      return tryParse(annotation).has_value();
    }

    return false;
  }

  std::optional<double> MzPAF::calculateTheoreticalMZ(
    const MzPAFAnnotation& ann, const AASequence& sequence)
  {
    if (!isStandardFragmentIon(ann.ion_series))
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

    // Use proper Residue types for accurate mass calculation
    switch (ann.ion_series)
    {
      case MzPAFIonSeries::A:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::AIon);
        break;
      case MzPAFIonSeries::B:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::BIon);
        break;
      case MzPAFIonSeries::C:
        mass = sequence.getPrefix(pos).getMonoWeight(Residue::CIon);
        break;
      case MzPAFIonSeries::X:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::XIon);
        break;
      case MzPAFIonSeries::Y:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::YIon);
        break;
      case MzPAFIonSeries::Z:
        mass = sequence.getSuffix(pos).getMonoWeight(Residue::ZIon);
        break;
      default:
        return std::nullopt;
    }

    // Subtract neutral losses
    for (const auto& loss : ann.neutral_losses)
    {
      mass -= loss.formula.getMonoWeight();
    }

    // Add isotope offset (C13-C12 mass difference per isotope)
    if (ann.isotope_offset.has_value())
    {
      mass += ann.isotope_offset.value() * Constants::C13C12_MASSDIFF_U;
    }

    // Calculate m/z: (mass + z * proton_mass) / z
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
