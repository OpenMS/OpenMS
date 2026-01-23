// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ProFormaDataJson.h>
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
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <type_traits>

namespace OpenMS
{

//============================================================================
// Internal ProFormaTokenizer class
//============================================================================
namespace detail
{

  /// Token types for ProForma lexical analysis
  enum class TokenType
  {
    LBRACKET,   // [
    RBRACKET,   // ]
    LPAREN,     // (
    RPAREN,     // )
    LBRACE,     // {
    RBRACE,     // }
    LANGLE,     // <
    RANGLE,     // >
    PLUS,       // +
    MINUS,      // -
    SLASH,      // /
    PIPE,       // |
    HASH,       // #
    COLON,      // :
    COMMA,      // ,
    CARET,      // ^
    QUESTION,   // ?
    AT,         // @
    NUMBER,     // Numeric literal (including decimals and signed values)
    IDENTIFIER, // Letter sequence
    END         // End of input
  };

  /// A single token from the ProForma input
  struct Token
  {
    TokenType type;
    std::string_view text;
    size_t position;
  };

  /// Zero-copy tokenizer for ProForma strings
  class ProFormaTokenizer
  {
  public:
    using TokenType = detail::TokenType;
    using Token = detail::Token;

    explicit ProFormaTokenizer(std::string_view input, size_t start_pos = 0)
      : input_(input),
        pos_(std::min(start_pos, input.size())),
        peeked_(std::nullopt)
    {
    }

    Token next()
    {
      if (peeked_.has_value())
      {
        Token token = *peeked_;
        peeked_.reset();
        return token;
      }
      return scanToken_();
    }

    Token peek()
    {
      if (peeked_.has_value())
      {
        return *peeked_;
      }
      peeked_ = scanToken_();
      return *peeked_;
    }

    bool hasMore() const
    {
      if (peeked_.has_value())
      {
        return peeked_->type != TokenType::END;
      }
      return pos_ < input_.size();
    }

    size_t position() const
    {
      return pos_;
    }

    std::string_view getContext(size_t pos, size_t before = 20, size_t after = 20) const
    {
      size_t start = (pos > before) ? (pos - before) : 0;
      size_t end = std::min(pos + after, input_.size());
      return input_.substr(start, end - start);
    }

    static const char* tokenTypeName(TokenType type)
    {
      switch (type)
      {
        case TokenType::LBRACKET:   return "LBRACKET";
        case TokenType::RBRACKET:   return "RBRACKET";
        case TokenType::LPAREN:     return "LPAREN";
        case TokenType::RPAREN:     return "RPAREN";
        case TokenType::LBRACE:     return "LBRACE";
        case TokenType::RBRACE:     return "RBRACE";
        case TokenType::LANGLE:     return "LANGLE";
        case TokenType::RANGLE:     return "RANGLE";
        case TokenType::PLUS:       return "PLUS";
        case TokenType::MINUS:      return "MINUS";
        case TokenType::SLASH:      return "SLASH";
        case TokenType::PIPE:       return "PIPE";
        case TokenType::HASH:       return "HASH";
        case TokenType::COLON:      return "COLON";
        case TokenType::COMMA:      return "COMMA";
        case TokenType::CARET:      return "CARET";
        case TokenType::QUESTION:   return "QUESTION";
        case TokenType::AT:         return "AT";
        case TokenType::NUMBER:     return "NUMBER";
        case TokenType::IDENTIFIER: return "IDENTIFIER";
        case TokenType::END:        return "END";
        default:                    return "UNKNOWN";
      }
    }

  private:
    std::string_view input_;
    size_t pos_;
    std::optional<Token> peeked_;

    Token scanToken_()
    {
      if (isAtEnd_()) return Token{TokenType::END, std::string_view{}, pos_};

      char c = current_();
      size_t start_pos = pos_;

      switch (c)
      {
        case '[': advance_(); return Token{TokenType::LBRACKET, input_.substr(start_pos, 1), start_pos};
        case ']': advance_(); return Token{TokenType::RBRACKET, input_.substr(start_pos, 1), start_pos};
        case '(': advance_(); return Token{TokenType::LPAREN, input_.substr(start_pos, 1), start_pos};
        case ')': advance_(); return Token{TokenType::RPAREN, input_.substr(start_pos, 1), start_pos};
        case '{': advance_(); return Token{TokenType::LBRACE, input_.substr(start_pos, 1), start_pos};
        case '}': advance_(); return Token{TokenType::RBRACE, input_.substr(start_pos, 1), start_pos};
        case '<': advance_(); return Token{TokenType::LANGLE, input_.substr(start_pos, 1), start_pos};
        case '>': advance_(); return Token{TokenType::RANGLE, input_.substr(start_pos, 1), start_pos};
        case '/': advance_(); return Token{TokenType::SLASH, input_.substr(start_pos, 1), start_pos};
        case '|': advance_(); return Token{TokenType::PIPE, input_.substr(start_pos, 1), start_pos};
        case '#': advance_(); return Token{TokenType::HASH, input_.substr(start_pos, 1), start_pos};
        case ':': advance_(); return Token{TokenType::COLON, input_.substr(start_pos, 1), start_pos};
        case ',': advance_(); return Token{TokenType::COMMA, input_.substr(start_pos, 1), start_pos};
        case '^': advance_(); return Token{TokenType::CARET, input_.substr(start_pos, 1), start_pos};
        case '?': advance_(); return Token{TokenType::QUESTION, input_.substr(start_pos, 1), start_pos};
        case '@': advance_(); return Token{TokenType::AT, input_.substr(start_pos, 1), start_pos};
        default: break;
      }

      if (c == '+' || c == '-')
      {
        char next_char = peek_(1);
        if (isDigit_(next_char) || next_char == '.')
        {
          return scanNumber_();
        }
        else
        {
          advance_();
          return Token{c == '+' ? TokenType::PLUS : TokenType::MINUS, input_.substr(start_pos, 1), start_pos};
        }
      }

      if (isDigit_(c) || (c == '.' && isDigit_(peek_(1))))
      {
        return scanNumber_();
      }

      if (isLetter_(c))
      {
        return scanIdentifier_();
      }

      advance_();
      return Token{TokenType::IDENTIFIER, input_.substr(start_pos, 1), start_pos};
    }

    Token scanNumber_()
    {
      size_t start_pos = pos_;
      if (current_() == '+' || current_() == '-') advance_();
      if (current_() == '.' && isDigit_(peek_(1)))
      {
        advance_();
        while (!isAtEnd_() && isDigit_(current_())) advance_();
        return Token{TokenType::NUMBER, input_.substr(start_pos, pos_ - start_pos), start_pos};
      }
      while (!isAtEnd_() && isDigit_(current_())) advance_();
      if (!isAtEnd_() && current_() == '.' && isDigit_(peek_(1)))
      {
        advance_();
        while (!isAtEnd_() && isDigit_(current_())) advance_();
      }
      return Token{TokenType::NUMBER, input_.substr(start_pos, pos_ - start_pos), start_pos};
    }

    Token scanIdentifier_()
    {
      size_t start_pos = pos_;
      while (!isAtEnd_() && isLetter_(current_())) advance_();
      return Token{TokenType::IDENTIFIER, input_.substr(start_pos, pos_ - start_pos), start_pos};
    }

    bool isAtEnd_() const { return pos_ >= input_.size(); }
    char current_() const { return isAtEnd_() ? '\0' : input_[pos_]; }
    char peek_(size_t offset) const { size_t i = pos_ + offset; return i >= input_.size() ? '\0' : input_[i]; }
    char advance_() { return isAtEnd_() ? '\0' : input_[pos_++]; }
    static bool isLetter_(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }
    static bool isDigit_(char c) { return c >= '0' && c <= '9'; }
  };

} // namespace detail

//============================================================================
// Internal ProFormaWriter class
//============================================================================
namespace detail
{
  // Type aliases for convenience within detail namespace
  using Peptidoform = ProForma::Peptidoform;
  using PeptidoformIon = ProForma::PeptidoformIon;
  using WriteMode = ProForma::WriteMode;
  using GlobalModEntry = ProForma::GlobalModEntry;
  using IsotopeReplacement = ProForma::IsotopeReplacement;
  using GlobalModification = ProForma::GlobalModification;
  using Modification = ProForma::Modification;
  using ModificationTag = ProForma::ModificationTag;
  using CvAccession = ProForma::CvAccession;
  using NamedMod = ProForma::NamedMod;
  using MassDelta = ProForma::MassDelta;
  using FormulaTag = ProForma::FormulaTag;
  using GlycanComposition = ProForma::GlycanComposition;
  using InfoTag = ProForma::InfoTag;
  using PositionConstraint = ProForma::PositionConstraint;
  using Label = ProForma::Label;
  using UnlocalisedMod = ProForma::UnlocalisedMod;
  using LabileModification = ProForma::LabileModification;
  using SequenceSection = ProForma::SequenceSection;
  using SequenceElement = ProForma::SequenceElement;
  using AmbiguousRegion = ProForma::AmbiguousRegion;
  using ModifiedRange = ProForma::ModifiedRange;
  using ChargeState = ProForma::ChargeState;
  using CvDatabase = ProForma::CvDatabase;
  using AdductIon = ProForma::AdductIon;
  using ConversionIssue = ProForma::ConversionIssue;
  using ConversionIssueType = ProForma::ConversionIssueType;

  class ProFormaWriter
  {
  public:
    static String toString(const Peptidoform& peptidoform, WriteMode mode);
    static String toString(const PeptidoformIon& ion, WriteMode mode);

  private:
    static void writeGlobalMods_(std::ostream& os, const std::vector<GlobalModEntry>& mods, WriteMode mode);
    static void writeIsotopeReplacement_(std::ostream& os, const IsotopeReplacement& isotope);
    static void writeGlobalModification_(std::ostream& os, const GlobalModification& mod, WriteMode mode);
    static void writeModification_(std::ostream& os, const Modification& mod, WriteMode mode);
    static void writeModificationTag_(std::ostream& os, const ModificationTag& tag, WriteMode mode);
    static void writeCvAccession_(std::ostream& os, const CvAccession& cv);
    static void writeNamedMod_(std::ostream& os, const NamedMod& named);
    static void writeMassDelta_(std::ostream& os, const MassDelta& delta, WriteMode mode);
    static void writeFormulaTag_(std::ostream& os, const FormulaTag& formula);
    static void writeGlycanComposition_(std::ostream& os, const GlycanComposition& glycan);
    static void writeInfoTag_(std::ostream& os, const InfoTag& info);
    static void writePositionConstraint_(std::ostream& os, const PositionConstraint& pc);
    static void writeLabel_(std::ostream& os, const Label& label, WriteMode mode);
    static void writeUnlocalisedMods_(std::ostream& os, const std::vector<UnlocalisedMod>& mods, WriteMode mode);
    static void writeLabileModifications_(std::ostream& os, const std::vector<LabileModification>& mods, WriteMode mode);
    static void writeNTermMods_(std::ostream& os, const std::vector<Modification>& mods, WriteMode mode);
    static void writeCTermMods_(std::ostream& os, const std::vector<Modification>& mods, WriteMode mode);
    static void writeSequence_(std::ostream& os, const std::vector<SequenceSection>& seq, WriteMode mode);
    static void writeSequenceElement_(std::ostream& os, const SequenceElement& elem, WriteMode mode);
    static void writeAmbiguousRegion_(std::ostream& os, const AmbiguousRegion& region, WriteMode mode);
    static void writeModifiedRange_(std::ostream& os, const ModifiedRange& range, WriteMode mode);
    static void writeChargeState_(std::ostream& os, const ChargeState& charge);
    static const char* cvDatabaseToString_(CvDatabase db);
    static const char* cvDatabaseToHint_(CvDatabase db);
    static const char* massSourceToString_(MassDelta::Source source);
  };

  String ProFormaWriter::toString(const Peptidoform& peptidoform, WriteMode mode)
  {
    std::ostringstream os;
    if (peptidoform.name.has_value()) os << "(>" << peptidoform.name.value() << ")";
    writeGlobalMods_(os, peptidoform.global_mods, mode);
    writeUnlocalisedMods_(os, peptidoform.unlocalised_mods, mode);
    writeLabileModifications_(os, peptidoform.labile_mods, mode);
    writeNTermMods_(os, peptidoform.n_term_mods, mode);
    writeSequence_(os, peptidoform.sequence, mode);
    writeCTermMods_(os, peptidoform.c_term_mods, mode);
    return String(os.str());
  }

  String ProFormaWriter::toString(const PeptidoformIon& ion, WriteMode mode)
  {
    std::ostringstream os;
    const char* separator = ion.is_chimeric ? "+" : "//";
    bool first = true;
    for (const auto& chain : ion.chains)
    {
      if (!first) os << separator;
      first = false;
      os << toString(chain, mode);
      if (ion.is_chimeric && chain.charge.has_value()) writeChargeState_(os, chain.charge.value());
    }
    if (ion.charge.has_value()) writeChargeState_(os, ion.charge.value());
    return String(os.str());
  }

  void ProFormaWriter::writeGlobalMods_(std::ostream& os, const std::vector<GlobalModEntry>& mods, WriteMode mode)
  {
    for (const auto& entry : mods)
    {
      std::visit([&os, mode](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, IsotopeReplacement>) writeIsotopeReplacement_(os, arg);
        else if constexpr (std::is_same_v<T, GlobalModification>) writeGlobalModification_(os, arg, mode);
      }, entry);
    }
  }

  void ProFormaWriter::writeIsotopeReplacement_(std::ostream& os, const IsotopeReplacement& isotope)
  {
    os << '<' << isotope.isotope << '>';
  }

  void ProFormaWriter::writeGlobalModification_(std::ostream& os, const GlobalModification& mod, WriteMode mode)
  {
    os << '<';
    writeModification_(os, mod.modification, mode);
    if (!mod.locations.empty())
    {
      os << '@';
      bool first = true;
      for (const auto& loc : mod.locations) { if (!first) os << ','; first = false; os << loc; }
    }
    os << '>';
  }

  void ProFormaWriter::writeModification_(std::ostream& os, const Modification& mod, WriteMode mode)
  {
    os << '[';
    bool first = true;
    for (const auto& alt : mod.alternatives)
    {
      if (!first) os << '|';
      first = false;
      bool is_label_only = false;
      if (const auto* info = std::get_if<InfoTag>(&alt.first))
      {
        if (info->text.empty() && alt.second.has_value()) is_label_only = true;
      }
      if (!is_label_only) writeModificationTag_(os, alt.first, mode);
      if (alt.second.has_value()) writeLabel_(os, alt.second.value(), mode);
    }
    os << ']';
  }

  void ProFormaWriter::writeModificationTag_(std::ostream& os, const ModificationTag& tag, WriteMode mode)
  {
    std::visit([&os, mode](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, CvAccession>) writeCvAccession_(os, arg);
      else if constexpr (std::is_same_v<T, NamedMod>) writeNamedMod_(os, arg);
      else if constexpr (std::is_same_v<T, MassDelta>) writeMassDelta_(os, arg, mode);
      else if constexpr (std::is_same_v<T, FormulaTag>) writeFormulaTag_(os, arg);
      else if constexpr (std::is_same_v<T, GlycanComposition>) writeGlycanComposition_(os, arg);
      else if constexpr (std::is_same_v<T, InfoTag>) writeInfoTag_(os, arg);
      else if constexpr (std::is_same_v<T, PositionConstraint>) writePositionConstraint_(os, arg);
    }, tag);
  }

  void ProFormaWriter::writeCvAccession_(std::ostream& os, const CvAccession& cv)
  {
    os << cvDatabaseToString_(cv.database) << ':' << cv.accession;
  }

  void ProFormaWriter::writeNamedMod_(std::ostream& os, const NamedMod& named)
  {
    if (named.cv_hint.has_value()) os << cvDatabaseToHint_(named.cv_hint.value()) << ':';
    os << named.name;
  }

  void ProFormaWriter::writeMassDelta_(std::ostream& os, const MassDelta& delta, WriteMode mode)
  {
    if (delta.source != MassDelta::Source::NONE) os << massSourceToString_(delta.source) << ':';
    if (mode == WriteMode::LOSSLESS && !delta.original_text.empty()) os << delta.original_text;
    else { if (delta.mass >= 0) os << '+'; os << std::fixed << std::setprecision(4) << delta.mass; }
  }

  void ProFormaWriter::writeFormulaTag_(std::ostream& os, const FormulaTag& formula)
  {
    os << "Formula:" << formula.formula_string;
    if (formula.charge.has_value()) { os << ":z"; int c = formula.charge.value(); if (c >= 0) os << '+'; os << c; }
  }

  void ProFormaWriter::writeGlycanComposition_(std::ostream& os, const GlycanComposition& glycan)
  {
    os << "Glycan:";
    for (const auto& component : glycan.components)
    {
      std::visit([&os](auto&& mono) {
        using T = std::decay_t<decltype(mono)>;
        if constexpr (std::is_same_v<T, String>) os << mono;
        else if constexpr (std::is_same_v<T, FormulaTag>) {
          os << "Formula:" << mono.formula_string;
          if (mono.charge.has_value()) { os << ":z"; int c = mono.charge.value(); if (c >= 0) os << '+'; os << c; }
        }
      }, component.first);
      os << component.second;
    }
  }

  void ProFormaWriter::writeInfoTag_(std::ostream& os, const InfoTag& info) { os << "INFO:" << info.text; }

  void ProFormaWriter::writePositionConstraint_(std::ostream& os, const PositionConstraint& pc)
  {
    os << "Position:";
    bool first = true;
    if (pc.n_term) { os << "N-term"; first = false; }
    if (pc.c_term) { if (!first) os << ','; os << "C-term"; first = false; }
    if (!pc.residues.empty()) { if (!first) os << ','; os << std::string(pc.residues.begin(), pc.residues.end()); }
  }

  void ProFormaWriter::writeLabel_(std::ostream& os, const Label& label, WriteMode mode)
  {
    os << '#' << label.identifier;
    if (label.score.has_value())
    {
      if (mode == WriteMode::CANONICAL) os << '(' << std::fixed << std::setprecision(2) << label.score.value() << ')';
      else os << '(' << label.score.value() << ')';
    }
  }

  void ProFormaWriter::writeUnlocalisedMods_(std::ostream& os, const std::vector<UnlocalisedMod>& mods, WriteMode mode)
  {
    for (const auto& unloc : mods)
    {
      for (const auto& mod : unloc.modifications) writeModification_(os, mod, mode);
      if (unloc.occurrence.has_value()) os << '^' << unloc.occurrence.value();
      os << '?';
    }
  }

  void ProFormaWriter::writeLabileModifications_(std::ostream& os, const std::vector<LabileModification>& mods, WriteMode mode)
  {
    for (const auto& labile : mods)
    {
      os << '{';
      bool first = true;
      for (const auto& alt : labile.modification.alternatives)
      {
        if (!first) os << '|';
        first = false;
        writeModificationTag_(os, alt.first, mode);
        if (alt.second.has_value()) writeLabel_(os, alt.second.value(), mode);
      }
      os << '}';
    }
  }

  void ProFormaWriter::writeNTermMods_(std::ostream& os, const std::vector<Modification>& mods, WriteMode mode)
  {
    if (!mods.empty()) { for (const auto& mod : mods) writeModification_(os, mod, mode); os << '-'; }
  }

  void ProFormaWriter::writeCTermMods_(std::ostream& os, const std::vector<Modification>& mods, WriteMode mode)
  {
    if (!mods.empty()) { os << '-'; for (const auto& mod : mods) writeModification_(os, mod, mode); }
  }

  void ProFormaWriter::writeSequence_(std::ostream& os, const std::vector<SequenceSection>& seq, WriteMode mode)
  {
    for (const auto& section : seq)
    {
      std::visit([&os, mode](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, SequenceElement>) writeSequenceElement_(os, arg, mode);
        else if constexpr (std::is_same_v<T, AmbiguousRegion>) writeAmbiguousRegion_(os, arg, mode);
        else if constexpr (std::is_same_v<T, ModifiedRange>) writeModifiedRange_(os, arg, mode);
      }, section);
    }
  }

  void ProFormaWriter::writeSequenceElement_(std::ostream& os, const SequenceElement& elem, WriteMode mode)
  {
    os << elem.amino_acid;
    for (const auto& mod : elem.modifications) writeModification_(os, mod, mode);
  }

  void ProFormaWriter::writeAmbiguousRegion_(std::ostream& os, const AmbiguousRegion& region, WriteMode mode)
  {
    os << "(?";
    for (const auto& elem : region.elements) writeSequenceElement_(os, elem, mode);
    os << ')';
  }

  void ProFormaWriter::writeModifiedRange_(std::ostream& os, const ModifiedRange& range, WriteMode mode)
  {
    os << '(';
    for (const auto& elem : range.elements) writeSequenceElement_(os, elem, mode);
    os << ')';
    for (const auto& mod : range.modifications) writeModification_(os, mod, mode);
  }

  void ProFormaWriter::writeChargeState_(std::ostream& os, const ChargeState& charge)
  {
    os << '/';
    std::visit([&os](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, int>) os << arg;
      else if constexpr (std::is_same_v<T, std::vector<AdductIon>>)
      {
        os << '[';
        int total_charge = 0;
        bool first = true;
        for (const auto& adduct : arg)
        {
          if (!first) os << ',';
          first = false;
          os << adduct.formula << ":z";
          if (adduct.charge >= 0) os << '+';
          os << adduct.charge;
          if (adduct.occurrence.has_value()) os << '^' << adduct.occurrence.value();
          total_charge += adduct.charge * adduct.occurrence.value_or(1);
        }
        os << ']' << std::abs(total_charge) << (total_charge >= 0 ? '+' : '-');
      }
    }, charge);
  }

  const char* ProFormaWriter::cvDatabaseToString_(CvDatabase db)
  {
    switch (db) {
      case CvDatabase::UNIMOD: return "UNIMOD";
      case CvDatabase::MOD:    return "MOD";
      case CvDatabase::RESID:  return "RESID";
      case CvDatabase::XLMOD:  return "XLMOD";
      case CvDatabase::GNO:    return "GNO";
      default: return "UNIMOD";
    }
  }

  const char* ProFormaWriter::cvDatabaseToHint_(CvDatabase db)
  {
    switch (db) {
      case CvDatabase::UNIMOD: return "U";
      case CvDatabase::MOD:    return "M";
      case CvDatabase::RESID:  return "R";
      case CvDatabase::XLMOD:  return "X";
      case CvDatabase::GNO:    return "G";
      default: return "U";
    }
  }

  const char* ProFormaWriter::massSourceToString_(MassDelta::Source source)
  {
    switch (source) {
      case MassDelta::Source::OBS: return "Obs";
      case MassDelta::Source::U:   return "U";
      case MassDelta::Source::M:   return "M";
      case MassDelta::Source::R:   return "R";
      case MassDelta::Source::X:   return "X";
      case MassDelta::Source::G:   return "G";
      default: return "";
    }
  }

} // namespace detail

//============================================================================
// ProForma::ParseError implementation
//============================================================================

const char* ProForma::errorCodeToString(ErrorCode code)
{
  switch (code)
  {
    case ErrorCode::UNEXPECTED_CHARACTER: return "Unexpected character";
    case ErrorCode::UNCLOSED_BRACKET: return "Unclosed bracket";
    case ErrorCode::UNMATCHED_BRACKET: return "Unmatched closing bracket";
    case ErrorCode::INVALID_CV_PREFIX: return "Invalid controlled vocabulary prefix";
    case ErrorCode::INVALID_CV_ACCESSION: return "Invalid CV accession number";
    case ErrorCode::INVALID_AMINO_ACID: return "Invalid amino acid";
    case ErrorCode::INVALID_MASS_VALUE: return "Invalid mass value";
    case ErrorCode::INVALID_FORMULA: return "Invalid chemical formula";
    case ErrorCode::UNKNOWN_MONOSACCHARIDE: return "Unknown monosaccharide";
    case ErrorCode::DANGLING_CROSSLINK_LABEL: return "Dangling crosslink label";
    case ErrorCode::EMPTY_SEQUENCE: return "Empty sequence";
    case ErrorCode::INVALID_CHARGE: return "Invalid charge state";
    case ErrorCode::INVALID_OCCURRENCE_SPECIFIER: return "Invalid occurrence specifier";
    case ErrorCode::UNEXPECTED_END_OF_INPUT: return "Unexpected end of input";
    case ErrorCode::INTERNAL_ERROR: return "Internal parser error";
    default: return "Unknown error";
  }
}

ProForma::ParseError::ParseError(
  const char* file, int line, const char* function,
  ErrorCode error_code, size_t error_position,
  const String& input, const String& message) noexcept :
  Exception::ParseError(file, line, function, input, message),
  code_(error_code),
  position_(std::min(error_position, input.size()))
{
  extractContext_(input, position_);
  Exception::GlobalExceptionHandler::getInstance().setMessage(what());
}

void ProForma::ParseError::extractContext_(const String& input, size_t pos)
{
  const size_t context_length = 20;
  if (pos > 0)
  {
    size_t start = (pos > context_length) ? pos - context_length : 0;
    context_before_ = input.substr(start, pos - start);
  }
  else context_before_ = "";
  if (pos < input.size())
  {
    size_t length = std::min(context_length, input.size() - pos);
    context_after_ = input.substr(pos, length);
  }
  else context_after_ = "";
}

String ProForma::ParseError::getFormattedMessage() const
{
  std::ostringstream oss;
  oss << "ProForma parse error at position " << position_ << ": " << ProForma::errorCodeToString(code_);
  oss << "\nContext: ";
  if (position_ > context_before_.size()) oss << "...";
  oss << context_before_;
  if (!context_after_.empty()) { oss << ">>>" << context_after_.substr(0, 1) << "<<<"; if (context_after_.size() > 1) oss << context_after_.substr(1); }
  else oss << ">>><END OF INPUT><<<";
  if (context_after_.size() >= 20) oss << "...";
  if (!expected_.empty() || !found_.empty())
  {
    if (!expected_.empty()) oss << "\nExpected: " << expected_;
    if (!found_.empty()) oss << "\nFound: " << found_;
  }
  return String(oss.str());
}

void ProForma::ParseError::setExpectedFound(const String& expected, const String& found)
{
  expected_ = expected;
  found_ = found;
}

//============================================================================
// JSON implementation (delegates to ProFormaDataJson.h inline functions)
//============================================================================

String ProForma::toJSON(const Peptidoform& pf)
{
  nlohmann::json j = pf;
  return String(j.dump());
}

ProForma::Peptidoform ProForma::peptidoformFromJSON(const String& json_str)
{
  try
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<Peptidoform>();
  }
  catch (const nlohmann::json::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, json_str, String("JSON parsing failed: ") + e.what());
  }
  catch (const std::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, json_str, String("JSON deserialization failed: ") + e.what());
  }
}

String ProForma::toJSON(const PeptidoformIon& pfi)
{
  nlohmann::json j = pfi;
  return String(j.dump());
}

ProForma::PeptidoformIon ProForma::peptidoformIonFromJSON(const String& json_str)
{
  try
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<PeptidoformIon>();
  }
  catch (const nlohmann::json::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, json_str, String("JSON parsing failed: ") + e.what());
  }
  catch (const std::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, json_str, String("JSON deserialization failed: ") + e.what());
  }
}

//============================================================================
// Internal ProFormaParserImpl class
//============================================================================
namespace detail
{
  using TokenType = detail::TokenType;
  using Token = detail::Token;
  using Tokenizer = detail::ProFormaTokenizer;

  // Type aliases for convenience within detail namespace
  using Peptidoform = ProForma::Peptidoform;
  using PeptidoformIon = ProForma::PeptidoformIon;
  using WriteMode = ProForma::WriteMode;
  using ErrorCode = ProForma::ErrorCode;
  using GlobalModEntry = ProForma::GlobalModEntry;
  using IsotopeReplacement = ProForma::IsotopeReplacement;
  using GlobalModification = ProForma::GlobalModification;
  using Modification = ProForma::Modification;
  using ModificationTag = ProForma::ModificationTag;
  using CvAccession = ProForma::CvAccession;
  using NamedMod = ProForma::NamedMod;
  using MassDelta = ProForma::MassDelta;
  using FormulaTag = ProForma::FormulaTag;
  using GlycanComposition = ProForma::GlycanComposition;
  using InfoTag = ProForma::InfoTag;
  using PositionConstraint = ProForma::PositionConstraint;
  using Label = ProForma::Label;
  using UnlocalisedMod = ProForma::UnlocalisedMod;
  using LabileModification = ProForma::LabileModification;
  using SequenceSection = ProForma::SequenceSection;
  using SequenceElement = ProForma::SequenceElement;
  using AmbiguousRegion = ProForma::AmbiguousRegion;
  using ModifiedRange = ProForma::ModifiedRange;
  using ChargeState = ProForma::ChargeState;
  using CvDatabase = ProForma::CvDatabase;
  using AdductIon = ProForma::AdductIon;
  using ConversionIssue = ProForma::ConversionIssue;
  using ConversionIssueType = ProForma::ConversionIssueType;
  using ConversionPolicy = ProForma::ConversionPolicy;

  class ProFormaParserImpl
  {
  public:
    explicit ProFormaParserImpl(std::string_view input) : tokenizer_(input), input_(input), current_token_{TokenType::END, {}, 0}, has_current_(false) {}

    PeptidoformIon parsePeptidoformIon();
    Peptidoform parsePeptidoform();

  private:
    Tokenizer tokenizer_;
    std::string_view input_;
    Token current_token_;
    bool has_current_;

    Peptidoform parsePeptidoformWithCharge_(bool is_chimeric_context);
    Peptidoform parsePeptidoform_();
    Tokenizer createLookahead_() const;
    bool hasNTerminalModPattern_();
    std::vector<GlobalModEntry> parseGlobalMods_();
    GlobalModEntry parseGlobalModEntry_();
    IsotopeReplacement parseIsotopeReplacement_();
    GlobalModification parseGlobalModification_();
    std::vector<UnlocalisedMod> parseUnlocalisedMods_();
    std::vector<LabileModification> parseLabileModifications_();
    std::vector<SequenceSection> parseSequence_();
    SequenceElement parseSequenceElement_();
    AmbiguousRegion parseAmbiguousRegion_();
    ModifiedRange parseModifiedRange_();
    std::vector<Modification> parseTerminalMods_();
    std::vector<Modification> parseModificationList_();
    Modification parseModification_();
    std::pair<ModificationTag, std::optional<Label>> parseModificationTagWithLabel_();
    ModificationTag parseModificationTag_();
    NamedMod parseNamedMod_();
    NamedMod parseNamedMod_(char cv_hint);
    CvAccession parseCvAccession_();
    MassDelta parseMassDelta_();
    FormulaTag parseFormulaTag_();
    GlycanComposition parseGlycanComposition_();
    InfoTag parseInfoTag_();
    PositionConstraint parsePositionConstraint_();
    Label parseLabel_();
    std::optional<ChargeState> parseChargeState_();
    std::vector<AdductIon> parseAdductIons_();
    AdductIon parseAdductIon_();

    Token current_();
    Token advance_();
    bool check_(TokenType type);
    bool match_(TokenType type);
    Token expect_(TokenType type, const char* expected_desc);
    bool isAtEnd_();
    void error_(ErrorCode code, const char* message);
    void errorAt_(ErrorCode code, size_t pos, const char* message);
    static bool isAminoAcid_(char c);
  };

  PeptidoformIon ProFormaParserImpl::parsePeptidoformIon()
  {
    PeptidoformIon ion;
    ion.chains.push_back(parsePeptidoformWithCharge_(false));
    while (!isAtEnd_())
    {
      if (match_(TokenType::SLASH))
      {
        if (!match_(TokenType::SLASH)) { ion.charge = parseChargeState_(); break; }
        ion.chains.push_back(parsePeptidoformWithCharge_(false));
      }
      else if (match_(TokenType::PLUS)) { ion.is_chimeric = true; ion.chains.push_back(parsePeptidoformWithCharge_(true)); }
      else break;
    }
    if (!isAtEnd_()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Unexpected characters after peptidoform ion");
    return ion;
  }

  Peptidoform ProFormaParserImpl::parsePeptidoform()
  {
    Peptidoform pf = parsePeptidoform_();
    if (!isAtEnd_()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Unexpected characters after peptidoform");
    return pf;
  }

  Peptidoform ProFormaParserImpl::parsePeptidoformWithCharge_(bool is_chimeric_context)
  {
    Peptidoform pf = parsePeptidoform_();
    if (check_(TokenType::SLASH))
    {
      Tokenizer lookahead = createLookahead_();
      lookahead.next();
      Token next = lookahead.peek();
      if (next.type == TokenType::SLASH) return pf;
      if (next.type == TokenType::PLUS || next.type == TokenType::MINUS || next.type == TokenType::NUMBER || next.type == TokenType::LBRACKET)
      {
        if (next.type == TokenType::LBRACKET) { int depth = 1; lookahead.next(); while (lookahead.hasMore() && depth > 0) { Token tok = lookahead.next(); if (tok.type == TokenType::LBRACKET) depth++; else if (tok.type == TokenType::RBRACKET) depth--; } }
        else { if (next.type == TokenType::PLUS || next.type == TokenType::MINUS) lookahead.next(); if (lookahead.peek().type == TokenType::NUMBER) lookahead.next(); }
        Token after_charge = lookahead.peek();
        bool followed_by_chimeric = (after_charge.type == TokenType::PLUS);
        bool at_end = (after_charge.type == TokenType::END);
        if (followed_by_chimeric || (at_end && is_chimeric_context)) { advance_(); pf.charge = parseChargeState_(); }
      }
    }
    return pf;
  }

  Peptidoform ProFormaParserImpl::parsePeptidoform_()
  {
    Peptidoform pf;
    std::string combined_name;
    while (check_(TokenType::LPAREN))
    {
      Tokenizer lookahead = createLookahead_();
      lookahead.next();
      if (lookahead.peek().type != TokenType::RANGLE) break;
      advance_();
      while (check_(TokenType::RANGLE)) advance_();
      std::string name;
      int paren_depth = 0;
      while (!isAtEnd_())
      {
        Token tok = current_();
        if (tok.type == TokenType::LPAREN) { paren_depth++; name += std::string(tok.text); advance_(); }
        else if (tok.type == TokenType::RPAREN) { if (paren_depth > 0) { paren_depth--; name += std::string(tok.text); advance_(); } else break; }
        else { name += std::string(tok.text); advance_(); }
      }
      expect_(TokenType::RPAREN, "')' to close gene/protein prefix");
      if (!combined_name.empty()) combined_name += " / ";
      combined_name += name;
    }
    if (!combined_name.empty()) pf.name = String(combined_name);
    while (check_(TokenType::LANGLE)) { auto mods = parseGlobalMods_(); for (auto& m : mods) pf.global_mods.push_back(std::move(m)); }
    pf.unlocalised_mods = parseUnlocalisedMods_();
    pf.labile_mods = parseLabileModifications_();
    if (check_(TokenType::LBRACKET))
    {
      if (hasNTerminalModPattern_()) { pf.n_term_mods = parseTerminalMods_(); expect_(TokenType::MINUS, "'-' after N-terminal modification"); }
    }
    pf.sequence = parseSequence_();
    if (pf.sequence.empty()) error_(ErrorCode::EMPTY_SEQUENCE, "Empty sequence");
    if (match_(TokenType::MINUS))
    {
      pf.c_term_mods = parseTerminalMods_();
      if (pf.c_term_mods.empty()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected modification after '-' for C-terminal modification");
    }
    return pf;
  }

  Tokenizer ProFormaParserImpl::createLookahead_() const
  {
    size_t start_pos = has_current_ ? current_token_.position : tokenizer_.position();
    return Tokenizer(input_, start_pos);
  }

  bool ProFormaParserImpl::hasNTerminalModPattern_()
  {
    Tokenizer lookahead = createLookahead_();
    if (lookahead.peek().type != TokenType::LBRACKET) return false;
    while (lookahead.peek().type == TokenType::LBRACKET)
    {
      int depth = 1; lookahead.next();
      while (lookahead.hasMore() && depth > 0) { Token tok = lookahead.next(); if (tok.type == TokenType::LBRACKET) depth++; else if (tok.type == TokenType::RBRACKET) depth--; }
      if (depth != 0) return false;
    }
    return lookahead.peek().type == TokenType::MINUS;
  }

  std::vector<GlobalModEntry> ProFormaParserImpl::parseGlobalMods_()
  {
    std::vector<GlobalModEntry> entries;
    if (!match_(TokenType::LANGLE)) return entries;
    entries.push_back(parseGlobalModEntry_());
    while (match_(TokenType::COMMA)) entries.push_back(parseGlobalModEntry_());
    expect_(TokenType::RANGLE, "'>' to close global modifications");
    return entries;
  }

  GlobalModEntry ProFormaParserImpl::parseGlobalModEntry_()
  {
    if (check_(TokenType::LBRACKET)) return parseGlobalModification_();
    else return parseIsotopeReplacement_();
  }

  IsotopeReplacement ProFormaParserImpl::parseIsotopeReplacement_()
  {
    IsotopeReplacement ir;
    std::string isotope_str;
    if (check_(TokenType::NUMBER)) { Token num = advance_(); isotope_str = std::string(num.text); }
    if (check_(TokenType::IDENTIFIER)) { Token id = advance_(); isotope_str += std::string(id.text); }
    if (isotope_str.empty()) error_(ErrorCode::INVALID_FORMULA, "Expected isotope specification");
    ir.isotope = isotope_str;
    return ir;
  }

  GlobalModification ProFormaParserImpl::parseGlobalModification_()
  {
    GlobalModification gm;
    expect_(TokenType::LBRACKET, "'['");
    gm.modification = parseModification_();
    expect_(TokenType::RBRACKET, "']'");
    for (const auto& alt : gm.modification.alternatives) { if (alt.second.has_value()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Labels are not allowed on global modifications"); }
    expect_(TokenType::AT, "'@' for global modification locations");
    std::string location;
    if (check_(TokenType::IDENTIFIER))
    {
      Token id = advance_();
      location = std::string(id.text);
      if (match_(TokenType::MINUS)) { Token term = expect_(TokenType::IDENTIFIER, "term identifier"); location += "-" + std::string(term.text); }
      if (match_(TokenType::COLON)) { Token suffix = expect_(TokenType::IDENTIFIER, "location suffix"); location += ":" + std::string(suffix.text); }
      gm.locations.push_back(location);
    }
    while (check_(TokenType::COMMA))
    {
      Tokenizer lookahead = tokenizer_;
      Token after_comma = lookahead.peek();
      if (after_comma.type == TokenType::LBRACKET || after_comma.type == TokenType::NUMBER) break;
      advance_();
      location.clear();
      if (check_(TokenType::IDENTIFIER))
      {
        Token id = advance_();
        location = std::string(id.text);
        if (match_(TokenType::MINUS)) { Token term = expect_(TokenType::IDENTIFIER, "term identifier"); location += "-" + std::string(term.text); }
        if (match_(TokenType::COLON)) { Token suffix = expect_(TokenType::IDENTIFIER, "location suffix"); location += ":" + std::string(suffix.text); }
        gm.locations.push_back(location);
      }
    }
    return gm;
  }

  std::vector<UnlocalisedMod> ProFormaParserImpl::parseUnlocalisedMods_()
  {
    std::vector<UnlocalisedMod> mods;
    while (check_(TokenType::LBRACKET))
    {
      Tokenizer lookahead = createLookahead_();
      int group_count = 0;
      bool found_question = false;
      while (true)
      {
        Token tok = lookahead.peek();
        if (tok.type != TokenType::LBRACKET) break;
        int depth = 0;
        tok = lookahead.next();
        depth = 1;
        group_count++;
        while (lookahead.hasMore() && depth > 0) { tok = lookahead.next(); if (tok.type == TokenType::LBRACKET) depth++; else if (tok.type == TokenType::RBRACKET) depth--; }
        tok = lookahead.peek();
        if (tok.type == TokenType::CARET) { lookahead.next(); tok = lookahead.peek(); if (tok.type == TokenType::NUMBER) { lookahead.next(); tok = lookahead.peek(); } }
        if (tok.type == TokenType::QUESTION) { found_question = true; break; }
        if (tok.type != TokenType::LBRACKET) break;
      }
      if (!found_question) break;
      while (check_(TokenType::LBRACKET))
      {
        UnlocalisedMod um;
        expect_(TokenType::LBRACKET, "'['");
        um.modifications.push_back(parseModification_());
        expect_(TokenType::RBRACKET, "']'");
        if (match_(TokenType::CARET))
        {
          Token num = expect_(TokenType::NUMBER, "occurrence count");
          try { um.occurrence = std::stoi(std::string(num.text)); }
          catch (const std::exception&) { errorAt_(ErrorCode::INVALID_MASS_VALUE, num.position, "Invalid occurrence count"); }
        }
        mods.push_back(std::move(um));
        if (check_(TokenType::QUESTION)) break;
      }
      expect_(TokenType::QUESTION, "'?' for unlocalised modification");
    }
    return mods;
  }

  std::vector<LabileModification> ProFormaParserImpl::parseLabileModifications_()
  {
    std::vector<LabileModification> mods;
    while (match_(TokenType::LBRACE))
    {
      LabileModification lm;
      lm.modification = parseModification_();
      for (const auto& alt : lm.modification.alternatives) { if (alt.second.has_value()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Labels are not allowed on labile modifications"); }
      expect_(TokenType::RBRACE, "'}'");
      mods.push_back(std::move(lm));
    }
    return mods;
  }

  std::vector<SequenceSection> ProFormaParserImpl::parseSequence_()
  {
    std::vector<SequenceSection> sections;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::MINUS || tok.type == TokenType::SLASH || tok.type == TokenType::END) break;
      if (tok.type == TokenType::LPAREN)
      {
        Tokenizer lookahead = createLookahead_();
        lookahead.next();
        if (lookahead.peek().type == TokenType::QUESTION) sections.push_back(parseAmbiguousRegion_());
        else sections.push_back(parseModifiedRange_());
        continue;
      }
      if (tok.type == TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();
        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c)) errorAt_(ErrorCode::INVALID_AMINO_ACID, tok.position, "Invalid amino acid character");
          SequenceElement elem;
          elem.amino_acid = c;
          if (i == text.size() - 1 && check_(TokenType::LBRACKET)) elem.modifications = parseModificationList_();
          sections.push_back(elem);
        }
        continue;
      }
      break;
    }
    return sections;
  }

  SequenceElement ProFormaParserImpl::parseSequenceElement_()
  {
    SequenceElement elem;
    Token tok = expect_(TokenType::IDENTIFIER, "amino acid");
    if (tok.text.size() != 1 || !isAminoAcid_(tok.text[0])) errorAt_(ErrorCode::INVALID_AMINO_ACID, tok.position, "Expected single amino acid letter");
    elem.amino_acid = tok.text[0];
    if (check_(TokenType::LBRACKET)) elem.modifications = parseModificationList_();
    return elem;
  }

  AmbiguousRegion ProFormaParserImpl::parseAmbiguousRegion_()
  {
    AmbiguousRegion region;
    expect_(TokenType::LPAREN, "'('");
    expect_(TokenType::QUESTION, "'?' for ambiguous region");
    while (!check_(TokenType::RPAREN) && !isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();
        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c)) errorAt_(ErrorCode::INVALID_AMINO_ACID, tok.position, "Invalid amino acid in ambiguous region");
          SequenceElement elem;
          elem.amino_acid = c;
          if (i == text.size() - 1 && check_(TokenType::LBRACKET)) elem.modifications = parseModificationList_();
          region.elements.push_back(elem);
        }
      }
      else break;
    }
    expect_(TokenType::RPAREN, "')'");
    return region;
  }

  ModifiedRange ProFormaParserImpl::parseModifiedRange_()
  {
    ModifiedRange range;
    expect_(TokenType::LPAREN, "'('");
    while (!check_(TokenType::RPAREN) && !isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::IDENTIFIER)
      {
        std::string_view text = tok.text;
        advance_();
        for (size_t i = 0; i < text.size(); ++i)
        {
          char c = text[i];
          if (!isAminoAcid_(c)) errorAt_(ErrorCode::INVALID_AMINO_ACID, tok.position + i, "Invalid amino acid in range");
          SequenceElement elem;
          elem.amino_acid = c;
          if (i == text.size() - 1 && check_(TokenType::LBRACKET)) elem.modifications = parseModificationList_();
          range.elements.push_back(elem);
        }
      }
      else if (tok.type == TokenType::LBRACKET)
      {
        if (!range.elements.empty()) { auto mods = parseModificationList_(); for (auto& m : mods) range.elements.back().modifications.push_back(std::move(m)); }
        else break;
      }
      else break;
    }
    expect_(TokenType::RPAREN, "')'");
    if (range.elements.empty()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Empty range is not allowed");
    if (check_(TokenType::LBRACKET)) range.modifications = parseModificationList_();
    return range;
  }

  std::vector<Modification> ProFormaParserImpl::parseTerminalMods_()
  {
    std::vector<Modification> mods;
    while (check_(TokenType::LBRACKET)) { expect_(TokenType::LBRACKET, "'['"); mods.push_back(parseModification_()); expect_(TokenType::RBRACKET, "']'"); }
    return mods;
  }

  std::vector<Modification> ProFormaParserImpl::parseModificationList_()
  {
    std::vector<Modification> mods;
    while (match_(TokenType::LBRACKET)) { mods.push_back(parseModification_()); expect_(TokenType::RBRACKET, "']'"); }
    return mods;
  }

  Modification ProFormaParserImpl::parseModification_()
  {
    Modification mod;
    auto [tag, label] = parseModificationTagWithLabel_();
    mod.alternatives.emplace_back(std::move(tag), std::move(label));
    while (match_(TokenType::PIPE)) { auto [alt_tag, alt_label] = parseModificationTagWithLabel_(); mod.alternatives.emplace_back(std::move(alt_tag), std::move(alt_label)); }
    return mod;
  }

  std::pair<ModificationTag, std::optional<Label>> ProFormaParserImpl::parseModificationTagWithLabel_()
  {
    if (check_(TokenType::HASH)) { Label label = parseLabel_(); InfoTag empty_tag; empty_tag.text = ""; return {std::move(empty_tag), std::move(label)}; }
    ModificationTag tag = parseModificationTag_();
    std::optional<Label> label;
    if (check_(TokenType::HASH)) label = parseLabel_();
    return {std::move(tag), std::move(label)};
  }

  ModificationTag ProFormaParserImpl::parseModificationTag_()
  {
    Token tok = current_();
    if (tok.type == TokenType::PLUS || tok.type == TokenType::MINUS || tok.type == TokenType::NUMBER) return parseMassDelta_();
    if (tok.type != TokenType::IDENTIFIER) error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected modification");
    std::string_view id = tok.text;
    if (id == "Formula") { advance_(); expect_(TokenType::COLON, "':' after Formula"); return parseFormulaTag_(); }
    else if (id == "Glycan") { advance_(); expect_(TokenType::COLON, "':' after Glycan"); return parseGlycanComposition_(); }
    else if (id == "INFO" || id == "info" || id == "Info") { advance_(); expect_(TokenType::COLON, "':' after INFO"); return parseInfoTag_(); }
    else if (id == "Position" || id == "position" || id == "POSITION") { advance_(); expect_(TokenType::COLON, "':' after Position"); return parsePositionConstraint_(); }
    else if (id == "Cation")
    {
      advance_(); expect_(TokenType::COLON, "':' after Cation");
      NamedMod nm; std::string cation_str;
      if (check_(TokenType::IDENTIFIER)) { Token elem = advance_(); cation_str = std::string(elem.text); }
      if (match_(TokenType::LBRACKET)) { cation_str += "["; while (!check_(TokenType::RBRACKET) && !isAtEnd_()) { Token t = advance_(); cation_str += std::string(t.text); } expect_(TokenType::RBRACKET, "']' for cation charge"); cation_str += "]"; }
      nm.name = "Cation:" + cation_str;
      return nm;
    }
    else if (id == "Obs") { advance_(); expect_(TokenType::COLON, "':' after Obs"); MassDelta md = parseMassDelta_(); md.source = MassDelta::Source::OBS; return md; }
    else if (id == "UNIMOD" || id == "MOD" || id == "RESID" || id == "XLMOD" || id == "GNO" || id == "unimod" || id == "mod" || id == "resid" || id == "xlmod" || id == "gno" || id == "xlink") return parseCvAccession_();
    else if (id.size() == 1)
    {
      char prefix = id[0];
      if (prefix == 'U' || prefix == 'M' || prefix == 'R' || prefix == 'X' || prefix == 'G')
      {
        Tokenizer lookahead = tokenizer_;
        if (lookahead.peek().type == TokenType::COLON)
        {
          advance_(); advance_();
          tok = current_();
          if (tok.type == TokenType::PLUS || tok.type == TokenType::MINUS || tok.type == TokenType::NUMBER)
          {
            MassDelta md = parseMassDelta_();
            switch (prefix) { case 'U': md.source = MassDelta::Source::U; break; case 'M': md.source = MassDelta::Source::M; break; case 'R': md.source = MassDelta::Source::R; break; case 'X': md.source = MassDelta::Source::X; break; case 'G': md.source = MassDelta::Source::G; break; }
            return md;
          }
          else return parseNamedMod_(prefix);
        }
      }
    }
    return parseNamedMod_();
  }

  NamedMod ProFormaParserImpl::parseNamedMod_()
  {
    NamedMod nm;
    nm.cv_hint = std::nullopt;
    std::string name;
    int paren_depth = 0, bracket_depth = 0;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::IDENTIFIER) { name += std::string(tok.text); advance_(); }
      else if (tok.type == TokenType::NUMBER) { name += std::string(tok.text); advance_(); }
      else if (tok.type == TokenType::MINUS)
      {
        Tokenizer lookahead = tokenizer_;
        Token next = lookahead.peek();
        if (next.type == TokenType::RANGLE) { name += "->"; advance_(); advance_(); }
        else if (next.type == TokenType::IDENTIFIER || next.type == TokenType::NUMBER || next.type == TokenType::LPAREN || next.type == TokenType::LBRACKET) { name += "-"; advance_(); }
        else if (paren_depth > 0 || bracket_depth > 0) { name += "-"; advance_(); }
        else break;
      }
      else if (tok.type == TokenType::LPAREN) { if (!name.empty() || paren_depth > 0 || bracket_depth > 0) { name += "("; paren_depth++; advance_(); } else break; }
      else if (tok.type == TokenType::RPAREN) { if (paren_depth > 0) { name += ")"; paren_depth--; advance_(); } else break; }
      else if (tok.type == TokenType::LBRACKET) { if (!name.empty() || paren_depth > 0 || bracket_depth > 0) { name += "["; bracket_depth++; advance_(); } else break; }
      else if (tok.type == TokenType::RBRACKET) { if (bracket_depth > 0) { name += "]"; bracket_depth--; advance_(); } else break; }
      else if (tok.type == TokenType::RANGLE) break;
      else break;
    }
    if (name.empty()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected modification name");
    nm.name = String(name);
    return nm;
  }

  NamedMod ProFormaParserImpl::parseNamedMod_(char cv_hint)
  {
    NamedMod nm = parseNamedMod_();
    switch (cv_hint) { case 'U': nm.cv_hint = CvDatabase::UNIMOD; break; case 'M': nm.cv_hint = CvDatabase::MOD; break; case 'R': nm.cv_hint = CvDatabase::RESID; break; case 'X': nm.cv_hint = CvDatabase::XLMOD; break; case 'G': nm.cv_hint = CvDatabase::GNO; break; }
    return nm;
  }

  CvAccession ProFormaParserImpl::parseCvAccession_()
  {
    CvAccession cv;
    Token prefix_tok = expect_(TokenType::IDENTIFIER, "CV prefix");
    std::string_view prefix = prefix_tok.text;
    if (prefix == "UNIMOD" || prefix == "unimod") cv.database = CvDatabase::UNIMOD;
    else if (prefix == "MOD" || prefix == "mod") cv.database = CvDatabase::MOD;
    else if (prefix == "RESID" || prefix == "resid") cv.database = CvDatabase::RESID;
    else if (prefix == "XLMOD" || prefix == "xlmod" || prefix == "xlink") cv.database = CvDatabase::XLMOD;
    else if (prefix == "GNO" || prefix == "gno") cv.database = CvDatabase::GNO;
    else errorAt_(ErrorCode::INVALID_CV_PREFIX, prefix_tok.position, "Invalid CV prefix");
    expect_(TokenType::COLON, "':' after CV prefix");
    if (cv.database == CvDatabase::GNO || cv.database == CvDatabase::RESID || cv.database == CvDatabase::XLMOD)
    {
      std::string accession;
      int bracket_depth = 0;
      while (!isAtEnd_())
      {
        Token tok = current_();
        if (tok.type == TokenType::IDENTIFIER || tok.type == TokenType::NUMBER) { accession += std::string(tok.text); advance_(); }
        else if (tok.type == TokenType::LBRACKET) { accession += "["; bracket_depth++; advance_(); }
        else if (tok.type == TokenType::RBRACKET) { if (bracket_depth > 0) { accession += "]"; bracket_depth--; advance_(); } else break; }
        else break;
      }
      if (accession.empty()) error_(ErrorCode::INVALID_CV_ACCESSION, "Expected accession");
      cv.accession = accession;
    }
    else { Token num = expect_(TokenType::NUMBER, "accession number"); cv.accession = String(num.text); }
    return cv;
  }

  MassDelta ProFormaParserImpl::parseMassDelta_()
  {
    MassDelta md;
    md.source = MassDelta::Source::NONE;
    std::string original;
    Token tok = current_();
    if (tok.type == TokenType::PLUS) { original = "+"; advance_(); tok = current_(); }
    else if (tok.type == TokenType::MINUS) { original = "-"; advance_(); tok = current_(); }
    if (tok.type == TokenType::NUMBER) { original += std::string(tok.text); advance_(); }
    else error_(ErrorCode::INVALID_MASS_VALUE, "Expected mass value");
    md.original_text = original;
    try { md.mass = std::stod(original); }
    catch (const std::exception&) { error_(ErrorCode::INVALID_MASS_VALUE, "Invalid mass value format"); }
    return md;
  }

  FormulaTag ProFormaParserImpl::parseFormulaTag_()
  {
    FormulaTag ft;
    std::string formula;
    int bracket_depth = 0;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::IDENTIFIER) { formula += std::string(tok.text); advance_(); }
      else if (tok.type == TokenType::NUMBER) { formula += std::string(tok.text); advance_(); }
      else if (tok.type == TokenType::LPAREN) { formula += "("; advance_(); }
      else if (tok.type == TokenType::RPAREN) { formula += ")"; advance_(); }
      else if (tok.type == TokenType::LBRACKET) { formula += "["; bracket_depth++; advance_(); }
      else if (tok.type == TokenType::RBRACKET) { if (bracket_depth > 0) { formula += "]"; bracket_depth--; advance_(); } else break; }
      else if (tok.type == TokenType::MINUS)
      {
        Tokenizer lookahead = tokenizer_;
        Token next = lookahead.peek();
        if (bracket_depth > 0 || next.type == TokenType::NUMBER) { formula += "-"; advance_(); }
        else break;
      }
      else if (tok.type == TokenType::COLON)
      {
        Tokenizer lookahead = tokenizer_;
        Token next = lookahead.peek();
        if (next.type == TokenType::IDENTIFIER && next.text == "z")
        {
          advance_(); advance_();
          int sign = 1;
          if (match_(TokenType::PLUS)) sign = 1;
          else if (match_(TokenType::MINUS)) sign = -1;
          Token num = expect_(TokenType::NUMBER, "charge value");
          try { ft.charge = sign * std::stoi(std::string(num.text)); }
          catch (const std::exception&) { errorAt_(ErrorCode::INVALID_CHARGE, num.position, "Invalid charge value"); }
          break;
        }
        else break;
      }
      else break;
    }
    if (formula.empty()) error_(ErrorCode::INVALID_FORMULA, "Empty formula");
    ft.formula_string = formula;
    return ft;
  }

  GlycanComposition ProFormaParserImpl::parseGlycanComposition_()
  {
    GlycanComposition gc;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type != TokenType::IDENTIFIER) break;
      std::string mono_name(tok.text);
      advance_();
      int count = 1;
      if (check_(TokenType::NUMBER))
      {
        Token num = advance_();
        try { count = std::stoi(std::string(num.text)); }
        catch (const std::exception&) { errorAt_(ErrorCode::INVALID_MASS_VALUE, num.position, "Invalid monosaccharide count"); }
      }
      gc.components.emplace_back(String(mono_name), count);
    }
    if (gc.components.empty()) error_(ErrorCode::UNKNOWN_MONOSACCHARIDE, "Empty glycan composition");
    return gc;
  }

  InfoTag ProFormaParserImpl::parseInfoTag_()
  {
    InfoTag it;
    std::string text;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::RBRACKET || tok.type == TokenType::PIPE || tok.type == TokenType::HASH || tok.type == TokenType::COMMA) break;
      text += std::string(tok.text);
      advance_();
    }
    it.text = text;
    return it;
  }

  PositionConstraint ProFormaParserImpl::parsePositionConstraint_()
  {
    PositionConstraint pc;
    while (!isAtEnd_())
    {
      Token tok = current_();
      if (tok.type == TokenType::RBRACKET || tok.type == TokenType::PIPE || tok.type == TokenType::HASH || tok.type == TokenType::RBRACE) break;
      if (tok.type == TokenType::IDENTIFIER && (tok.text == "N" || tok.text == "C"))
      {
        Tokenizer lookahead = createLookahead_();
        lookahead.next();
        Token next1 = lookahead.peek();
        if (next1.type == TokenType::MINUS)
        {
          lookahead.next();
          Token next2 = lookahead.peek();
          if (next2.type == TokenType::IDENTIFIER && next2.text == "term")
          {
            if (tok.text == "N") pc.n_term = true; else pc.c_term = true;
            advance_(); advance_(); advance_();
            if (check_(TokenType::COMMA)) advance_();
            continue;
          }
        }
      }
      if (tok.type == TokenType::COMMA) { advance_(); continue; }
      if (tok.type == TokenType::IDENTIFIER)
      {
        for (char c : tok.text) { if (isAminoAcid_(c)) pc.residues.push_back(c); else errorAt_(ErrorCode::INVALID_AMINO_ACID, tok.position, "Invalid amino acid in position constraint"); }
        advance_();
      }
      else error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected amino acid residues or terminal position after Position:");
    }
    if (pc.residues.empty() && !pc.n_term && !pc.c_term) error_(ErrorCode::UNEXPECTED_CHARACTER, "Position constraint requires at least one residue or terminal position");
    return pc;
  }

  Label ProFormaParserImpl::parseLabel_()
  {
    Label label;
    expect_(TokenType::HASH, "'#'");
    std::string label_str;
    if (check_(TokenType::IDENTIFIER)) { Token id = advance_(); label_str = std::string(id.text); if (check_(TokenType::NUMBER)) { Token num = advance_(); label_str += std::string(num.text); } }
    else if (check_(TokenType::NUMBER)) { Token num = advance_(); label_str = std::string(num.text); }
    else error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected label identifier");
    label.identifier = String(label_str);
    if (label.identifier == "BRANCH") label.type = Label::Type::BRANCH;
    else if (label.identifier.hasPrefix("XL")) label.type = Label::Type::CROSSLINK;
    else label.type = Label::Type::AMBIGUOUS;
    if (match_(TokenType::LPAREN))
    {
      Token score_tok = expect_(TokenType::NUMBER, "score value");
      try { label.score = std::stod(std::string(score_tok.text)); }
      catch (const std::exception&) { errorAt_(ErrorCode::INVALID_MASS_VALUE, score_tok.position, "Invalid score value"); }
      expect_(TokenType::RPAREN, "')'");
    }
    return label;
  }

  std::optional<ChargeState> ProFormaParserImpl::parseChargeState_()
  {
    Token tok = current_();
    if (tok.type == TokenType::LBRACKET) return parseAdductIons_();
    else
    {
      int sign = 1;
      if (tok.type == TokenType::PLUS) { sign = 1; advance_(); tok = current_(); }
      else if (tok.type == TokenType::MINUS) { sign = -1; advance_(); tok = current_(); }
      if (tok.type == TokenType::NUMBER)
      {
        int charge;
        try { charge = sign * std::stoi(std::string(tok.text)); }
        catch (const std::exception&) { errorAt_(ErrorCode::INVALID_CHARGE, tok.position, "Invalid charge value"); }
        advance_();
        return charge;
      }
      else error_(ErrorCode::INVALID_CHARGE, "Expected charge value");
    }
    return std::nullopt;
  }

  std::vector<AdductIon> ProFormaParserImpl::parseAdductIons_()
  {
    std::vector<AdductIon> adducts;
    expect_(TokenType::LBRACKET, "'['");
    adducts.push_back(parseAdductIon_());
    while (match_(TokenType::COMMA)) adducts.push_back(parseAdductIon_());
    expect_(TokenType::RBRACKET, "']'");
    return adducts;
  }

  AdductIon ProFormaParserImpl::parseAdductIon_()
  {
    AdductIon adduct;
    std::string formula;
    while (!isAtEnd_())
    {
      if (check_(TokenType::COLON)) { Tokenizer lookahead = createLookahead_(); lookahead.next(); Token next = lookahead.peek(); if (next.type == TokenType::IDENTIFIER && next.text == "z") break; }
      Token tok = current_();
      if (tok.type == TokenType::RBRACKET || tok.type == TokenType::COMMA) error_(ErrorCode::UNEXPECTED_CHARACTER, "Unexpected token in adduct formula, expected ':z'");
      formula += std::string(tok.text);
      advance_();
    }
    if (formula.empty()) error_(ErrorCode::UNEXPECTED_CHARACTER, "Expected adduct formula");
    adduct.formula = String(formula);
    expect_(TokenType::COLON, "':'");
    expect_(TokenType::IDENTIFIER, "'z'");
    int sign = 1;
    if (match_(TokenType::PLUS)) sign = 1;
    else if (match_(TokenType::MINUS)) sign = -1;
    Token charge_tok = expect_(TokenType::NUMBER, "charge value");
    try { adduct.charge = sign * std::stoi(std::string(charge_tok.text)); }
    catch (const std::exception&) { errorAt_(ErrorCode::INVALID_CHARGE, charge_tok.position, "Invalid adduct charge value"); }
    if (match_(TokenType::CARET))
    {
      Token occ = expect_(TokenType::NUMBER, "occurrence count");
      try { adduct.occurrence = std::stoi(std::string(occ.text)); }
      catch (const std::exception&) { errorAt_(ErrorCode::INVALID_MASS_VALUE, occ.position, "Invalid adduct occurrence count"); }
    }
    return adduct;
  }

  Token ProFormaParserImpl::current_() { if (!has_current_) { current_token_ = tokenizer_.next(); has_current_ = true; } return current_token_; }
  Token ProFormaParserImpl::advance_() { Token tok = current_(); has_current_ = false; return tok; }
  bool ProFormaParserImpl::check_(TokenType type) { return current_().type == type; }
  bool ProFormaParserImpl::match_(TokenType type) { if (check_(type)) { advance_(); return true; } return false; }
  Token ProFormaParserImpl::expect_(TokenType type, const char* expected_desc) { Token tok = current_(); if (tok.type != type) errorAt_(ErrorCode::UNEXPECTED_CHARACTER, tok.position, (std::string("Expected ") + expected_desc).c_str()); return advance_(); }
  bool ProFormaParserImpl::isAtEnd_() { return current_().type == TokenType::END; }
  void ProFormaParserImpl::error_(ErrorCode code, const char* message) { errorAt_(code, current_().position, message); }
  void ProFormaParserImpl::errorAt_(ErrorCode code, size_t pos, const char* message) { throw ProForma::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, code, pos, input_, message); }
  bool ProFormaParserImpl::isAminoAcid_(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }

} // namespace detail

//============================================================================
// Public ProFormaParser static methods
//============================================================================

ProForma::Peptidoform ProForma::parse(const String& input)
{
  detail::ProFormaParserImpl parser(input);
  return parser.parsePeptidoform();
}

ProForma::PeptidoformIon ProForma::parseIon(const String& input)
{
  detail::ProFormaParserImpl parser(input);
  return parser.parsePeptidoformIon();
}

String ProForma::toString(const Peptidoform& pf, WriteMode mode)
{
  return detail::ProFormaWriter::toString(pf, mode);
}

String ProForma::toString(const PeptidoformIon& pfi, WriteMode mode)
{
  return detail::ProFormaWriter::toString(pfi, mode);
}

//============================================================================
// Modification resolution helpers
//============================================================================
namespace
{
  // Type aliases for convenience within anonymous namespace
  using Peptidoform = ProForma::Peptidoform;
  using PeptidoformIon = ProForma::PeptidoformIon;
  using Modification = ProForma::Modification;
  using ModificationTag = ProForma::ModificationTag;
  using CvAccession = ProForma::CvAccession;
  using NamedMod = ProForma::NamedMod;
  using MassDelta = ProForma::MassDelta;
  using FormulaTag = ProForma::FormulaTag;
  using GlycanComposition = ProForma::GlycanComposition;
  using InfoTag = ProForma::InfoTag;
  using PositionConstraint = ProForma::PositionConstraint;
  using Label = ProForma::Label;
  using UnlocalisedMod = ProForma::UnlocalisedMod;
  using LabileModification = ProForma::LabileModification;
  using SequenceSection = ProForma::SequenceSection;
  using SequenceElement = ProForma::SequenceElement;
  using AmbiguousRegion = ProForma::AmbiguousRegion;
  using ModifiedRange = ProForma::ModifiedRange;
  using GlobalModEntry = ProForma::GlobalModEntry;
  using GlobalModification = ProForma::GlobalModification;
  using IsotopeReplacement = ProForma::IsotopeReplacement;
  using CvDatabase = ProForma::CvDatabase;
  using ConversionIssue = ProForma::ConversionIssue;
  using ConversionIssueType = ProForma::ConversionIssueType;
  using ConversionPolicy = ProForma::ConversionPolicy;

  // Helper to resolve a single modification tag to a ResidueModification
  const ResidueModification* resolveModificationTag_(
    const ModificationTag& tag,
    char residue = '\0',
    ResidueModification::TermSpecificity term_spec = ResidueModification::NUMBER_OF_TERM_SPECIFICITY)
  {
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    return std::visit([&](auto&& arg) -> const ResidueModification* {
      using T = std::decay_t<decltype(arg)>;

      if constexpr (std::is_same_v<T, CvAccession>)
      {
        String full_accession;
        switch (arg.database)
        {
          case CvDatabase::UNIMOD: full_accession = "UNIMOD:" + arg.accession; break;
          case CvDatabase::MOD:    full_accession = "MOD:" + arg.accession; break;
          case CvDatabase::RESID:  full_accession = "RESID:" + arg.accession; break;
          case CvDatabase::XLMOD:  full_accession = "XLMOD:" + arg.accession; break;
          case CvDatabase::GNO:    full_accession = "GNO:" + arg.accession; break;
        }
        try
        {
          String residue_str = (residue != '\0') ? String(1, residue) : "";
          return mod_db->getModification(full_accession, residue_str, term_spec);
        }
        catch (const Exception::ElementNotFound&) { return nullptr; }
      }
      else if constexpr (std::is_same_v<T, NamedMod>)
      {
        bool multiple_matches = false;
        String residue_str = (residue != '\0') ? String(1, residue) : "";
        return mod_db->searchModificationsFast(arg.name, multiple_matches, residue_str, term_spec);
      }
      else if constexpr (std::is_same_v<T, MassDelta>)
      {
        String residue_str = (residue != '\0') ? String(1, residue) : "";
        return mod_db->getBestModificationByDiffMonoMass(arg.mass, 0.01, residue_str, term_spec);
      }
      else if constexpr (std::is_same_v<T, FormulaTag>) { return nullptr; }
      else if constexpr (std::is_same_v<T, GlycanComposition>) { return nullptr; }
      else if constexpr (std::is_same_v<T, InfoTag>) { return nullptr; }
      else { return nullptr; }
    }, tag);
  }

  void resolveModification_(Modification& mod, char residue, ResidueModification::TermSpecificity term_spec)
  {
    if (mod.alternatives.empty()) return;
    const auto& [tag, label] = mod.alternatives[0];
    mod.resolved_mod = resolveModificationTag_(tag, residue, term_spec);
  }

  // Helper to get modification mass from a Modification struct
  std::pair<bool, double> getModificationMass_(const Modification& mod)
  {
    if (mod.resolved_mod != nullptr) return {true, mod.resolved_mod->getDiffMonoMass()};
    if (mod.alternatives.empty()) return {false, 0.0};

    const auto& tag = mod.alternatives[0].first;

    if (const auto* md = std::get_if<MassDelta>(&tag)) return {true, md->mass};
    if (const auto* ft = std::get_if<FormulaTag>(&tag))
    {
      try { EmpiricalFormula formula(ft->formula_string); return {true, formula.getMonoWeight()}; }
      catch (...) { return {false, 0.0}; }
    }
    if (std::holds_alternative<InfoTag>(tag) || std::holds_alternative<PositionConstraint>(tag)) return {true, 0.0};
    return {false, 0.0};
  }

  void checkModificationForMass_(const Modification& mod, size_t position, std::vector<ConversionIssue>& issues)
  {
    auto [has_mass, mass] = getModificationMass_(mod);
    (void)mass;
    if (!has_mass)
    {
      issues.push_back({ConversionIssueType::UNRESOLVED_MOD,
        "Modification at position " + std::to_string(position) + " has no resolvable mass", position});
    }
  }

  double calculateChainMass_(const Peptidoform& pf_resolved, std::set<String>& counted_crosslinks)
  {
    double mass = 0.0;

    auto addModMass = [&](const Modification& mod) {
      if (!mod.alternatives.empty() && mod.alternatives[0].second.has_value())
      {
        const auto& label = mod.alternatives[0].second.value();
        if (label.type == Label::Type::CROSSLINK)
        {
          if (counted_crosslinks.count(label.identifier) > 0) return;
          counted_crosslinks.insert(label.identifier);
        }
      }
      auto [has_mass, mod_mass] = getModificationMass_(mod);
      if (has_mass) mass += mod_mass;
    };

    for (const auto& section : pf_resolved.sequence)
    {
      if (const auto* elem = std::get_if<SequenceElement>(&section))
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(elem->amino_acid);
        mass += res->getMonoWeight(Residue::Internal);
        for (const auto& mod : elem->modifications) addModMass(mod);
      }
      else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
      {
        if (!region->elements.empty())
        {
          const Residue* res = ResidueDB::getInstance()->getResidue(region->elements[0].amino_acid);
          mass += res->getMonoWeight(Residue::Internal);
          for (const auto& mod : region->elements[0].modifications) addModMass(mod);
        }
      }
      else if (const auto* range = std::get_if<ModifiedRange>(&section))
      {
        for (const auto& elem : range->elements)
        {
          const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
          mass += res->getMonoWeight(Residue::Internal);
        }
        for (const auto& mod : range->modifications) addModMass(mod);
      }
    }

    mass += EmpiricalFormula("H2O").getMonoWeight();

    for (const auto& mod : pf_resolved.n_term_mods) addModMass(mod);
    for (const auto& mod : pf_resolved.c_term_mods) addModMass(mod);

    for (const auto& um : pf_resolved.unlocalised_mods)
    {
      int occurrence = um.occurrence.value_or(1);
      for (const auto& mod : um.modifications)
      {
        auto [has_mass, mod_mass] = getModificationMass_(mod);
        if (has_mass) mass += mod_mass * occurrence;
      }
    }

    for (const auto& lm : pf_resolved.labile_mods)
    {
      auto [has_mass, mod_mass] = getModificationMass_(lm.modification);
      if (has_mass) mass += mod_mass;
    }

    for (const auto& entry : pf_resolved.global_mods)
    {
      if (const auto* gm = std::get_if<GlobalModification>(&entry))
      {
        auto [has_mass, mod_mass] = getModificationMass_(gm->modification);
        if (has_mass)
        {
          int count = 0;
          for (const auto& section : pf_resolved.sequence)
          {
            if (const auto* elem = std::get_if<SequenceElement>(&section))
            {
              for (const String& loc : gm->locations)
              {
                if (loc.size() == 1 && elem->amino_acid == loc[0]) { ++count; break; }
              }
            }
          }
          mass += mod_mass * count;
        }
      }
    }

    return mass;
  }

  // Spectrum generation helpers
  std::tuple<bool, size_t, double, String> findCrossLink(const Peptidoform& chain)
  {
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
              double mass = 0.0;
              if (const auto* md = std::get_if<MassDelta>(&mod.alternatives[0].first)) mass = md->mass;
              else if (mod.resolved_mod != nullptr) mass = mod.resolved_mod->getDiffMonoMass();
              return {true, position, mass, label.identifier};
            }
          }
        }
        ++position;
      }
    }
    return {false, 0, 0.0, ""};
  }

  std::vector<ConversionIssue> collectPeptidoformSpectrumIssues(const Peptidoform& pf)
  {
    std::vector<ConversionIssue> issues;
    if (!ProForma::isRepresentableAsAASequence(pf))
    {
      auto conv_issues = ProForma::getAASequenceConversionIssues(pf);
      for (auto& issue : conv_issues) issues.push_back(std::move(issue));
    }
    return issues;
  }

  std::vector<ConversionIssue> collectPeptidoformIonSpectrumIssues(const PeptidoformIon& pfi)
  {
    std::vector<ConversionIssue> issues;
    if (pfi.chains.empty())
    {
      issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "No peptide chains to fragment", 0});
      return issues;
    }
    if (pfi.is_chimeric)
    {
      issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
        "Theoretical spectrum generation not supported for chimeric spectra.", 0});
      return issues;
    }
    if (pfi.chains.size() == 1) return collectPeptidoformSpectrumIssues(pfi.chains[0]);
    if (pfi.chains.size() != 2)
    {
      issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
        "Only two-chain cross-links are currently supported for spectrum generation", 0});
      return issues;
    }

    auto [alpha_found, alpha_pos, alpha_mass, alpha_label] = findCrossLink(pfi.chains[0]);
    auto [beta_found, beta_pos, beta_mass, beta_label] = findCrossLink(pfi.chains[1]);

    if (!alpha_found || !beta_found)
    {
      issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "Cross-link label not found in both chains", 0});
      return issues;
    }
    if (alpha_label != beta_label)
    {
      issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "Cross-link labels don't match between chains", 0});
      return issues;
    }
    return issues;
  }
} // anonymous namespace

//============================================================================
// Modification resolution
//============================================================================

void ProForma::resolveModifications(Peptidoform& pf)
{
  for (auto& section : pf.sequence)
  {
    if (auto* elem = std::get_if<SequenceElement>(&section))
    {
      for (auto& mod : elem->modifications)
        resolveModification_(mod, elem->amino_acid, ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    }
    else if (auto* region = std::get_if<AmbiguousRegion>(&section))
    {
      for (auto& elem : region->elements)
        for (auto& mod : elem.modifications)
          resolveModification_(mod, elem.amino_acid, ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    }
    else if (auto* range = std::get_if<ModifiedRange>(&section))
    {
      for (auto& mod : range->modifications)
        resolveModification_(mod, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
    }
  }

  for (auto& mod : pf.n_term_mods) resolveModification_(mod, '\0', ResidueModification::N_TERM);
  for (auto& mod : pf.c_term_mods) resolveModification_(mod, '\0', ResidueModification::C_TERM);

  for (auto& um : pf.unlocalised_mods)
    for (auto& mod : um.modifications)
      resolveModification_(mod, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

  for (auto& lm : pf.labile_mods)
    resolveModification_(lm.modification, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);

  for (auto& entry : pf.global_mods)
    if (auto* gm = std::get_if<GlobalModification>(&entry))
      resolveModification_(gm->modification, '\0', ResidueModification::NUMBER_OF_TERM_SPECIFICITY);
}

//============================================================================
// AASequence conversion
//============================================================================

std::vector<ProForma::ConversionIssue> ProForma::getAASequenceConversionIssues(const Peptidoform& pf)
{
  std::vector<ConversionIssue> issues;
  Peptidoform pf_copy = pf;
  resolveModifications(pf_copy);

  if (!pf_copy.unlocalised_mods.empty())
    issues.push_back({ConversionIssueType::UNLOCALISED_MOD, "Peptidoform contains unlocalised modifications", SIZE_MAX});

  if (!pf_copy.labile_mods.empty())
    issues.push_back({ConversionIssueType::LABILE_MOD, "Peptidoform contains labile modifications", SIZE_MAX});

  for (const auto& entry : pf_copy.global_mods)
    if (std::holds_alternative<GlobalModification>(entry))
    {
      issues.push_back({ConversionIssueType::GLOBAL_MOD, "Peptidoform contains global modifications", SIZE_MAX});
      break;
    }

  size_t position = 0;
  for (const auto& section : pf_copy.sequence)
  {
    if (std::holds_alternative<AmbiguousRegion>(section))
      issues.push_back({ConversionIssueType::AMBIGUOUS_REGION,
        "Peptidoform contains ambiguous region at position " + std::to_string(position), position});
    else if (std::holds_alternative<ModifiedRange>(section))
      issues.push_back({ConversionIssueType::MODIFIED_RANGE,
        "Peptidoform contains modified range at position " + std::to_string(position), position});
    else if (const auto* elem = std::get_if<SequenceElement>(&section))
    {
      for (const auto& mod : elem->modifications)
      {
        if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
        {
          const auto& [tag, label] = mod.alternatives[0];
          bool is_empty_info = std::holds_alternative<InfoTag>(tag) && std::get<InfoTag>(tag).text.empty();
          if (!is_empty_info)
            issues.push_back({ConversionIssueType::UNRESOLVED_MOD,
              "Modification at position " + std::to_string(position) + " could not be resolved", position});
        }
        if (mod.alternatives.size() > 1)
          issues.push_back({ConversionIssueType::ALTERNATIVE_MODS,
            "Modification at position " + std::to_string(position) + " has multiple alternatives", position});

        for (const auto& [tag, label] : mod.alternatives)
          if (label.has_value() && label->type == Label::Type::CROSSLINK)
          {
            issues.push_back({ConversionIssueType::CROSS_LINK,
              "Modification at position " + std::to_string(position) + " is part of a cross-link", position});
            break;
          }
      }
    }
    position++;
  }

  for (const auto& mod : pf_copy.n_term_mods)
  {
    if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
    {
      const auto& [tag, label] = mod.alternatives[0];
      bool is_empty_info = std::holds_alternative<InfoTag>(tag) && std::get<InfoTag>(tag).text.empty();
      if (!is_empty_info)
        issues.push_back({ConversionIssueType::UNRESOLVED_MOD, "N-terminal modification could not be resolved", SIZE_MAX});
    }
    if (mod.alternatives.size() > 1)
      issues.push_back({ConversionIssueType::ALTERNATIVE_MODS, "N-terminal modification has multiple alternatives", SIZE_MAX});
  }

  for (const auto& mod : pf_copy.c_term_mods)
  {
    if (mod.resolved_mod == nullptr && !mod.alternatives.empty())
    {
      const auto& [tag, label] = mod.alternatives[0];
      bool is_empty_info = std::holds_alternative<InfoTag>(tag) && std::get<InfoTag>(tag).text.empty();
      if (!is_empty_info)
        issues.push_back({ConversionIssueType::UNRESOLVED_MOD, "C-terminal modification could not be resolved", SIZE_MAX});
    }
    if (mod.alternatives.size() > 1)
      issues.push_back({ConversionIssueType::ALTERNATIVE_MODS, "C-terminal modification has multiple alternatives", SIZE_MAX});
  }

  return issues;
}

bool ProForma::isRepresentableAsAASequence(const Peptidoform& pf)
{
  return getAASequenceConversionIssues(pf).empty();
}

AASequence ProForma::toAASequence(const Peptidoform& pf, ConversionPolicy policy)
{
  std::vector<ConversionIssue> issues = getAASequenceConversionIssues(pf);

  if (policy == ConversionPolicy::FAIL_ON_LOSS && !issues.empty())
  {
    std::string error_msg = "Cannot convert Peptidoform to AASequence: ";
    for (const auto& issue : issues) error_msg += issue.description + "; ";
    throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
  }

  Peptidoform pf_copy = pf;
  resolveModifications(pf_copy);

  std::string unmod_seq;
  for (const auto& section : pf_copy.sequence)
  {
    if (const auto* elem = std::get_if<SequenceElement>(&section)) unmod_seq += elem->amino_acid;
    else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
    {
      if (!region->elements.empty()) unmod_seq += region->elements[0].amino_acid;
    }
    else if (const auto* range = std::get_if<ModifiedRange>(&section))
      for (const auto& elem : range->elements) unmod_seq += elem.amino_acid;
  }

  AASequence seq = AASequence::fromString(unmod_seq);

  size_t seq_pos = 0;
  for (const auto& section : pf_copy.sequence)
  {
    if (const auto* elem = std::get_if<SequenceElement>(&section))
    {
      for (const auto& mod : elem->modifications)
      {
        if (mod.resolved_mod != nullptr) seq.setModification(seq_pos, mod.resolved_mod);
        else if (policy == ConversionPolicy::FAIL_ON_LOSS)
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Unresolved modification at position " + std::to_string(seq_pos));
      }
      seq_pos++;
    }
    else if (std::holds_alternative<AmbiguousRegion>(section)) seq_pos++;
    else if (const auto* range = std::get_if<ModifiedRange>(&section)) seq_pos += range->elements.size();
  }

  if (!pf_copy.n_term_mods.empty() && pf_copy.n_term_mods[0].resolved_mod != nullptr)
    seq.setNTerminalModification(pf_copy.n_term_mods[0].resolved_mod);

  if (!pf_copy.c_term_mods.empty() && pf_copy.c_term_mods[0].resolved_mod != nullptr)
    seq.setCTerminalModification(pf_copy.c_term_mods[0].resolved_mod);

  return seq;
}

ProForma::Peptidoform ProForma::fromAASequence(const AASequence& seq)
{
  Peptidoform pf;

  for (Size i = 0; i < seq.size(); ++i)
  {
    SequenceElement elem;
    elem.amino_acid = seq[i].getOneLetterCode()[0];

    if (seq[i].isModified())
    {
      const ResidueModification* mod = seq[i].getModification();
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
        elem.modifications.push_back(std::move(proforma_mod));
      }
    }
    pf.sequence.push_back(std::move(elem));
  }

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

//============================================================================
// Mass calculation methods
//============================================================================

std::vector<ProForma::ConversionIssue> ProForma::getMassCalculationIssues(const Peptidoform& pf)
{
  std::vector<ConversionIssue> issues;
  Peptidoform pf_copy = pf;
  resolveModifications(pf_copy);

  size_t position = 0;
  for (const auto& section : pf_copy.sequence)
  {
    if (const auto* elem = std::get_if<SequenceElement>(&section))
    {
      const Residue* res = ResidueDB::getInstance()->getResidue(elem->amino_acid);
      if (res == nullptr)
        issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
          String("Unknown amino acid '") + elem->amino_acid + "' at position " + String(position), position});
      for (const auto& mod : elem->modifications) checkModificationForMass_(mod, position, issues);
      ++position;
    }
    else if (const auto* region = std::get_if<AmbiguousRegion>(&section))
    {
      std::set<double> masses;
      for (const auto& elem : region->elements)
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
        if (res != nullptr) masses.insert(res->getMonoWeight(Residue::Internal));
        else issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
          String("Unknown amino acid '") + elem.amino_acid + "' in ambiguous region", position});
      }
      if (masses.size() > 1)
        issues.push_back({ConversionIssueType::AMBIGUOUS_REGION,
          "Ambiguous region contains amino acids with different masses", position});
      ++position;
    }
    else if (const auto* range = std::get_if<ModifiedRange>(&section))
    {
      for (const auto& elem : range->elements)
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(elem.amino_acid);
        if (res == nullptr)
          issues.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
            String("Unknown amino acid '") + elem.amino_acid + "' in range", position});
        ++position;
      }
      for (const auto& mod : range->modifications)
        checkModificationForMass_(mod, position - range->elements.size(), issues);
    }
  }

  for (const auto& mod : pf_copy.n_term_mods) checkModificationForMass_(mod, 0, issues);
  for (const auto& mod : pf_copy.c_term_mods) checkModificationForMass_(mod, position > 0 ? position - 1 : 0, issues);
  for (const auto& um : pf_copy.unlocalised_mods)
    for (const auto& mod : um.modifications) checkModificationForMass_(mod, SIZE_MAX, issues);
  for (const auto& lm : pf_copy.labile_mods) checkModificationForMass_(lm.modification, SIZE_MAX, issues);
  for (const auto& entry : pf_copy.global_mods)
    if (const auto* gm = std::get_if<GlobalModification>(&entry)) checkModificationForMass_(gm->modification, SIZE_MAX, issues);

  return issues;
}

std::vector<ProForma::ConversionIssue> ProForma::getMassCalculationIssues(const PeptidoformIon& pfi)
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

bool ProForma::canCalculateMass(const Peptidoform& pf) { return getMassCalculationIssues(pf).empty(); }
bool ProForma::canCalculateMass(const PeptidoformIon& pfi) { return getMassCalculationIssues(pfi).empty(); }

double ProForma::getMonoWeight(const Peptidoform& pf)
{
  auto issues = getMassCalculationIssues(pf);
  if (!issues.empty())
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate mass: " + issues[0].description, "");

  Peptidoform pf_copy = pf;
  resolveModifications(pf_copy);
  std::set<String> counted_crosslinks;
  return calculateChainMass_(pf_copy, counted_crosslinks);
}

double ProForma::getMonoWeight(const PeptidoformIon& pfi)
{
  auto issues = getMassCalculationIssues(pfi);
  if (!issues.empty())
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate mass: " + issues[0].description, "");

  if (pfi.chains.empty()) return 0.0;

  if (pfi.is_chimeric)
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate single mass for chimeric spectra.", "");

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

double ProForma::getMZ(const PeptidoformIon& pfi)
{
  if (!pfi.charge.has_value())
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate m/z: no charge state specified", "");

  int charge = 0;
  if (const int* simple_charge = std::get_if<int>(&pfi.charge.value())) charge = *simple_charge;
  else if (const auto* adducts = std::get_if<std::vector<AdductIon>>(&pfi.charge.value()))
  {
    for (const auto& adduct : *adducts) charge += adduct.charge * adduct.occurrence.value_or(1);
  }

  if (charge == 0)
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate m/z: charge state is zero", "");

  double mass = getMonoWeight(pfi);
  return (mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
}

double ProForma::getMZ(const Peptidoform& pf, int charge)
{
  if (charge == 0)
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "Cannot calculate m/z: charge state is zero", "");

  double mass = getMonoWeight(pf);
  return (mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
}

std::optional<double> ProForma::tryGetMonoWeight(const Peptidoform& pf)
{
  std::vector<ConversionIssue> issues;
  return tryGetMonoWeight(pf, issues);
}

std::optional<double> ProForma::tryGetMonoWeight(const Peptidoform& pf, std::vector<ConversionIssue>& issues_out)
{
  issues_out.clear();
  Peptidoform pf_copy = pf;
  resolveModifications(pf_copy);
  issues_out = getMassCalculationIssues(pf_copy);
  if (!issues_out.empty()) return std::nullopt;
  std::set<String> counted_crosslinks;
  return calculateChainMass_(pf_copy, counted_crosslinks);
}

std::optional<double> ProForma::tryGetMonoWeight(const PeptidoformIon& pfi)
{
  std::vector<ConversionIssue> issues;
  return tryGetMonoWeight(pfi, issues);
}

std::optional<double> ProForma::tryGetMonoWeight(const PeptidoformIon& pfi, std::vector<ConversionIssue>& issues_out)
{
  issues_out.clear();
  if (pfi.chains.empty()) return 0.0;

  if (pfi.is_chimeric)
  {
    issues_out.push_back({ConversionIssueType::UNSUPPORTED_FEATURE,
      "Cannot calculate single mass for chimeric spectra.", 0});
    return std::nullopt;
  }

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
  if (!issues_out.empty()) return std::nullopt;

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

std::optional<double> ProForma::tryGetMZ(const Peptidoform& pf, int charge)
{
  std::vector<ConversionIssue> issues;
  return tryGetMZ(pf, charge, issues);
}

std::optional<double> ProForma::tryGetMZ(const Peptidoform& pf, int charge, std::vector<ConversionIssue>& issues_out)
{
  issues_out.clear();
  if (charge == 0)
  {
    issues_out.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "Charge state is zero", 0});
    return std::nullopt;
  }
  auto mass = tryGetMonoWeight(pf, issues_out);
  if (!mass.has_value()) return std::nullopt;
  return (*mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
}

std::optional<double> ProForma::tryGetMZ(const PeptidoformIon& pfi)
{
  std::vector<ConversionIssue> issues;
  return tryGetMZ(pfi, issues);
}

std::optional<double> ProForma::tryGetMZ(const PeptidoformIon& pfi, std::vector<ConversionIssue>& issues_out)
{
  issues_out.clear();
  if (!pfi.charge.has_value())
  {
    issues_out.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "No charge state specified", 0});
    return std::nullopt;
  }

  int charge = 0;
  if (const int* simple_charge = std::get_if<int>(&pfi.charge.value())) charge = *simple_charge;
  else if (const auto* adducts = std::get_if<std::vector<AdductIon>>(&pfi.charge.value()))
  {
    for (const auto& adduct : *adducts) charge += adduct.charge * adduct.occurrence.value_or(1);
  }

  if (charge == 0)
  {
    issues_out.push_back({ConversionIssueType::UNSUPPORTED_FEATURE, "Charge state is zero", 0});
    return std::nullopt;
  }

  auto mass = tryGetMonoWeight(pfi, issues_out);
  if (!mass.has_value()) return std::nullopt;
  return (*mass + charge * Constants::PROTON_MASS_U) / std::abs(charge);
}

//============================================================================
// Theoretical spectrum generation
//============================================================================

bool ProForma::canGenerateSpectrum(const Peptidoform& pf) { return collectPeptidoformSpectrumIssues(pf).empty(); }
bool ProForma::canGenerateSpectrum(const PeptidoformIon& pfi) { return collectPeptidoformIonSpectrumIssues(pfi).empty(); }

std::vector<ProForma::ConversionIssue> ProForma::getSpectrumGenerationIssues(const Peptidoform& pf)
{
  return collectPeptidoformSpectrumIssues(pf);
}

std::vector<ProForma::ConversionIssue> ProForma::getSpectrumGenerationIssues(const PeptidoformIon& pfi)
{
  return collectPeptidoformIonSpectrumIssues(pfi);
}

MSSpectrum ProForma::generateSpectrum(
  const Peptidoform& pf,
  int min_charge,
  int max_charge,
  const std::string& ion_types,
  bool add_losses,
  bool add_metainfo)
{
  auto issues = collectPeptidoformSpectrumIssues(pf);
  if (!issues.empty())
  {
    std::string error_msg = "Spectrum generation failed: ";
    for (const auto& issue : issues) error_msg += issue.description + "; ";
    throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
  }

  AASequence seq = toAASequence(pf, ConversionPolicy::FAIL_ON_LOSS);

  TheoreticalSpectrumGenerator generator;
  Param param = generator.getParameters();
  param.setValue("add_a_ions", ion_types.find('a') != std::string::npos ? "true" : "false");
  param.setValue("add_b_ions", ion_types.find('b') != std::string::npos ? "true" : "false");
  param.setValue("add_c_ions", ion_types.find('c') != std::string::npos ? "true" : "false");
  param.setValue("add_x_ions", ion_types.find('x') != std::string::npos ? "true" : "false");
  param.setValue("add_y_ions", ion_types.find('y') != std::string::npos ? "true" : "false");
  param.setValue("add_z_ions", ion_types.find('z') != std::string::npos ? "true" : "false");
  param.setValue("add_precursor_peaks", ion_types.find('M') != std::string::npos ? "true" : "false");
  param.setValue("add_abundant_immonium_ions", ion_types.find('I') != std::string::npos ? "true" : "false");
  param.setValue("add_losses", add_losses ? "true" : "false");
  param.setValue("add_metainfo", add_metainfo ? "true" : "false");
  generator.setParameters(param);

  MSSpectrum spectrum;
  generator.getSpectrum(spectrum, seq, min_charge, max_charge);
  return spectrum;
}

MSSpectrum ProForma::generateSpectrum(
  const PeptidoformIon& pfi,
  int min_charge,
  int max_charge,
  const std::string& ion_types,
  bool add_losses,
  bool add_metainfo)
{
  auto issues = collectPeptidoformIonSpectrumIssues(pfi);
  if (!issues.empty())
  {
    std::string error_msg = "Spectrum generation failed: ";
    for (const auto& issue : issues) error_msg += issue.description + "; ";
    throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
  }

  if (pfi.chains.size() == 1)
    return generateSpectrum(pfi.chains[0], min_charge, max_charge, ion_types, add_losses, add_metainfo);

  const Peptidoform& alpha_pf = pfi.chains[0];
  const Peptidoform& beta_pf = pfi.chains[1];

  auto [alpha_found, alpha_pos, alpha_mass, alpha_label] = findCrossLink(alpha_pf);
  auto [beta_found, beta_pos, beta_mass, beta_label] = findCrossLink(beta_pf);

  AASequence alpha_seq = toAASequence(alpha_pf, ConversionPolicy::BEST_EFFORT);
  AASequence beta_seq = toAASequence(beta_pf, ConversionPolicy::BEST_EFFORT);

  double linker_mass = (alpha_mass > 0.001) ? alpha_mass : beta_mass;

  OPXLDataStructs::ProteinProteinCrossLink crosslink;
  crosslink.alpha = &alpha_seq;
  crosslink.beta = &beta_seq;
  crosslink.cross_link_position = {static_cast<SignedSize>(alpha_pos), static_cast<SignedSize>(beta_pos)};
  crosslink.cross_linker_mass = linker_mass;
  crosslink.cross_linker_name = alpha_label;

  TheoreticalSpectrumGeneratorXLMS generator;
  Param param = generator.getParameters();
  param.setValue("add_a_ions", ion_types.find('a') != std::string::npos ? "true" : "false");
  param.setValue("add_b_ions", ion_types.find('b') != std::string::npos ? "true" : "false");
  param.setValue("add_c_ions", ion_types.find('c') != std::string::npos ? "true" : "false");
  param.setValue("add_x_ions", ion_types.find('x') != std::string::npos ? "true" : "false");
  param.setValue("add_y_ions", ion_types.find('y') != std::string::npos ? "true" : "false");
  param.setValue("add_z_ions", ion_types.find('z') != std::string::npos ? "true" : "false");
  param.setValue("add_precursor_peaks", ion_types.find('M') != std::string::npos ? "true" : "false");
  param.setValue("add_losses", add_losses ? "true" : "false");
  param.setValue("add_metainfo", add_metainfo ? "true" : "false");
  generator.setParameters(param);

  MSSpectrum spectrum;
  generator.getXLinkIonSpectrum(spectrum, crosslink, true, min_charge, max_charge);
  generator.getXLinkIonSpectrum(spectrum, crosslink, false, min_charge, max_charge);
  spectrum.sortByPosition();

  return spectrum;
}

} // namespace OpenMS

