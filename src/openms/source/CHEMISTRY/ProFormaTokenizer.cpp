// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaTokenizer.h>

#include <algorithm>

namespace OpenMS
{

  ProFormaTokenizer::ProFormaTokenizer(std::string_view input)
    : input_(input),
      pos_(0),
      peeked_(std::nullopt)
  {
  }

  ProFormaTokenizer::Token ProFormaTokenizer::next()
  {
    // If we have a peeked token, consume and return it
    if (peeked_.has_value())
    {
      Token token = *peeked_;
      peeked_.reset();
      return token;
    }

    return scanToken_();
  }

  ProFormaTokenizer::Token ProFormaTokenizer::peek()
  {
    // If we already have a peeked token, return it
    if (peeked_.has_value())
    {
      return *peeked_;
    }

    // Scan the next token and cache it
    peeked_ = scanToken_();
    return *peeked_;
  }

  bool ProFormaTokenizer::hasMore() const
  {
    // If we have a peeked token, check if it's END
    if (peeked_.has_value())
    {
      return peeked_->type != TokenType::END;
    }

    // Otherwise check if we're at the end of input
    return pos_ < input_.size();
  }

  size_t ProFormaTokenizer::position() const
  {
    return pos_;
  }

  std::string_view ProFormaTokenizer::getContext(size_t pos, size_t before, size_t after) const
  {
    // Calculate start position (don't go before beginning)
    size_t start = (pos > before) ? (pos - before) : 0;

    // Calculate end position (don't go past end)
    size_t end = std::min(pos + after, input_.size());

    return input_.substr(start, end - start);
  }

  const char* ProFormaTokenizer::tokenTypeName(TokenType type)
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

  ProFormaTokenizer::Token ProFormaTokenizer::scanToken_()
  {
    // Check for end of input
    if (isAtEnd_())
    {
      return Token{TokenType::END, std::string_view{}, pos_};
    }

    char c = current_();
    size_t start_pos = pos_;

    // Single-character tokens
    switch (c)
    {
      case '[':
        advance_();
        return Token{TokenType::LBRACKET, input_.substr(start_pos, 1), start_pos};
      case ']':
        advance_();
        return Token{TokenType::RBRACKET, input_.substr(start_pos, 1), start_pos};
      case '(':
        advance_();
        return Token{TokenType::LPAREN, input_.substr(start_pos, 1), start_pos};
      case ')':
        advance_();
        return Token{TokenType::RPAREN, input_.substr(start_pos, 1), start_pos};
      case '{':
        advance_();
        return Token{TokenType::LBRACE, input_.substr(start_pos, 1), start_pos};
      case '}':
        advance_();
        return Token{TokenType::RBRACE, input_.substr(start_pos, 1), start_pos};
      case '<':
        advance_();
        return Token{TokenType::LANGLE, input_.substr(start_pos, 1), start_pos};
      case '>':
        advance_();
        return Token{TokenType::RANGLE, input_.substr(start_pos, 1), start_pos};
      case '/':
        advance_();
        return Token{TokenType::SLASH, input_.substr(start_pos, 1), start_pos};
      case '|':
        advance_();
        return Token{TokenType::PIPE, input_.substr(start_pos, 1), start_pos};
      case '#':
        advance_();
        return Token{TokenType::HASH, input_.substr(start_pos, 1), start_pos};
      case ':':
        advance_();
        return Token{TokenType::COLON, input_.substr(start_pos, 1), start_pos};
      case ',':
        advance_();
        return Token{TokenType::COMMA, input_.substr(start_pos, 1), start_pos};
      case '^':
        advance_();
        return Token{TokenType::CARET, input_.substr(start_pos, 1), start_pos};
      case '?':
        advance_();
        return Token{TokenType::QUESTION, input_.substr(start_pos, 1), start_pos};
      case '@':
        advance_();
        return Token{TokenType::AT, input_.substr(start_pos, 1), start_pos};
      default:
        break;
    }

    // Check for +/- which could be part of a number or standalone
    // If followed by a digit or '.', treat as number; otherwise standalone
    if (c == '+' || c == '-')
    {
      char next_char = peek_(1);
      if (isDigit_(next_char) || next_char == '.')
      {
        return scanNumber_();
      }
      else
      {
        // Return as standalone PLUS or MINUS token
        advance_();
        return Token{c == '+' ? TokenType::PLUS : TokenType::MINUS,
                     input_.substr(start_pos, 1), start_pos};
      }
    }

    // Number (starts with digit or decimal point)
    if (isDigit_(c) || (c == '.' && isDigit_(peek_(1))))
    {
      return scanNumber_();
    }

    // Identifier (letter sequence)
    if (isLetter_(c))
    {
      return scanIdentifier_();
    }

    // Unrecognized character - return it as a single-character identifier
    // This allows the parser to produce a better error message
    advance_();
    return Token{TokenType::IDENTIFIER, input_.substr(start_pos, 1), start_pos};
  }

  ProFormaTokenizer::Token ProFormaTokenizer::scanNumber_()
  {
    size_t start_pos = pos_;

    // Handle optional leading sign
    if (current_() == '+' || current_() == '-')
    {
      advance_();
    }

    // Check if we have a leading decimal point (e.g., ".5" or "+.5")
    if (current_() == '.' && isDigit_(peek_(1)))
    {
      advance_(); // consume the '.'

      // Consume fractional part
      while (!isAtEnd_() && isDigit_(current_()))
      {
        advance_();
      }

      size_t length = pos_ - start_pos;
      return Token{TokenType::NUMBER, input_.substr(start_pos, length), start_pos};
    }

    // Consume integer part
    while (!isAtEnd_() && isDigit_(current_()))
    {
      advance_();
    }

    // Check for decimal point followed by digits
    if (!isAtEnd_() && current_() == '.' && isDigit_(peek_(1)))
    {
      advance_(); // consume the '.'

      // Consume fractional part
      while (!isAtEnd_() && isDigit_(current_()))
      {
        advance_();
      }
    }

    size_t length = pos_ - start_pos;
    return Token{TokenType::NUMBER, input_.substr(start_pos, length), start_pos};
  }

  ProFormaTokenizer::Token ProFormaTokenizer::scanIdentifier_()
  {
    size_t start_pos = pos_;

    // Consume all letters
    while (!isAtEnd_() && isLetter_(current_()))
    {
      advance_();
    }

    size_t length = pos_ - start_pos;
    return Token{TokenType::IDENTIFIER, input_.substr(start_pos, length), start_pos};
  }

  bool ProFormaTokenizer::isAtEnd_() const
  {
    return pos_ >= input_.size();
  }

  char ProFormaTokenizer::current_() const
  {
    if (isAtEnd_())
    {
      return '\0';
    }
    return input_[pos_];
  }

  char ProFormaTokenizer::peek_(size_t offset) const
  {
    size_t index = pos_ + offset;
    if (index >= input_.size())
    {
      return '\0';
    }
    return input_[index];
  }

  char ProFormaTokenizer::advance_()
  {
    if (isAtEnd_())
    {
      return '\0';
    }
    return input_[pos_++];
  }

  bool ProFormaTokenizer::isLetter_(char c)
  {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
  }

  bool ProFormaTokenizer::isDigit_(char c)
  {
    return c >= '0' && c <= '9';
  }

  bool ProFormaTokenizer::isWhitespace_(char c)
  {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
  }

} // namespace OpenMS
