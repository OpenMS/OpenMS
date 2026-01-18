// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <optional>
#include <string>
#include <string_view>

namespace OpenMS
{

  /**
    @brief Tokenizer for ProForma v2 peptidoform notation

    This class provides lexical analysis (tokenization) for ProForma strings.
    It produces tokens suitable for parsing the ProForma grammar, supporting
    zero-copy operation via std::string_view for performance.

    The tokenizer handles:
    - Single-character tokens: [ ] ( ) { } < > + - / | # : , ^ ? @
    - Numbers: integers and decimals (e.g., "123", "15.9949", "+2", "-1")
    - Identifiers: letter sequences (e.g., "UNIMOD", "Oxidation", "PEPTIDE")
    - End of input marker

    @note The tokenizer does not skip whitespace - ProForma strings should not contain
          whitespace according to the specification.

    @note The input string must remain valid for the lifetime of the tokenizer
          since tokens reference slices of the input.

    Usage example:
    @code
    std::string input = "EM[UNIMOD:35]K";
    ProFormaTokenizer tokenizer(input);

    while (tokenizer.hasMore())
    {
        ProFormaTokenizer::Token token = tokenizer.next();
        std::cout << "Token: " << token.text << " at position " << token.position << std::endl;
    }
    @endcode

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProFormaTokenizer
  {
  public:

    /// Token types produced by the tokenizer
    enum class TokenType
    {
      LBRACKET,    ///< Left square bracket: [
      RBRACKET,    ///< Right square bracket: ]
      LPAREN,      ///< Left parenthesis: (
      RPAREN,      ///< Right parenthesis: )
      LBRACE,      ///< Left curly brace: {
      RBRACE,      ///< Right curly brace: }
      LANGLE,      ///< Left angle bracket: <
      RANGLE,      ///< Right angle bracket: >
      PLUS,        ///< Plus sign: +
      MINUS,       ///< Minus sign: -
      SLASH,       ///< Forward slash: /
      PIPE,        ///< Vertical bar (pipe): |
      HASH,        ///< Hash/pound sign: #
      COLON,       ///< Colon: :
      COMMA,       ///< Comma: ,
      CARET,       ///< Caret: ^
      QUESTION,    ///< Question mark: ?
      AT,          ///< At sign: @
      NUMBER,      ///< Numeric literal (integer or decimal, possibly with leading +/-)
      IDENTIFIER,  ///< Letter sequence (A-Za-z)
      END          ///< End of input
    };

    /**
      @brief A single token from the input stream

      Tokens store their type, a view of the original text, and the byte
      position in the input for error reporting.
    */
    struct Token
    {
      TokenType type;         ///< The type of this token
      std::string_view text;  ///< View into the original input (zero-copy)
      size_t position;        ///< Byte offset in the original input (0-indexed)

      /// Check if this is an end-of-input token
      bool isEnd() const { return type == TokenType::END; }

      /// Check if this token is a specific single-character type
      bool is(TokenType t) const { return type == t; }
    };

    /**
      @brief Construct a tokenizer for the given input string

      @param input The ProForma string to tokenize. Must remain valid for
                   the lifetime of this tokenizer.
      @param start_pos Optional starting position (default 0). Used for efficient
                       lookahead without re-scanning from the beginning.
    */
    explicit ProFormaTokenizer(std::string_view input, size_t start_pos = 0);

    /// Default destructor
    ~ProFormaTokenizer() = default;

    /// Copy constructor
    ProFormaTokenizer(const ProFormaTokenizer&) = default;

    /// Move constructor
    ProFormaTokenizer(ProFormaTokenizer&&) = default;

    /// Copy assignment operator
    ProFormaTokenizer& operator=(const ProFormaTokenizer&) = default;

    /// Move assignment operator
    ProFormaTokenizer& operator=(ProFormaTokenizer&&) = default;

    /**
      @brief Consume and return the next token

      Advances the tokenizer position past the returned token.

      @return The next token from the input stream
    */
    Token next();

    /**
      @brief Look at the next token without consuming it

      Multiple calls to peek() without intervening next() calls return the same token.

      @return The next token that would be returned by next()
    */
    Token peek();

    /**
      @brief Check if more tokens are available

      @return true if there are more tokens (i.e., next() would not return END)
    */
    bool hasMore() const;

    /**
      @brief Get the current position in the input

      @return The byte offset of the next character to be scanned
    */
    size_t position() const;

    /**
      @brief Get a context string around a position for error messages

      Returns a substring of the input centered around the given position,
      useful for providing context in error messages.

      @param pos The position to center the context around
      @param before Maximum number of characters to include before pos
      @param after Maximum number of characters to include after pos
      @return A view of the context substring
    */
    std::string_view getContext(size_t pos, size_t before = 20, size_t after = 20) const;

    /**
      @brief Get a human-readable name for a token type

      @param type The token type
      @return A string describing the token type
    */
    static const char* tokenTypeName(TokenType type);

  private:

    /// Scan and return the next token from the current position
    Token scanToken_();

    /// Scan a number token (integer, decimal, optionally signed)
    Token scanNumber_();

    /// Scan an identifier token (letter sequence)
    Token scanIdentifier_();

    /// Check if we have reached the end of input
    bool isAtEnd_() const;

    /// Get the current character (or '\0' if at end)
    char current_() const;

    /// Get the character at offset from current position (or '\0' if out of bounds)
    char peek_(size_t offset) const;

    /// Advance to the next character and return the previous one
    char advance_();

    /// Check if a character is a letter (A-Za-z)
    static bool isLetter_(char c);

    /// Check if a character is a digit (0-9)
    static bool isDigit_(char c);

    /// The input string (must remain valid for tokenizer lifetime)
    std::string_view input_;

    /// Current position in the input
    size_t pos_ = 0;

    /// Cached peeked token (if any)
    std::optional<Token> peeked_;
  };

} // namespace OpenMS
