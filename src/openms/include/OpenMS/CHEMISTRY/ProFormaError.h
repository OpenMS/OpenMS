// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

  /** @ingroup Chemistry

      @brief Error codes for programmatic handling of ProForma parse errors.

      These error codes provide machine-readable categorization of parsing failures,
      enabling downstream code to handle specific error types appropriately.
  */
  enum class ProFormaErrorCode
  {
    UNEXPECTED_CHARACTER,       ///< Unexpected character encountered during parsing
    UNCLOSED_BRACKET,           ///< Opening bracket without matching close bracket
    UNMATCHED_BRACKET,          ///< Closing bracket without matching open bracket
    INVALID_CV_PREFIX,          ///< Invalid controlled vocabulary prefix (e.g., not UNIMOD, MOD, etc.)
    INVALID_CV_ACCESSION,       ///< Invalid CV accession number format
    INVALID_AMINO_ACID,         ///< Invalid amino acid one-letter code
    INVALID_MASS_VALUE,         ///< Invalid mass value format or value
    INVALID_FORMULA,            ///< Invalid chemical formula
    UNKNOWN_MONOSACCHARIDE,     ///< Unknown monosaccharide abbreviation
    DANGLING_CROSSLINK_LABEL,   ///< Crosslink label without a matching partner
    EMPTY_SEQUENCE,             ///< Empty sequence provided
    INVALID_CHARGE,             ///< Invalid charge state specification
    INVALID_OCCURRENCE_SPECIFIER, ///< Invalid occurrence specifier (e.g., ^2)
    UNEXPECTED_END_OF_INPUT,    ///< Unexpected end of input string
    INTERNAL_ERROR              ///< Internal parser error (should not occur)
  };

  /**
      @brief Convert error code to human-readable string.

      @param[in] code The error code to convert.
      @return A human-readable string describing the error code.
  */
  OPENMS_DLLAPI const char* proFormaErrorCodeToString(ProFormaErrorCode code);

  /**
      @brief Structured parse error with context for ProForma parsing.

      This exception extends Exception::ParseError to provide more detailed
      information about ProForma parsing failures, including:
      - Machine-readable error code for programmatic handling
      - Position of the error in the input string
      - Context snippets showing surrounding text
      - Expected vs. found information for clearer diagnostics

      Usage example:
      @code
      try {
        parseProForma("PEPTIDE[Oxidation");
      } catch (const ProFormaParseError& e) {
        std::cerr << e.getFormattedMessage() << std::endl;
        // Output:
        // ProForma parse error at position 7: Unclosed bracket
        // Context: PEPTIDE[Oxidation
        //                 ^
        // Expected: ']'
        // Found: end of input
      }
      @endcode

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProFormaParseError :
    public Exception::ParseError
  {
  public:
    /**
        @brief Constructs a ProFormaParseError with full context information.

        @param[in] file Source file where the exception was thrown (use __FILE__).
        @param[in] line Source line where the exception was thrown (use __LINE__).
        @param[in] function Function name where the exception was thrown (use OPENMS_PRETTY_FUNCTION).
        @param[in] error_code Machine-readable error code categorizing the error.
        @param[in] error_position Byte position (0-indexed) in the input where the error occurred.
                                  Will be clamped to input.size() if out of bounds.
        @param[in] input The complete input string that was being parsed.
        @param[in] message Human-readable description of the error.
    */
    ProFormaParseError(
      const char* file,
      int line,
      const char* function,
      ProFormaErrorCode error_code,
      size_t error_position,
      const String& input,
      const String& message
    ) noexcept;

    /// Get the error code for programmatic handling
    ProFormaErrorCode getErrorCode() const noexcept
    {
      return code_;
    }

    /// Get byte position of error (0-indexed)
    size_t getPosition() const noexcept
    {
      return position_;
    }

    /// Get context before error (~20 chars)
    const String& getContextBefore() const noexcept
    {
      return context_before_;
    }

    /// Get context after error (~20 chars)
    const String& getContextAfter() const noexcept
    {
      return context_after_;
    }

    /// Get what was expected (if set)
    const String& getExpected() const noexcept
    {
      return expected_;
    }

    /// Get what was found (if set)
    const String& getFound() const noexcept
    {
      return found_;
    }

    /**
        @brief Get formatted error message with context visualization.

        Returns a multi-line string formatted for human reading, showing:
        - Error description with position
        - Context snippet with caret pointing to error location
        - Expected vs. found information (if available)

        @return Formatted error message.
    */
    String getFormattedMessage() const;

    /**
        @brief Set expected/found information for more detailed error messages.

        @param[in] expected Description of what was expected at the error position.
        @param[in] found Description of what was actually found.
    */
    void setExpectedFound(const String& expected, const String& found);

  private:
    ProFormaErrorCode code_;      ///< Machine-readable error code
    size_t position_;             ///< Position in input where error occurred
    String context_before_;       ///< Text before the error position
    String context_after_;        ///< Text after the error position
    String expected_;             ///< What was expected (optional)
    String found_;                ///< What was found (optional)

    /**
        @brief Extract context snippets from input around the error position.

        Extracts approximately 20 characters before and after the error position
        for display in error messages. The position is assumed to be already
        clamped to input.size().

        @param[in] input The complete input string.
        @param[in] pos The error position in the input (clamped to input.size()).
    */
    void extractContext_(const String& input, size_t pos);
  };

} // namespace OpenMS
