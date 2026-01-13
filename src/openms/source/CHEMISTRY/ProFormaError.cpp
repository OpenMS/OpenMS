// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaError.h>
#include <OpenMS/CONCEPT/GlobalExceptionHandler.h>

#include <sstream>
#include <algorithm>

namespace OpenMS
{

  const char* proFormaErrorCodeToString(ProFormaErrorCode code)
  {
    switch (code)
    {
      case ProFormaErrorCode::UNEXPECTED_CHARACTER:
        return "Unexpected character";
      case ProFormaErrorCode::UNCLOSED_BRACKET:
        return "Unclosed bracket";
      case ProFormaErrorCode::UNMATCHED_BRACKET:
        return "Unmatched closing bracket";
      case ProFormaErrorCode::INVALID_CV_PREFIX:
        return "Invalid controlled vocabulary prefix";
      case ProFormaErrorCode::INVALID_CV_ACCESSION:
        return "Invalid CV accession number";
      case ProFormaErrorCode::INVALID_AMINO_ACID:
        return "Invalid amino acid";
      case ProFormaErrorCode::INVALID_MASS_VALUE:
        return "Invalid mass value";
      case ProFormaErrorCode::INVALID_FORMULA:
        return "Invalid chemical formula";
      case ProFormaErrorCode::UNKNOWN_MONOSACCHARIDE:
        return "Unknown monosaccharide";
      case ProFormaErrorCode::DANGLING_CROSSLINK_LABEL:
        return "Dangling crosslink label";
      case ProFormaErrorCode::EMPTY_SEQUENCE:
        return "Empty sequence";
      case ProFormaErrorCode::INVALID_CHARGE:
        return "Invalid charge state";
      case ProFormaErrorCode::INVALID_OCCURRENCE_SPECIFIER:
        return "Invalid occurrence specifier";
      case ProFormaErrorCode::UNEXPECTED_END_OF_INPUT:
        return "Unexpected end of input";
      case ProFormaErrorCode::INTERNAL_ERROR:
        return "Internal parser error";
      default:
        return "Unknown error";
    }
  }

  ProFormaParseError::ProFormaParseError(
    const char* file,
    int line,
    const char* function,
    ProFormaErrorCode error_code,
    size_t error_position,
    const String& input,
    const String& message
  ) noexcept :
    ParseError(file, line, function, input, message),
    code_(error_code),
    position_(error_position)
  {
    extractContext_(input, error_position);
    GlobalExceptionHandler::getInstance().setMessage(what());
  }

  void ProFormaParseError::extractContext_(const String& input, size_t pos)
  {
    const size_t context_length = 20;

    // Extract context before the error position
    if (pos > 0)
    {
      size_t start = (pos > context_length) ? pos - context_length : 0;
      context_before_ = input.substr(start, pos - start);
    }
    else
    {
      context_before_ = "";
    }

    // Extract context after the error position
    if (pos < input.size())
    {
      size_t length = std::min(context_length, input.size() - pos);
      context_after_ = input.substr(pos, length);
    }
    else
    {
      context_after_ = "";
    }
  }

  String ProFormaParseError::getFormattedMessage() const
  {
    std::ostringstream oss;

    // First line: error description with position
    oss << "ProForma parse error at position " << position_ << ": "
        << proFormaErrorCodeToString(code_);

    // Second line: context with the error position highlighted
    oss << "\nContext: ";

    // Add ellipsis if there's more text before the context
    if (position_ > context_before_.size())
    {
      oss << "...";
    }

    oss << context_before_;

    // Mark the error position
    if (!context_after_.empty())
    {
      oss << ">>>" << context_after_.substr(0, 1) << "<<<";
      if (context_after_.size() > 1)
      {
        oss << context_after_.substr(1);
      }
    }
    else
    {
      oss << ">>><END OF INPUT><<<";
    }

    // Add ellipsis if there's more text after the context
    // (We compare against context_length which was used to extract the context)
    if (context_after_.size() >= 20)
    {
      oss << "...";
    }

    // Expected/Found information if available
    if (!expected_.empty() || !found_.empty())
    {
      if (!expected_.empty())
      {
        oss << "\nExpected: " << expected_;
      }
      if (!found_.empty())
      {
        oss << "\nFound: " << found_;
      }
    }

    return String(oss.str());
  }

  void ProFormaParseError::setExpectedFound(const String& expected, const String& found)
  {
    expected_ = expected;
    found_ = found;
  }

} // namespace OpenMS
