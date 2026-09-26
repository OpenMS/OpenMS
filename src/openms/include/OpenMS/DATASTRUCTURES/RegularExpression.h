// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <memory>
#include <string>

namespace OpenMS
{
namespace Internal
{
  struct RegularExpressionAccess;
  struct RegularExpressionMatch;
} // namespace Internal

/**
  @brief Compiled Perl-compatible regular expression, including named capture groups.

  Copies share immutable compiled state. assign() replaces only this expression.
  The regex engine and its headers are private to OpenMS.
*/
class OPENMS_DLLAPI RegularExpression
{
public:
  /// Default constructor: holds no pattern (empty() is true) and matches like the pattern ""
  RegularExpression();
  /**
    @brief Compiles @p pattern (Perl syntax)
    @throws std::runtime_error if @p pattern is not a valid regular expression
  */
  explicit RegularExpression(const std::string& pattern);
  /**
    @brief Replaces the expression by @p pattern; copies made earlier keep the old one
    @throws std::runtime_error if @p pattern is not a valid regular expression
  */
  void assign(const std::string& pattern);
  /// The pattern this expression was compiled from; "" if it holds none
  std::string str() const;
  /// True only if no pattern was set; a compiled "" is not empty
  bool empty() const;
  /// Search text, optionally returning the first complete match.
  bool search(const std::string& text, std::string* first_match = nullptr) const;

private:
  struct Impl;
  std::shared_ptr<const Impl> impl_;
  friend struct Internal::RegularExpressionAccess;
};
} // namespace OpenMS
