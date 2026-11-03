// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "RegularExpressionInternal.h"

namespace OpenMS
{
RegularExpression::RegularExpression() = default;
RegularExpression::RegularExpression(const std::string& pattern): impl_(std::make_shared<Impl>(pattern))
{
}
void RegularExpression::assign(const std::string& pattern)
{ impl_ = std::make_shared<Impl>(pattern); }
std::string RegularExpression::str() const
{ return impl_ ? impl_->expression.str() : std::string(); }
bool RegularExpression::empty() const
{ return ! impl_ || impl_->expression.empty(); }
bool RegularExpression::search(const std::string& text, std::string* first_match) const
{
  const boost::regex& expression = Internal::RegularExpressionAccess::get(*this);
  // Skip building a match_results when the caller does not want the matched text: that is the
  // overload our callers used before, and it also lets Boost add regex_constants::match_any.
  // Measured on the patterns we actually search with, this is worth ~10% on short inputs and
  // nothing at all on the anchored decoy-affix pattern -- it is not a hot-path fix, just the
  // cheaper way to answer a yes/no question.
  if (first_match == nullptr) { return boost::regex_search(text, expression); }

  boost::smatch match;
  const bool found = boost::regex_search(text, match, expression);
  if (found) { first_match->assign(match[0].first, match[0].second); }
  else { first_match->clear(); }
  return found;
}
} // namespace OpenMS
