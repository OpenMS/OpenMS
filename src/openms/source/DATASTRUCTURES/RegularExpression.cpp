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
  boost::smatch match;
  const bool found = boost::regex_search(text, match, Internal::RegularExpressionAccess::get(*this));
  if (first_match) *first_match = found ? match.str() : std::string();
  return found;
}
} // namespace OpenMS
