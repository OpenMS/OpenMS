// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/RegularExpression.h>
#include <boost/regex.hpp>

namespace OpenMS
{
struct RegularExpression::Impl
{
  boost::regex expression;
  Impl() = default;
  explicit Impl(const std::string& pattern): expression(pattern)
  {
  }
};

namespace Internal
{
  struct RegularExpressionAccess
  {
    static const boost::regex& get(const RegularExpression& expression)
    {
      // A default-constructed boost::regex holds no implementation and asserts (release: UB) as
      // soon as it is matched against, so the fallback for a pattern-less RegularExpression has
      // to be a compiled pattern. The empty pattern makes RegularExpression() behave exactly
      // like RegularExpression("").
      static const boost::regex empty("");
      return expression.impl_ ? expression.impl_->expression : empty;
    }
  };

  struct RegularExpressionMatch : boost::smatch
  {
  };
} // namespace Internal
} // namespace OpenMS
