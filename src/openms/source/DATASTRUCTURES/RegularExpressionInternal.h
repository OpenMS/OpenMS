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
      static const boost::regex empty;
      return expression.impl_ ? expression.impl_->expression : empty;
    }
  };

  struct RegularExpressionMatch : boost::smatch
  {
  };
} // namespace Internal
} // namespace OpenMS
