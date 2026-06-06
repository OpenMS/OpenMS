// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// All string utilities live here: toStr(), appendToStr(), trim(), split(),
// global operator+(std::string, numeric), QuotingMethod enum, etc.
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <string>

namespace OpenMS
{
  /// OpenMS::String is now a plain std::string.
  ///
  /// Migration guide for former OpenMS::String member calls:
  ///   s.trim()                  →  StringUtils::trim(s)
  ///   s.toLower()               →  StringUtils::toLower(s)
  ///   s.toUpper()               →  StringUtils::toUpper(s)
  ///   s.split(delim, parts)     →  StringUtils::split(s, delim, parts)
  ///   s.toInt() / toDouble()    →  StringUtils::toInt32(s) / toDouble(s)
  ///   String(42) / String(3.14) →  StringUtils::toStr(42) / toStr(3.14)
  ///   String::random(n)         →  StringUtils::random(n)
  ///   String::number(d, n)      →  StringUtils::number(d, n)
  ///   String::EMPTY             →  ""  or  std::string{}
  ///   String::ESCAPE            →  OpenMS::QuotingMethod::ESCAPE
  ///   String::NONE              →  OpenMS::QuotingMethod::NONE
  ///   String::DOUBLE            →  OpenMS::QuotingMethod::DOUBLE
}
