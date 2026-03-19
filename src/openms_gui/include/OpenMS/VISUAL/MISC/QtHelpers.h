// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QString>

namespace OpenMS
{
  /// Convert OpenMS::String to QString (free function replacing the former String::toQString() member)
  inline QString toQString(const String& s)
  {
    return QString::fromStdString(s);
  }

  /// Convert std::string to QString
  inline QString toQString(const std::string& s)
  {
    return QString::fromStdString(s);
  }
} // namespace OpenMS
