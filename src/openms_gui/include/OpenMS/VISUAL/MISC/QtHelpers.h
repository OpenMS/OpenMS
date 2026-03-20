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
#include <QStringList>

#include <vector>

namespace OpenMS
{
  /// Convert OpenMS::String (or std::string) to QString
  /// Free function replacing the former String::toQString() member
  inline QString toQString(const std::string& s)
  {
    return QString::fromStdString(s);
  }

  /// Convert QStringList to std::vector<String> (replaces StringListUtils::fromQStringList)
  inline std::vector<String> fromQStringList(const QStringList& rhs)
  {
    std::vector<String> sl;
    sl.reserve(rhs.size());
    for (const auto& item : rhs)
    {
      sl.push_back(item.toStdString());
    }
    return sl;
  }
} // namespace OpenMS
