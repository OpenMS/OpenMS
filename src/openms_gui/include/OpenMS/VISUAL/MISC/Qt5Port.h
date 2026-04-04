// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <QSet>
#include <QString>
#include <QStringList>

#include <vector>

namespace OpenMS
{

/// Drop-in for QT5's QStringList::toSet
template <typename T, template<typename> typename C>
QSet<T> toQSet(const C<T> &container)
{
  return QSet<T>(container.begin(), container.end());
}

/// Convert QStringList to StringList (replaces StringListUtils::fromQStringList)
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

/// Convert OpenMS::String to QString (replaces String::toQString())
inline QString toQString(const String& s) { return QString::fromStdString(s); }

/// Construct OpenMS::String from QString (replaces String(const QString&))
inline String fromQString(const QString& s) { return String(s.toStdString()); }

} // namespace OpenMS
