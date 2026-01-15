// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


namespace OpenMS
{

/// Drop-in for QT5's QStringList::toSet   
template <typename T, template<typename> typename C>
QSet<T> toQSet(const C<T> &container)
{
  return QSet<T>(container.begin(), container.end());
}

} // namespace OpenMS