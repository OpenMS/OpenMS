// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <boost/mem_fn.hpp>
#include <QtCore/QStringList>

using namespace std;

namespace OpenMS
{

  StringList StringListUtils::fromQStringList(const QStringList& rhs)
  {
    StringList sl;
    sl.reserve(rhs.size());

    for (QStringList::const_iterator it = rhs.begin(); it != rhs.end(); ++it)
    {
      sl.push_back(it->toStdString());
    }

    return sl;
  }

  void StringListUtils::toUpper(StringList& sl)
  {
    std::for_each(sl.begin(), sl.end(), boost::mem_fn(&String::toUpper));
  }

  void StringListUtils::toLower(StringList& sl)
  {
    std::for_each(sl.begin(), sl.end(), boost::mem_fn(&String::toLower));
  }

  StringListUtils::Iterator StringListUtils::searchPrefix(const Iterator& start, const Iterator& end, const String& text, bool trim)
  {
    return find_if(start, end, PrefixPredicate_(text, trim));
  }

  StringListUtils::ConstIterator StringListUtils::searchPrefix(const ConstIterator& start, const ConstIterator& end, const String& text, bool trim)
  {
    return find_if(start, end, PrefixPredicate_(text, trim));
  }

  StringListUtils::ConstIterator StringListUtils::searchPrefix(const StringList& container, const String& text, bool trim)
  {
    return searchPrefix(container.begin(), container.end(), text, trim);
  }

  StringListUtils::Iterator StringListUtils::searchPrefix(StringList& container, const String& text, bool trim)
  {
    return searchPrefix(container.begin(), container.end(), text, trim);
  }

  StringListUtils::Iterator StringListUtils::searchSuffix(const Iterator& start, const Iterator& end, const String& text, bool trim)
  {
    return find_if(start, end, SuffixPredicate_(text, trim));
  }

  StringListUtils::ConstIterator StringListUtils::searchSuffix(const ConstIterator& start, const ConstIterator& end, const String& text, bool trim)
  {
    return find_if(start, end, SuffixPredicate_(text, trim));
  }

  StringListUtils::ConstIterator StringListUtils::searchSuffix(const StringList& container, const String& text, bool trim)
  {
    return searchSuffix(container.begin(), container.end(), text, trim);
  }

  StringListUtils::Iterator StringListUtils::searchSuffix(StringList& container, const String& text, bool trim)
  {
    return searchSuffix(container.begin(), container.end(), text, trim);
  }

} // namespace OpenMS
