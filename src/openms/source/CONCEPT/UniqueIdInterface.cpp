// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/UniqueIdInterface.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>


namespace OpenMS
{
  Size UniqueIdInterface::setUniqueId()
  {
    unique_id_ = UniqueIdGenerator::getUniqueId();
    return 1;
  }

  Size UniqueIdInterface::ensureUniqueId()
  {
    if (!hasValidUniqueId())
    {
      unique_id_ = UniqueIdGenerator::getUniqueId();
      return 1;
    }
    else
      return 0;
  }

  void UniqueIdInterface::setUniqueId(const String & rhs)
  {
    clearUniqueId();

    String::size_type last_underscore = rhs.rfind('_');
    String s = rhs.substr(last_underscore + 1);

    for (String::const_iterator s_i = s.begin(); s_i < s.end(); ++s_i)
    {
      int i = (*s_i - '0');
      if (i < 0 || i > 9)
      {
        clearUniqueId();
        return;
      }
      unique_id_ = 10 * unique_id_ + i;
    }

  }

}
