// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


namespace OpenMS
{
  namespace Helpers
  {
    /// return the index of an element in a container
    /// useful for matching 'names_of_...' arrays to their enum value
    template <class ContainerType>
    Size indexOf(const ContainerType& cont, const typename ContainerType::value_type& val)
    {
      auto it = std::find(cont.begin(), cont.end(), val);
      if (it == cont.end())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, val);
      }
      return std::distance(cont.begin(), it);
    }

  }
}


