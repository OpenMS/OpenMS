// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  namespace Helpers 
  {

    /**
        @brief Helper function to add constness to a vector of shared pointers
    */
    template <class T>
    const std::vector<boost::shared_ptr<const T> >&
    constifyPointerVector(const std::vector<boost::shared_ptr<T> >& vec) 
    {
      return reinterpret_cast<const std::vector<boost::shared_ptr<const T> >&>(vec);
    }


    /**
      * @brief Helper comparing two pointers for equality (taking NULL into account)
    */
    template <class PtrType>
    inline bool cmpPtrSafe(const PtrType& a, const PtrType& b)
    {
       // We are not interested whether the pointers are equal but whether the
       // contents are equal
      if (a == nullptr && b == nullptr)
      {
        return true;
      }
      else if (a == nullptr || b == nullptr)
      {
        return false; // one is null the other is not
      }
      else
      {
        // compare the internal object
        return (*a == *b);
      }
    }

    /**
      * @brief Helper function to compare two pointer-containers for equality of all elements
    */
    template <class ContainerType>
    inline bool cmpPtrContainer(const ContainerType& a, const ContainerType& b)
    {
      if (a.size() != b.size()) return false;

      // check that all elements of a and b are equal using safe comparison
      // (taking NULL into account)
      for (typename ContainerType::size_type i = 0; i < a.size(); i++)
      {
        if (!cmpPtrSafe(a[i], b[i]))
        {
          return false;
        }
      }
      return true;
    }

  }
}


