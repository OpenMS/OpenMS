// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <boost/shared_ptr.hpp>

#include <compare>
#include <vector>

namespace OpenMS
{

  // missing in libc++
  template <typename InputIt1, typename InputIt2, typename Comp = std::compare_three_way>
  auto lexicographical_compare_three_way(InputIt1 first1, InputIt1 last1,
                                                 InputIt2 first2, InputIt2 last2,
                                                 Comp comp = {}) {
    for (; first1 != last1 && first2 != last2; ++first1, ++first2) {
        auto cmp = comp(*first1, *first2);
        if (cmp != 0)
            return cmp;
    }
    return (first1 == last1) ?
           ((first2 == last2) ? std::strong_ordering::equal : std::strong_ordering::less) :
           std::strong_ordering::greater;
  }

  // For ADL version of <=> for STL containers
  // This is needed to ensure that the std::vector<T> <=> operator is used
  template <typename T>
  auto operator<=>(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    return std::lexicographical_compare_three_way(
        lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::compare_three_way{}
    );
  }

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


