// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <iterator>
#include <ostream>
#include <vector>

// This header collects io relevant parts of ListUtils. Separating the from the
// rest avoids inclusion of ostream headers in a lot of classes.

namespace OpenMS
{
  /**
    @brief Output stream operator for std::vectors.

    @param os The target stream.
    @param v The vector to write to stream.
  */
  template <typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
  {
    os << "[";

    if (!v.empty())
    {
      for (auto it = v.begin(); it < v.end() - 1; ++it)
      { // convert to String manually, since this is much faster than ostream build-in conversion; 
        // If T is a String, the compiler will (hopefully) elide the copy
        os << String(*it) << ", ";
      }
      os << String(v.back());
    }

    os << "]";
    return os;
  }

  template<typename T>
  struct VecLowPrecision
  {
    const std::vector<T>& value;
    VecLowPrecision(const std::vector<T>& v) : value(v) {}
  };
  /// modified version of the stream operator (works for vectors of float, double, long double) which prints only
  /// three fractional digits; usage 'os << VecLowPrecision(my_vec);'
  template <typename T>
  inline std::ostream& operator<<(std::ostream& os, const VecLowPrecision<T>& val)
  {
    os << "[";
    const auto& v = val.value;
    if (!v.empty())
    {
      for (auto it = v.begin(); it < v.end() - 1; ++it)
      { // convert to String manually, since this is much faster than ostreams build-in conversion; 
        os << String(*it, false) << ", ";
      }
      os << String(v.back(), false);
    }

    os << "]";
    return os;
  }

  /// Operator for appending entries with less code
  template <typename TString>
  inline std::vector<String>& operator<<(std::vector<String>& sl, const TString& string)
  {
    sl.push_back(string);
    return sl;
  }

}

