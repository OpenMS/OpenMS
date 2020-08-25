// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
      { // convert to String manually, since this is much faster than ostreams build-in conversion; 
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

