// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <algorithm> // for "min"
#include <string>
#include <cstring>
#include <vector>

class QString;

namespace OpenMS
{

  /**
    *  Minimal replacement for boost::string_ref or std::experimental::string_view until we increase our min boost version
    *  @brief StringView provides a non-owning view on an existing string.
    */ 
  class OPENMS_DLLAPI StringView
  {
    public:

    // create view on string
    StringView() = default;

    // construct from other view
    StringView(const StringView&) = default;

    // copy assignment
    StringView& operator=(const StringView&) = default;

    // create view on string
    StringView(const std::string& s) : begin_(s.data()), size_(s.size())
    {
    }

    /// less operator
    bool operator<(const StringView other) const
    {
      if (size_ < other.size_) return true;

      if (size_ > other.size_) return false;

      // same size
      // same sequence, if both Views point to the same start
      if (begin_ == other.begin_) return false;

      return strncmp(begin_, other.begin_, size_) < 0;
    }

    bool operator==(const StringView other) const
    {
      if (size_ != other.size_) return false;

      //same size
      // same sequence, if both Views point to the same start
      if (begin_ == other.begin_) return true;

      return strncmp(begin_, other.begin_, size_) == 0;
    }

    /// create view that references a substring of the original string
    inline StringView substr(Size start, Size length) const
    {
      if (!size_) return *this;

      StringView sv(*this);
      sv.begin_ = begin_ + start;
      sv.size_ = std::min(length, sv.size_ - start);
      return sv;
    }
    
    /// size of view
    inline Size size() const
    {
      return size_;
    }

    /// create String object from view
    inline String getString() const
    {
      if (!size_) return String();
      return String(begin_, begin_ + size_);
    }

    private:
      const char* begin_;
      Size size_;
  };
	
} // namespace OpenMS

