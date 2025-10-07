// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

