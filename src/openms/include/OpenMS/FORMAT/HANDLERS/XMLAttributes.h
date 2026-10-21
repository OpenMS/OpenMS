// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS::Internal
{
  /**
    @brief Non-owning UTF-16 attribute view for an XML start-element callback.

    The parser adapter owns the attribute vector. This view and the pointers it
    returns must not outlive the callback that receives them.
  */
  class OPENMS_DLLAPI XMLAttributes
  {
  public:
    using Values = std::vector<std::pair<std::u16string, std::u16string>>;

    explicit XMLAttributes(const Values& attributes) noexcept : attributes_(&attributes) {}

    XMLAttributes(Values&&) = delete;
    XMLAttributes(const Values&&) = delete;

    /// UTF-16 value by name, or nullptr if absent.
    const char16_t* value(const char16_t* qname) const;

    /// UTF-16 value by UTF-8 name, or nullptr if absent.
    const char16_t* value(const char* qname) const;

    /// Index by UTF-16 name, or a negative value if absent.
    SignedSize index(const char16_t* qname) const;

    /// UTF-16 value by index, or nullptr if out of range.
    const char16_t* valueByIndex(Size i) const;

  private:
    const Values* attributes_;
  };
}
