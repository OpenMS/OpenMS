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
#include <string_view>
#include <utility>
#include <vector>

namespace OpenMS::Internal
{
  /**
    @brief Non-owning attribute view for an XML start-element callback.

    The view refers either to UTF-16 name/value pairs or to a parser-owned
    @ref NativeStore of UTF-8 slices. The latter converts a value to UTF-16
    only when it is requested through value() or valueByIndex(); valueUTF8()
    copies it without conversion. This view and the pointers it returns must
    not outlive the callback that receives them.
  */
  class OPENMS_DLLAPI XMLAttributes
  {
  public:
    using Values = std::vector<std::pair<std::u16string, std::u16string>>;

    /**
      @brief Reusable UTF-8 attribute storage of an XML parser adapter.

      Names and values are borrowed from the parser for one callback. The
      storage keeps its capacity across clear(), so a parse does not allocate
      per element once the largest element has been seen.
    */
    class OPENMS_DLLAPI NativeStore
    {
    public:
      /// Forgets the attributes of the previous element (keeps capacity).
      void clear() noexcept { count_ = 0; }
      /// Appends an attribute; @p prefix may be empty. The qualified name is "prefix:local".
      void add(std::string_view prefix, std::string_view local, std::string_view value);
      /// Number of attributes
      Size size() const noexcept { return count_; }

    private:
      friend class XMLAttributes;
      struct Entry
      {
        std::string_view local;         ///< UTF-8 name if unprefixed
        std::string_view value;         ///< UTF-8 value
        std::string joined;             ///< UTF-8 "prefix:local" if prefixed
        bool prefixed = false;
        mutable std::u16string value16; ///< UTF-16 value, valid if converted
        mutable bool converted = false;

        /// UTF-8 qualified name (a view into joined is only taken here, as entries move when the store grows)
        std::string_view name() const noexcept { return prefixed ? std::string_view(joined) : local; }
      };
      std::vector<Entry> entries_;
      Size count_ = 0;
    };

    explicit XMLAttributes(const Values& attributes) noexcept : attributes_(&attributes) {}
    explicit XMLAttributes(const NativeStore& attributes) noexcept : native_(&attributes) {}

    XMLAttributes(Values&&) = delete;
    XMLAttributes(const Values&&) = delete;
    XMLAttributes(NativeStore&&) = delete;
    XMLAttributes(const NativeStore&&) = delete;

    /// Number of attributes
    Size size() const noexcept;

    /// UTF-16 value by name, or nullptr if absent.
    const char16_t* value(const char16_t* qname) const;

    /// UTF-16 value by UTF-8 name, or nullptr if absent.
    const char16_t* value(const char* qname) const;

    /// Index by UTF-16 name, or a negative value if absent.
    SignedSize index(const char16_t* qname) const;

    /// UTF-16 value by index, or nullptr if out of range.
    const char16_t* valueByIndex(Size i) const;

    /// Assigns the UTF-8 value of the attribute named @p qname to @p value. Returns false (leaving @p value untouched) if absent.
    bool valueUTF8(const char16_t* qname, std::string& value) const;

    /// Assigns the UTF-8 value of the attribute with UTF-8 name @p qname to @p value. Returns false (leaving @p value untouched) if absent.
    bool valueUTF8(const char* qname, std::string& value) const;

  private:
    /// Index by UTF-8 name, or a negative value if absent.
    SignedSize indexUTF8_(const char* qname) const;
    /// UTF-8 value at index @p i (which must be valid)
    void valueUTF8ByIndex_(Size i, std::string& value) const;

    const Values* attributes_ = nullptr;
    const NativeStore* native_ = nullptr;
  };
}
