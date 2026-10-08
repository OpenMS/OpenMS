// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS::Internal
{
  /**
    @brief Xerces-free, lightweight view over the attribute list of an XML element.

    Wraps a live @c xercesc::Attributes (owned by the SAX parser and valid only
    for the duration of the @c startElement callback). It exposes only the small
    surface the OpenMS handlers actually use (value-by-name, index-by-name,
    value-by-index), so that Xerces types never appear in a public header.

    The underlying attribute list is held as an opaque pointer, so this header
    pulls in no Xerces declarations — consumers of libOpenMS need Xerces neither
    on their include path nor at link time. The accessor implementations
    (in @c XMLAttributes.cpp) recover the concrete Xerces type internally.

    The wrapper is a non-owning view: it stores a pointer to the parser-owned
    attribute list and must not outlive the callback it was created for.
  */
  class OPENMS_DLLAPI XMLAttributes
  {
  public:
    XMLAttributes(const XMLAttributes&) = default;
    XMLAttributes& operator=(const XMLAttributes&) = default;

    /**
      @brief Wrap a live attribute list (non-owning).

      Takes an opaque handle so this header need not name any Xerces type. The
      only caller is the internal SAX2HandlerAdapter, which passes the address of
      the live @c xercesc::Attributes (recovered internally in @c XMLAttributes.cpp).
    */
    explicit XMLAttributes(const void* attributes) noexcept :
      attributes_(attributes)
    {
    }

    /// Value of the attribute named @p qname (UTF-16 code units), or @c nullptr if absent.
    const char16_t* value(const char16_t* qname) const;

    /// Value of the attribute named @p qname (narrow name, transcoded internally), or @c nullptr if absent.
    const char16_t* value(const char* qname) const;

    /// Index of the attribute named @p qname, or a negative value if absent.
    SignedSize index(const char16_t* qname) const;

    /// Value of the attribute at position @p i (UTF-16 code units).
    const char16_t* valueByIndex(Size i) const;

    /// Opaque handle to the underlying attribute list (libOpenMS-internal use).
    const void* handle() const noexcept { return attributes_; }

  private:
    const void* attributes_;
  };

} // namespace OpenMS::Internal
