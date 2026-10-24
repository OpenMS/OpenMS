// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>

#include <OpenMS/FORMAT/HANDLERS/StringManager.h>

#include <cstring>

namespace OpenMS::Internal
{
  namespace
  {
    bool equalsASCII(const char16_t* wide, const char* ascii)
    {
      while (*ascii && *wide == static_cast<unsigned char>(*ascii))
      {
        ++wide;
        ++ascii;
      }
      return *ascii == 0 && *wide == 0;
    }

    bool isASCII(const char* str)
    {
      for (; *str; ++str)
      {
        if (static_cast<unsigned char>(*str) > 127) return false;
      }
      return true;
    }

    /// Compares a null-terminated UTF-16 name with a UTF-8 name without converting either.
    bool equalsUTF8(const char16_t* wide, std::string_view utf8)
    {
      const auto* in = reinterpret_cast<const unsigned char*>(utf8.data());
      const auto* end = in + utf8.size();
      while (in != end)
      {
        unsigned code = *in++;
        if (code >= 0x80)
        {
          // The parser delivers valid UTF-8; treat anything else as a mismatch.
          const int continuation = code >= 0xf0 ? 3 : code >= 0xe0 ? 2 : code >= 0xc0 ? 1 : -1;
          if (continuation < 0 || end - in < continuation) return false;
          code &= 0x3f >> continuation;
          for (int i = 0; i < continuation; ++i) code = (code << 6) | (*in++ & 0x3f);
          if (code > 0xffff)
          {
            code -= 0x10000;
            if (*wide++ != static_cast<char16_t>(0xd800 + (code >> 10))) return false;
            code = 0xdc00 + (code & 0x3ff);
          }
        }
        if (*wide++ != static_cast<char16_t>(code)) return false;
      }
      return *wide == 0;
    }
  }

  void XMLAttributes::NativeStore::add(std::string_view prefix, std::string_view local, std::string_view value)
  {
    if (count_ == entries_.size()) entries_.emplace_back();
    Entry& entry = entries_[count_++];
    entry.prefixed = !prefix.empty();
    if (entry.prefixed) { entry.joined.assign(prefix).append(1, ':').append(local); }
    else { entry.local = local; }
    entry.value = value;
    entry.converted = false;
  }

  Size XMLAttributes::size() const noexcept
  {
    return native_ ? native_->count_ : attributes_->size();
  }

  const char16_t* XMLAttributes::value(const char16_t* qname) const
  {
    const SignedSize i = index(qname);
    return i < 0 ? nullptr : valueByIndex(static_cast<Size>(i));
  }

  const char16_t* XMLAttributes::value(const char* qname) const
  {
    const SignedSize i = indexUTF8_(qname);
    return i < 0 ? nullptr : valueByIndex(static_cast<Size>(i));
  }

  SignedSize XMLAttributes::index(const char16_t* qname) const
  {
    if (!qname) return -1;
    if (native_)
    {
      for (Size i = 0; i < native_->count_; ++i)
      {
        if (equalsUTF8(qname, native_->entries_[i].name())) return static_cast<SignedSize>(i);
      }
      return -1;
    }
    for (Size i = 0; i < attributes_->size(); ++i)
    {
      if ((*attributes_)[i].first == qname) return static_cast<SignedSize>(i);
    }
    return -1;
  }

  SignedSize XMLAttributes::indexUTF8_(const char* qname) const
  {
    if (!qname) return -1;
    if (!isASCII(qname))
    {
      // Reject invalid UTF-8 names (throws Exception::ConversionError).
      const auto name = StringManager::convert(qname);
      if (!native_) return index(name.c_str());
    }
    if (native_)
    {
      // Valid UTF-8 strings are equal exactly if their bytes are.
      const std::string_view name(qname);
      for (Size i = 0; i < native_->count_; ++i)
      {
        if (native_->entries_[i].name() == name) return static_cast<SignedSize>(i);
      }
      return -1;
    }
    // XML vocabulary names are usually ASCII; avoid allocating a transcoder.
    for (Size i = 0; i < attributes_->size(); ++i)
    {
      if (equalsASCII((*attributes_)[i].first.c_str(), qname)) return static_cast<SignedSize>(i);
    }
    return -1;
  }

  const char16_t* XMLAttributes::valueByIndex(Size i) const
  {
    if (i >= size()) return nullptr;
    if (!native_) return (*attributes_)[i].second.c_str();
    const auto& entry = native_->entries_[i];
    if (!entry.converted)
    {
      StringManager::convert(entry.value.data(), entry.value.size(), entry.value16);
      entry.converted = true;
    }
    return entry.value16.c_str();
  }

  void XMLAttributes::valueUTF8ByIndex_(Size i, std::string& value) const
  {
    if (native_) value.assign(native_->entries_[i].value);
    else value = StringManager::convert((*attributes_)[i].second.c_str());
  }

  bool XMLAttributes::valueUTF8(const char16_t* qname, std::string& value) const
  {
    const SignedSize i = index(qname);
    if (i < 0) return false;
    valueUTF8ByIndex_(static_cast<Size>(i), value);
    return true;
  }

  bool XMLAttributes::valueUTF8(const char* qname, std::string& value) const
  {
    const SignedSize i = indexUTF8_(qname);
    if (i < 0) return false;
    valueUTF8ByIndex_(static_cast<Size>(i), value);
    return true;
  }
}
