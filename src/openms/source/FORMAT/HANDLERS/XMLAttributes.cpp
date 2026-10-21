// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>

#include <OpenMS/FORMAT/HANDLERS/StringManager.h>

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
  }

  const char16_t* XMLAttributes::value(const char16_t* qname) const
  {
    const SignedSize i = index(qname);
    return i < 0 ? nullptr : (*attributes_)[i].second.c_str();
  }

  const char16_t* XMLAttributes::value(const char* qname) const
  {
    if (!qname) return nullptr;
    bool ascii = true;
    for (const auto* ch = qname; *ch; ++ch)
    {
      if (static_cast<unsigned char>(*ch) > 127) { ascii = false; break; }
    }
    if (ascii)
    {
      // XML vocabulary names are usually ASCII; avoid allocating a transcoder.
      for (const auto& item : *attributes_)
      {
        if (equalsASCII(item.first.c_str(), qname)) return item.second.c_str();
      }
      return nullptr;
    }
    const auto name = StringManager::convert(qname);
    return value(name.c_str());
  }

  SignedSize XMLAttributes::index(const char16_t* qname) const
  {
    if (!qname) return -1;
    for (Size i = 0; i < attributes_->size(); ++i)
    {
      if ((*attributes_)[i].first == qname) return static_cast<SignedSize>(i);
    }
    return -1;
  }

  const char16_t* XMLAttributes::valueByIndex(Size i) const
  {
    return i < attributes_->size() ? (*attributes_)[i].second.c_str() : nullptr;
  }
}
