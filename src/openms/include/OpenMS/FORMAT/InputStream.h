// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS contributors $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  /** @brief Internal byte stream shared by the XML compression adapters. */
  class OPENMS_DLLAPI InputStream
  {
  public:
    virtual ~InputStream() = default;
    virtual UInt64 curPos() const = 0;
    virtual Size readBytes(unsigned char* to_fill, Size max_to_read) = 0;
    virtual const char16_t* getContentType() const { return nullptr; }
  };
}
