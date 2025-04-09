// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Guillaume Belz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>

using namespace std;

#ifdef OPENMS_BIG_ENDIAN
template <typename T>
T ByteReverse(const T in)
{
  T out;
  const char * pin = (const char *) &in;
  char * pout = (char *) (&out + 1) - 1;

  int i;
  for (i = sizeof(T); i > 0; --i)
  {
    *pout-- = *pin++;
  }
  return out;
}

#endif

namespace OpenMS::Internal
{
  FidHandler::FidHandler(const String & filename) :
    ifstream(filename.c_str(), ios_base::binary | ios_base::in)
  {
    index_ = 0;
    seekg(0, ios::beg);
  }

  FidHandler::~FidHandler() = default;

  Size FidHandler::getIndex() const
  {
    return index_;
  }

  Size FidHandler::getIntensity()
  {
    // intensity is coded in 32 bits little-endian integer format
    Int32 result = 0;
    read((char *) &result, 4);
#ifdef OPENMS_BIG_ENDIAN
    result = ByteReverse<Int32>(result);
#endif
    index_++;
    return (result > 0) ? result : 0;
  }
} // namespace OpenMS  // namespace Internal
