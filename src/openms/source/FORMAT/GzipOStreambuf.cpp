// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/GzipOStreambuf.h>
#include <cstring>

namespace OpenMS
{
  GzipOStreambuf::GzipOStreambuf(const char* filename) :
    gzfile_(filename)
  {
    // Set the buffer for the streambuf
    setp(buffer_, buffer_ + buffer_size_ - 1);
  }

  GzipOStreambuf::~GzipOStreambuf()
  {
    sync(); // Flush any remaining data
  }

  bool GzipOStreambuf::isOpen() const
  {
    return gzfile_.isOpen();
  }

  GzipOStreambuf::int_type GzipOStreambuf::overflow(int_type c)
  {
    if (c != traits_type::eof())
    {
      *pptr() = static_cast<char>(c);
      pbump(1);
    }
    
    if (flushBuffer_() == -1)
    {
      return traits_type::eof();
    }
    
    return c;
  }

  int GzipOStreambuf::sync()
  {
    if (flushBuffer_() == -1)
    {
      return -1;
    }
    return 0;
  }

  void GzipOStreambuf::flushBuffer_()
  {
    int num = static_cast<int>(pptr() - pbase());
    if (num > 0)
    {
      gzfile_.write(pbase(), num);
      pbump(-num); // Reset put pointer
    }
  }

} // namespace OpenMS
