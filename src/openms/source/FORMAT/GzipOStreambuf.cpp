// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/GzipOStreambuf.h>

namespace OpenMS
{
  GzipOStreambuf::GzipOStreambuf(const char* filename) :
    gzfile_(filename)
  {
    // Set the buffer for the streambuf
    setp(buffer_, buffer_ + buffer_size_);
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
    
    flushBuffer_();
    
    return c;
  }

  int GzipOStreambuf::sync()
  {
    flushBuffer_();
    return 0;
  }

  void GzipOStreambuf::flushBuffer_()
  {
    std::ptrdiff_t num = pptr() - pbase();
    if (num > 0)
    {
      gzfile_.write(pbase(), static_cast<size_t>(num));
      setp(buffer_, buffer_ + buffer_size_); // Reset buffer pointers
    }
  }

} // namespace OpenMS
