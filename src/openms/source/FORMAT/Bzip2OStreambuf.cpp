// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Bzip2OStreambuf.h>

namespace OpenMS
{
  Bzip2OStreambuf::Bzip2OStreambuf(const char* filename) :
    bz2file_(filename)
  {
    // Set the buffer for the streambuf
    setp(buffer_, buffer_ + buffer_size_);
  }

  Bzip2OStreambuf::~Bzip2OStreambuf()
  {
    // Flush any remaining data
    // Catch and suppress exceptions since destructors should not throw
    try
    {
      sync();
    }
    catch (...)
    {
      // Suppress exceptions in destructor
    }
  }

  bool Bzip2OStreambuf::isOpen() const
  {
    return bz2file_.isOpen();
  }

  Bzip2OStreambuf::int_type Bzip2OStreambuf::overflow(int_type c)
  {
    if (c != traits_type::eof())
    {
      *pptr() = static_cast<char>(c);
      pbump(1);
    }
    
    flushBuffer_();
    
    return c;
  }

  int Bzip2OStreambuf::sync()
  {
    flushBuffer_();
    return 0;
  }

  void Bzip2OStreambuf::flushBuffer_()
  {
    std::ptrdiff_t num = pptr() - pbase();
    if (num > 0)
    {
      // Check if file is still open before attempting to write
      if (bz2file_.isOpen())
      {
        bz2file_.write(pbase(), static_cast<size_t>(num));
      }
      setp(buffer_, buffer_ + buffer_size_); // Reset buffer pointers
    }
  }

} // namespace OpenMS
