// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/GzipOfstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  GzipOfstream::GzipOfstream(const char * filename) :
    gzfile_(nullptr)
  {
    open(filename);
  }

  GzipOfstream::GzipOfstream() :
    gzfile_(nullptr)
  {
  }

  GzipOfstream::~GzipOfstream()
  {
    close();
  }

  size_t GzipOfstream::write(const char * s, size_t n)
  {
    if (gzfile_ != nullptr)
    {
      int bytes_written = gzwrite(gzfile_, s, (unsigned int) n);
      if (bytes_written < 0)
      {
        int gzerror_code = 0;
        const char* err_string = gzerror(gzfile_, &gzerror_code);
        std::string error_message = err_string ? err_string : "unknown error (code: " + std::to_string(gzerror_code) + ")";
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "error writing to gzip file: " + error_message);
      }
      return bytes_written;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no file for compression initialized");
    }
  }

  void GzipOfstream::open(const char * filename)
  {
    if (gzfile_ != nullptr)
    {
      close();
    }
    
    // Use compression level 6 (default) with binary mode
    gzfile_ = gzopen(filename, "wb6");
    gzbuffer(gzfile_, 524288); // set buffer size to 512kb to improve write performance
    
    if (gzfile_ == nullptr)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void GzipOfstream::close()
  {
    if (gzfile_ != nullptr)
    {
      gzclose(gzfile_);
      gzfile_ = nullptr;
    }
  }

} // namespace OpenMS
