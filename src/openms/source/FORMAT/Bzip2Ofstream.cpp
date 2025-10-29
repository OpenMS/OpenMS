// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Bzip2Ofstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  Bzip2Ofstream::Bzip2Ofstream(const char * filename) :
    file_(nullptr),
    bzip2file_(nullptr),
    bzerror_(0)
  {
    open(filename);
  }

  Bzip2Ofstream::Bzip2Ofstream() :
    file_(nullptr),
    bzip2file_(nullptr),
    bzerror_(0)
  {
  }

  Bzip2Ofstream::~Bzip2Ofstream()
  {
    close();
  }

  size_t Bzip2Ofstream::write(const char * s, size_t n)
  {
    if (file_ != nullptr && bzip2file_ != nullptr)
    {
      BZ2_bzWrite(&bzerror_, bzip2file_, const_cast<char*>(s), (int) n);
      
      if (bzerror_ == BZ_IO_ERROR)
      {
        BZ2_bzWriteClose(&bzerror_, bzip2file_, 0, nullptr, nullptr);
        fclose(file_);
        file_ = nullptr;
        bzip2file_ = nullptr;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IO error during bzip2 compression");
      }
      else if (bzerror_ != BZ_OK)
      {
        BZ2_bzWriteClose(&bzerror_, bzip2file_, 0, nullptr, nullptr);
        fclose(file_);
        file_ = nullptr;
        bzip2file_ = nullptr;
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "error during bzip2 compression");
      }
      
      return n; // bzip2 doesn't return the number of bytes written
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "no file for compression initialized");
    }
  }

  void Bzip2Ofstream::open(const char * filename)
  {
    if (file_ != nullptr)
    {
      close();
    }
    
    file_ = fopen(filename, "wb");
    if (file_ == nullptr)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    
    // blockSize100k: 9 gives the best compression (1-9 are valid)
    // verbosity: 0 is silent
    // workFactor: 0 is default (30)
    bzip2file_ = BZ2_bzWriteOpen(&bzerror_, file_, 9, 0, 0);
    
    if (bzerror_ != BZ_OK || bzip2file_ == nullptr)
    {
      if (file_ != nullptr)
      {
        fclose(file_);
        file_ = nullptr;
      }
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void Bzip2Ofstream::close()
  {
    if (bzip2file_ != nullptr)
    {
      BZ2_bzWriteClose(&bzerror_, bzip2file_, 0, nullptr, nullptr);
      bzip2file_ = nullptr;
    }
    
    if (file_ != nullptr)
    {
      fclose(file_);
      file_ = nullptr;
    }
  }

} // namespace OpenMS
