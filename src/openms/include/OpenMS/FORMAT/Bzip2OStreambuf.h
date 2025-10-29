// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/Bzip2Ofstream.h>
#include <streambuf>
#include <ostream>

namespace OpenMS
{
  /**
    @brief std::streambuf implementation for bzip2 compression
    
    This class provides a std::streambuf interface that writes to a Bzip2Ofstream.
    It can be used to create a std::ostream that automatically compresses data.
  */
  class OPENMS_DLLAPI Bzip2OStreambuf : public std::streambuf
  {
  public:
    /// Constructor
    explicit Bzip2OStreambuf(const char* filename);
    
    /// Destructor
    virtual ~Bzip2OStreambuf();
    
    /// Returns whether the stream is open
    bool isOpen() const;
    
  protected:
    /// Called when the buffer is full or when sync() is called
    virtual int_type overflow(int_type c = traits_type::eof()) override;
    
    /// Synchronize the stream buffer
    virtual int sync() override;
    
  private:
    Bzip2Ofstream bz2file_;
    static const size_t buffer_size_ = 4096;
    char buffer_[4096];
    
    /// Flush the buffer to the bzip2 file
    void flushBuffer_();
  };

} // namespace OpenMS
