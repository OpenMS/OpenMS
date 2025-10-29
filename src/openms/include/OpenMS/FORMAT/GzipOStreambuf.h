// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/GzipOfstream.h>
#include <streambuf>
#include <ostream>

namespace OpenMS
{
  /**
    @brief std::streambuf implementation for gzip compression
    
    This class provides a std::streambuf interface that writes to a GzipOfstream.
    It can be used to create a std::ostream that automatically compresses data.
  */
  class OPENMS_DLLAPI GzipOStreambuf : public std::streambuf
  {
  public:
    /// Constructor
    explicit GzipOStreambuf(const char* filename);
    
    /// Destructor
    virtual ~GzipOStreambuf();
    
    /// Returns whether the stream is open
    bool isOpen() const;
    
  protected:
    /// Called when the buffer is full or when sync() is called
    virtual int_type overflow(int_type c = traits_type::eof()) override;
    
    /// Synchronize the stream buffer
    virtual int sync() override;
    
  private:
    GzipOfstream gzfile_;
    static const size_t buffer_size_ = 4096;
    char buffer_[4096];
    
    /// Flush the buffer to the gzip file
    void flushBuffer_();
  };

} // namespace OpenMS
