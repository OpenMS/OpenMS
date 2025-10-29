// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <zlib.h>

namespace OpenMS
{
/**
    @brief Compresses files in the gzip format (*.gz)
*/
  class OPENMS_DLLAPI GzipOfstream
  {
public:
    /// Default Constructor
    GzipOfstream();

    /// Detailed constructor with filename
    explicit GzipOfstream(const char * filename);

    /// Destructor
    virtual ~GzipOfstream();

    /**
      * @brief Writes n bytes from buffer s to the gzip compressed file
      * 
      * @param s Buffer containing the data to be written
      * @param n The number of bytes to write from buffer s
      * @return The number of actually written bytes
      * 
      * @exception Exception::ConversionError is thrown if compression fails
      * @exception Exception::IllegalArgument is thrown if no file for compression is given
    */
    size_t write(const char * s, size_t n);

    /**
      * @brief returns whether a file is open
    */
    bool isOpen() const;

    /**
      * @brief opens a file for writing (compression)
      *
      * @note any previous open files will be closed first!
    */
    void open(const char * filename);

    /**
      * @brief closes current file
    */
    void close();

protected:

    /// a gzFile object (void*). Necessary for compression
    gzFile gzfile_;

    /// not implemented
    GzipOfstream(const GzipOfstream & gzip);
    GzipOfstream & operator=(const GzipOfstream & gzip);
  };

  inline bool GzipOfstream::isOpen() const
  {
    return gzfile_ != nullptr;
  }

} // namespace OpenMS
