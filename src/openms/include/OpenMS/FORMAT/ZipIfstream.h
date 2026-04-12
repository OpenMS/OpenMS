// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstddef>

namespace OpenMS
{
/**
    @brief Decompresses a single-entry ZIP archive for streaming reading.

    Wraps libzip to stream bytes from the single non-directory entry in a ZIP file.
    Directory entries (names ending with '/') are silently skipped. Exactly one
    non-directory entry must be present; otherwise open() throws.
*/
  class OPENMS_DLLAPI ZipIfstream
  {
public:
    /// Default Constructor
    ZipIfstream();

    /// Detailed constructor with filename
    explicit ZipIfstream(const char* filename);

    /// Destructor
    virtual ~ZipIfstream();

    /**
      * @brief Reads n bytes from the zip entry into buffer s
      *
      * @param[out] s Buffer to be filled with the output
      * @param[in] n The size of the buffer s
      * @return The number of actually read bytes. If it is less than n, the end of the entry was reached and the stream is closed.
      *
      * @exception Exception::ParseError is thrown if decompression fails
      * @exception Exception::IllegalArgument is thrown if no file for decompression is given
    */
    size_t read(char* s, size_t n);

    /**
      * @brief indicates whether the read function can be used safely
      *
      * @return true if end of entry was reached. Otherwise false.
    */
    bool streamEnd() const;

    /**
      * @brief returns whether a file is open.
    */
    bool isOpen() const;

    /**
      * @brief opens a ZIP archive for reading from its single non-directory entry
      *
      * @note any previous open files will be closed first!
      * @exception Exception::FileNotFound if file does not exist
      * @exception Exception::FileNotReadable if file cannot be opened
      * @exception Exception::FileEmpty if ZIP contains no non-directory entries
      * @exception Exception::ParseError if ZIP contains more than one non-directory entry
    */
    void open(const char* filename);

    /**
      * @brief closes current file.
    */
    void close();

    ZipIfstream(const ZipIfstream&) = delete;
    ZipIfstream& operator=(const ZipIfstream&) = delete;
    ZipIfstream(ZipIfstream&&) = delete;
    ZipIfstream& operator=(ZipIfstream&&) = delete;

protected:
    /// opaque pointers to libzip handles (void* to avoid including zip.h in header)
    void* zip_archive_ = nullptr;
    void* zip_entry_ = nullptr;
    /// true if end of entry is reached
    bool stream_at_end_ = true;
  };

  inline bool ZipIfstream::isOpen() const
  {
    return zip_entry_ != nullptr;
  }

  inline bool ZipIfstream::streamEnd() const
  {
    return stream_at_end_;
  }

} // namespace OpenMS
