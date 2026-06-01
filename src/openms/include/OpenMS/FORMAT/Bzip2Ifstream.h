// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <bzlib.h>
#include <istream>

namespace OpenMS
{
  /**
    @brief Streaming decompressor for bzip2 (@c .bz2) files.

    Reads raw bytes out of a bzip2-compressed file on demand. The class
    owns one underlying @c FILE handle and one bzip2 read handle; both are
    released by @ref close (also called from the destructor). Copy and
    assignment are deliberately suppressed.

    Typical lifecycle: construct (or call @ref open) with the path of a
    @c .bz2 file, repeatedly call @ref read until it returns fewer bytes
    than requested (at which point the stream is closed automatically,
    @ref isOpen returns @c false and @ref streamEnd returns @c true).

    @ingroup FileIO
  */
  class OPENMS_DLLAPI Bzip2Ifstream
  {
public:
    /// Default constructor; leaves the instance in a closed state (@ref isOpen returns @c false, @ref streamEnd returns @c true).
    Bzip2Ifstream();

    /**
      @brief Construct and open @p filename for reading in one step.

      @param[in] filename Path of the bzip2-compressed file to open
                          (binary mode is used so the contents are not
                          translated on Windows).

      @throws Exception::FileNotFound    when @p filename does not exist.
      @throws Exception::FileNotReadable when @p filename exists but cannot be read.
      @throws Exception::IOException     when opening the file fails for any other reason.
      @throws Exception::ConversionError when initialising the bzip2 reader on the opened file fails.
    */
    explicit Bzip2Ifstream(const char * filename);

    /// Destructor; closes the bzip2 reader and the underlying file handle.
    virtual ~Bzip2Ifstream();

    /**
      @brief Read up to @p n decompressed bytes into @p s.

      Returns the actual number of bytes written into @p s. A return value
      smaller than @p n indicates that the end of the compressed stream was
      reached during this call; the stream is then closed automatically
      (@ref isOpen becomes @c false and @ref streamEnd becomes @c true).

      @param[out] s Buffer that receives the decompressed bytes; must be at
                    least @p n bytes long.
      @param[in]  n Number of bytes to attempt to read.
      @return Number of decompressed bytes actually written into @p s.

      @note The output is a raw byte stream and is @b not null-terminated.
      @note When the return value is less than @p n the stream has been
            closed; check @ref isOpen before issuing further reads.

      @throws Exception::ParseError       when bzip2 reports a decompression failure mid-stream.
      @throws Exception::IllegalArgument  when called on an instance that has no open file (default-constructed, or already closed -- including after the previous call read to end of file).
    */
    size_t read(char * s, size_t n);

    /**
      @brief Whether the stream has been closed (end of file reached or no file opened).

      @return @c true if no further @ref read calls should be issued, @c false otherwise.
    */
    bool streamEnd() const;

    /**
      @brief Whether a file is currently open for reading.

      @return @c true if a file is open, @c false otherwise.
    */
    bool isOpen() const;

    /**
      @brief Open @p filename for reading; closes the currently open file first.

      @param[in] filename Path of the bzip2-compressed file to open.

      @throws Exception::FileNotFound    when @p filename does not exist.
      @throws Exception::FileNotReadable when @p filename exists but cannot be read.
      @throws Exception::IOException     when opening the file fails for any other reason.
      @throws Exception::ConversionError when initialising the bzip2 reader on the opened file fails.
    */
    void open(const char * filename);

    /// Close the currently open file (no-op if none is open). After the call, @ref isOpen returns @c false and @ref streamEnd returns @c true.
    void close();

protected:
    /// Underlying @c FILE handle used to feed the bzip2 reader; @c nullptr when no file is open.
    FILE * file_;
    /// bzip2 read handle attached to @ref file_; @c nullptr when no file is open.
    BZFILE * bzip2file_;
    /// Number of bytes produced by the most recent @ref read call.
    size_t     n_buffer_;
    /// Most recent bzip2 status code returned to the read API.
    int     bzerror_;
    /// @c true once the stream has been closed (either explicitly via @ref close or implicitly because the end of the compressed file was reached).
    bool stream_at_end_;

    //not implemented
    Bzip2Ifstream(const Bzip2Ifstream & bzip2);
    Bzip2Ifstream & operator=(const Bzip2Ifstream & bzip2);
  };

  //return bzip2file???!!!!????
  inline bool Bzip2Ifstream::isOpen() const
  {
    return file_ != nullptr;
  }

  inline bool Bzip2Ifstream::streamEnd() const
  {
    return stream_at_end_;
  }

} //namespace OpenMS
