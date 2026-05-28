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
    @brief Streaming reader for the single non-directory entry of a ZIP archive.

    Backed by libzip. Directory entries (names ending in @c '/') are
    silently ignored when scanning the archive; the archive must
    contain exactly one non-directory entry, or @ref open will throw.
    Bytes are then streamed lazily from that entry; the stream auto-
    closes once the entry has been fully consumed.

    @note When OpenMS is built without libzip, every call that
          actually touches an archive (@ref open and consequently the
          single-argument constructor, plus @ref read) throws
          @c Exception::NotImplemented. The default constructor,
          @ref close, @ref isOpen and @ref streamEnd remain usable.

    Copy and move are deliberately suppressed.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ZipIfstream
  {
public:
    /// Default constructor; leaves the instance in a closed state (@ref isOpen returns @c false, @ref streamEnd returns @c true).
    ZipIfstream();

    /**
      @brief Construct and open @p filename in one step; see @ref open for the exact behaviour and exception set.

      @param[in] filename Path of the ZIP archive to open.
    */
    explicit ZipIfstream(const char* filename);

    /// Destructor; closes the archive and the entry handle.
    virtual ~ZipIfstream();

    /**
      @brief Read up to @p n decompressed bytes from the open entry into @p s.

      Returns the actual number of bytes written into @p s. A return
      value smaller than @p n indicates that the end of the entry was
      reached during this call and that the stream has been closed
      automatically (@ref isOpen becomes @c false and @ref streamEnd
      becomes @c true). Passing @c n == @c 0 is a no-op that returns
      @c 0 without touching the entry.

      @param[out] s Buffer that receives the decompressed bytes;
                    must be at least @p n bytes long.
      @param[in]  n Number of bytes to attempt to read.
      @return Number of decompressed bytes actually written into @p s.

      @throws Exception::ParseError      when @c zip_fread reports a
                                         decompression error.
      @throws Exception::IllegalArgument when called on an instance
                                         that has no open entry
                                         (default-constructed, or
                                         already closed -- including
                                         after a previous @ref read
                                         drained the entry).
      @throws Exception::NotImplemented  when OpenMS was built without
                                         libzip.
    */
    size_t read(char* s, size_t n);

    /**
      @brief Whether the stream has been closed (end of entry reached or no archive opened).

      @return @c true if no further @ref read calls should be issued, @c false otherwise.
    */
    bool streamEnd() const;

    /**
      @brief Whether an archive entry is currently open for reading.

      @return @c true if an entry is open, @c false otherwise.
    */
    bool isOpen() const;

    /**
      @brief Open @p filename and select its single non-directory entry for reading.

      Closes any previously open archive first, scans the new archive
      for entries and selects the single non-directory one as the
      streaming source.

      @param[in] filename Path of the ZIP archive to open.

      @throws Exception::FileNotFound    when @p filename does not exist.
      @throws Exception::FileNotReadable when @p filename exists but
                                         the archive cannot be opened
                                         (or the selected entry cannot
                                         be read).
      @throws Exception::FileEmpty       when the archive contains no
                                         non-directory entries.
      @throws Exception::ParseError      when the archive contains more
                                         than one non-directory entry.
                                         The exception message lists
                                         the conflicting entry names.
      @throws Exception::NotImplemented  when OpenMS was built without
                                         libzip.
    */
    void open(const char* filename);

    /// Close the currently open archive (no-op if none is open). After the call, @ref isOpen returns @c false and @ref streamEnd returns @c true.
    void close();

    ZipIfstream(const ZipIfstream&) = delete;
    ZipIfstream& operator=(const ZipIfstream&) = delete;
    ZipIfstream(ZipIfstream&&) = delete;
    ZipIfstream& operator=(ZipIfstream&&) = delete;

protected:
    void* zip_archive_ = nullptr; ///< Opaque libzip @c zip_t* handle; @c nullptr when no archive is open. Typed as @c void* to keep @c zip.h out of the header.
    void* zip_entry_ = nullptr;   ///< Opaque libzip @c zip_file_t* handle for the currently streamed entry; @c nullptr when no entry is open.
    bool stream_at_end_ = true;   ///< @c true once the stream has been closed (either explicitly via @ref close or implicitly because the entry was fully consumed by @ref read).
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
