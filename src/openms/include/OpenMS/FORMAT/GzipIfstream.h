// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <zlib.h>

namespace OpenMS
{
/**
    @brief Decompresses files which are compressed in the gzip format (*.gzip)
*/
  class OPENMS_DLLAPI GzipIfstream
  {
public:
    ///Default Constructor
    GzipIfstream();

    /// Detailed constructor with filename
    explicit GzipIfstream(const char * filename);

    ///Destructor
    virtual ~GzipIfstream();

    /**
      * @brief Reads n bytes from the bzip2 compressed file into buffer s
      * 
      * @param s Buffer to be filled with the output 
      * @param n The size of the buffer s
      * @return The number of actually read bytes. If it is less than n, the end of the file was reached and the stream is closed
      * 
      * @note This returns a raw byte stream that is *not* null-terminated. Be careful here.
      * @note The length of the buffer needs to at least n
      * @note Closes the stream if the end of file is reached. Check isOpen before reading from the file again
      * 
      * @exception Exception::ConversionError is thrown if decompression fails
      * @exception Exception::IllegalArgument is thrown if no file for decompression is given. This can happen even happen if a file was already open but read until the end.
    */
    size_t read(char * s, size_t n);

    /**
      * @brief indicates whether the read function can be used safely
      *
      * @return true if end of file was reached. Otherwise false.
    */
    bool streamEnd() const;

    /**
      * @brief returns whether a file is open.
    */
    bool isOpen() const;

    /**
      * @brief opens a file for reading (decompression)
      *
      * @note any previous open files will be closed first!
    */
    void open(const char * filename);

    /**
      * @brief closes current file.
    */
    void close();

    /*
        @brief updates crc32 check sum whether the buffer is corrupted
        @note if this function is used it has to be called after every call of function read
        @param s the buffer which will be checked
        @param n the size of the buffer
    *
    //void updateCRC32(const char* s,const size_t n);

    *
        @brief	checks if data is corrupted after crc32 was computed
        @note   it can only be used if updateCRC32 was called after every call of function read
        @return true if the buffer and hence the file is corrupted; no decompression is possible
    *
    //bool isCorrupted();

    //unsigned long Crc32_ComputeBuf( unsigned long inCrc32, const void *buf,
//                                 size_t bufLen );*/

protected:

    ///a gzFile object(void*) . Necessary for decompression
    gzFile gzfile_;
    ///counts the last read duffer
    int n_buffer_;
    ///saves the last returned error by the read function
    int gzerror_;
    ///true if end of file is reached
    bool stream_at_end_;

    //needed if one wants to know whether file is okay
    //unsigned long original_crc;
    //needed if one wants to know whether file is okay
    //unsigned long crc;

    ///not implemented
    GzipIfstream(const GzipIfstream & bzip2);
    GzipIfstream & operator=(const GzipIfstream & bzip2);
  };

  inline bool GzipIfstream::isOpen() const
  {
    return gzfile_ != nullptr;
  }

  inline bool GzipIfstream::streamEnd() const
  {
    return stream_at_end_;
  }

/*	inline bool GzipIfstream::isCorrupted()
    {
        std::cout<<"CRC"<<crc<<std::endl;
        return (crc != original_crc);
    }*/

} //namespace OpenMS
