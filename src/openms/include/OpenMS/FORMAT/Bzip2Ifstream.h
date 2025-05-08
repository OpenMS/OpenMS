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
    @brief Decompresses files which are compressed in the bzip2 format (*.bz2)
*/
  class OPENMS_DLLAPI Bzip2Ifstream
  {
public:
    ///Default Constructor
    Bzip2Ifstream();
    /// Detailed constructor with filename
    explicit Bzip2Ifstream(const char * filename);
    ///Destructor
    virtual ~Bzip2Ifstream();

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
      * @note any previous open files will be closed first!
    */
    void open(const char * filename);

    /**
      * @brief closes current file.
    */
    void close();

protected:
    /// pointer to a FILE object. Necessary for opening the file
    FILE * file_;
    /// a pointer to a BZFILE object. Necessary for decompression
    BZFILE * bzip2file_;
    ///counts the last read buffer
    size_t     n_buffer_;
    ///saves the last returned error by the read function
    int     bzerror_;
    ///true if end of file is reached
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
