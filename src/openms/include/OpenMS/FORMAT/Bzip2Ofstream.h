// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: GitHub Copilot $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <bzlib.h>

namespace OpenMS
{
/**
    @brief Compresses files in the bzip2 format (*.bz2)
*/
  class OPENMS_DLLAPI Bzip2Ofstream
  {
public:
    /// Default Constructor
    Bzip2Ofstream();
    
    /// Detailed constructor with filename
    explicit Bzip2Ofstream(const char * filename);
    
    /// Destructor
    virtual ~Bzip2Ofstream();

    /**
      * @brief Writes n bytes from buffer s to the bzip2 compressed file
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
    /// pointer to a FILE object. Necessary for opening the file
    FILE * file_;
    
    /// a pointer to a BZFILE object. Necessary for compression
    BZFILE * bzip2file_;
    
    /// saves the last returned error by the write function
    int bzerror_;

    /// not implemented
    Bzip2Ofstream(const Bzip2Ofstream & bzip2);
    Bzip2Ofstream & operator=(const Bzip2Ofstream & bzip2);
  };

  inline bool Bzip2Ofstream::isOpen() const
  {
    return file_ != nullptr;
  }

} // namespace OpenMS
