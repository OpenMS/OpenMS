// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_BZIP2IFSTREAM_H
#define OPENMS_FORMAT_BZIP2IFSTREAM_H

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
#endif //OPENMS_FORMAT_BZIP2IFSTREAM_H
