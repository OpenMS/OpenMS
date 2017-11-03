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

#ifndef OPENMS_FORMAT_GZIPINPUTSTREAM_H
#define OPENMS_FORMAT_GZIPINPUTSTREAM_H

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/GzipIfstream.h>

#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/PlatformUtils.hpp>

namespace OpenMS
{
  class String;

  /**
    * @brief Implements the BinInputStream class of the xerces-c library in order to read gzip compressed XML files.
    * 
  */
  class OPENMS_DLLAPI GzipInputStream :
    public xercesc::BinInputStream
  {
public:
    ///Constructor
    explicit GzipInputStream(const String& file_name);

    explicit GzipInputStream(const char* const file_name);

    ///Destructor
    virtual ~GzipInputStream();

    ///returns true if file is open
    bool getIsOpen() const;

    /**
      * @brief returns the current position in the file
      *
      * @note Implementation of the xerces-c input stream interface
    */
    virtual XMLFilePos curPos() const;

    /**
      * @brief writes bytes into buffer from file
      *
      * @note Implementation of the xerces-c input stream interface
      *
      * @param to_fill is the buffer which is written to
      * @param max_to_read is the size of the buffer
      *
      * @return returns the number of bytes which were actually read
      *
    */
    virtual XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read);

    /**
      * @brief returns 0
      *
      * @note Implementation of the xerces-c input stream interface
      *
      * If no content type is provided for the data, 0 is returned (as is the
      * case here, see xerces docs).
      *
      *
    */
    virtual const XMLCh* getContentType() const;


private:
    ///pointer to an compression stream
    GzipIfstream* gzip_;
    ///current index of the actual file
    XMLSize_t file_current_index_;

    //not implemented
    GzipInputStream();
    GzipInputStream(const GzipInputStream& stream);
    GzipInputStream& operator=(const GzipInputStream& stream);
  };

  inline XMLFilePos GzipInputStream::curPos() const
  {
    return file_current_index_;
  }

  inline bool GzipInputStream::getIsOpen() const
  {
    return gzip_->isOpen();
  }

} // namespace OpenMS

#endif // OPENMS_FORMAT_GZIPInputStream_H
