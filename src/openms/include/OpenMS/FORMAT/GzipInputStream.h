// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

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
    ~GzipInputStream() override;

    ///returns true if file is open
    bool getIsOpen() const;

    /**
      * @brief returns the current position in the file
      *
      * @note Implementation of the xerces-c input stream interface
    */
    XMLFilePos curPos() const override;

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
    XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read) override;

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
    const XMLCh* getContentType() const override;

    GzipInputStream() = delete;
    GzipInputStream(const GzipInputStream& stream) = delete;
    GzipInputStream& operator=(const GzipInputStream& stream) = delete;

private:
    ///pointer to an compression stream
    GzipIfstream* gzip_ = nullptr;
    ///current index of the actual file
    XMLSize_t file_current_index_;
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

