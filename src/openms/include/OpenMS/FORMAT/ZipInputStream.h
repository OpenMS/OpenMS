// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <xercesc/util/BinInputStream.hpp>

namespace OpenMS
{
  class ZipIfstream;
  class String;

  /**
    * @brief Implements the BinInputStream class of the xerces-c library in order to read ZIP compressed XML files.
  */
  class OPENMS_DLLAPI ZipInputStream :
    public xercesc::BinInputStream
  {
public:
    /// Constructor
    explicit ZipInputStream(const String& file_name);

    explicit ZipInputStream(const char* const file_name);

    /// Destructor
    ~ZipInputStream() override;

    /// returns true if file is open
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
      * @param[out] to_fill is the buffer which is written to
      * @param[in] max_to_read is the size of the buffer
      *
      * @return returns the number of bytes which were actually read
    */
    XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read) override;

    /**
      * @brief returns 0
      *
      * @note Implementation of the xerces-c input stream interface
    */
    const XMLCh* getContentType() const override;

    ZipInputStream() = delete;
    ZipInputStream(const ZipInputStream& stream) = delete;
    ZipInputStream& operator=(const ZipInputStream& stream) = delete;

private:
    /// pointer to the zip decompression stream
    ZipIfstream* zip_ = nullptr;
    /// current index of the actual file
    XMLSize_t file_current_index_;
  };

  inline XMLFilePos ZipInputStream::curPos() const
  {
    return file_current_index_;
  }

} // namespace OpenMS
