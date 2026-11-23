// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/InputStream.h>
#include <string>

namespace OpenMS
{
  class GzipIfstream;

  /**
    * @brief Byte stream for reading gzip compressed XML files.
    * 
  */
  class OPENMS_DLLAPI GzipInputStream :
    public InputStream
  {
public:
    ///Constructor
    explicit GzipInputStream(const std::string& file_name);

    explicit GzipInputStream(const char* const file_name);

    ///Destructor
    ~GzipInputStream() override;

    ///returns true if file is open
    bool getIsOpen() const;

    /**
      * @brief returns the current position in the file
      *
      * @note OpenMS byte stream interface
    */
    UInt64 curPos() const override;

    /**
      * @brief writes bytes into buffer from file
      *
      * @note OpenMS byte stream interface
      *
      * @param[out] to_fill is the buffer which is written to
      * @param[in] max_to_read is the size of the buffer
      *
      * @return returns the number of bytes which were actually read
      *
    */
    Size readBytes(unsigned char* const to_fill, const Size max_to_read) override;


    GzipInputStream() = delete;
    GzipInputStream(const GzipInputStream& stream) = delete;
    GzipInputStream& operator=(const GzipInputStream& stream) = delete;

private:
    ///pointer to an compression stream
    GzipIfstream* gzip_ = nullptr;
    ///current index of the actual file
    UInt64 file_current_index_;
  };

  inline UInt64 GzipInputStream::curPos() const
  {
    return file_current_index_;
  }

} // namespace OpenMS

