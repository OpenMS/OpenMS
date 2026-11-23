// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/InputStream.h>
#include <string>

namespace OpenMS
{
  class ZipIfstream;

  /**
    * @brief Byte stream for reading ZIP compressed XML files.
  */
  class OPENMS_DLLAPI ZipInputStream :
    public InputStream
  {
public:
    /// Constructor
    explicit ZipInputStream(const std::string& file_name);

    explicit ZipInputStream(const char* const file_name);

    /// Destructor
    ~ZipInputStream() override;

    /// returns true if file is open
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
    */
    Size readBytes(unsigned char* const to_fill, const Size max_to_read) override;


    ZipInputStream() = delete;
    ZipInputStream(const ZipInputStream& stream) = delete;
    ZipInputStream& operator=(const ZipInputStream& stream) = delete;

private:
    /// pointer to the zip decompression stream
    ZipIfstream* zip_ = nullptr;
    /// current index of the actual file
    UInt64 file_current_index_;
  };

  inline UInt64 ZipInputStream::curPos() const
  {
    return file_current_index_;
  }

} // namespace OpenMS
