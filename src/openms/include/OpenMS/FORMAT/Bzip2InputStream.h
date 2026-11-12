// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/InputStream.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

namespace OpenMS
{
  /**
    * @brief Byte stream for reading bzip2 compressed XML files.
    *
  */
  class OPENMS_DLLAPI Bzip2InputStream :
    public InputStream
  {
public:
    ///Constructor
    explicit Bzip2InputStream(const std::string& file_name);

    explicit Bzip2InputStream(const char* const file_name);

    ///Destructor
    ~Bzip2InputStream() override;

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
    Size readBytes(unsigned char* const  to_fill, const Size max_to_read) override;



private:
    ///pointer to an compression stream
    Bzip2Ifstream* bzip2_;
    ///current index of the actual file
    UInt64     file_current_index_;

    //not implemented
    Bzip2InputStream();
    Bzip2InputStream(const Bzip2InputStream& stream);
    Bzip2InputStream& operator=(const Bzip2InputStream& stream);
  };

  inline UInt64 Bzip2InputStream::curPos() const
  {
    return file_current_index_;
  }

  inline bool Bzip2InputStream::getIsOpen() const
  {
    return bzip2_->isOpen();
  }

} // namespace OpenMS

