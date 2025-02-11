// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>


namespace OpenMS
{
  class String;
  /**
    * @brief Implements the BinInputStream class of the xerces-c library in order to read bzip2 compressed XML files.
    *
  */
  class OPENMS_DLLAPI Bzip2InputStream :
    public xercesc::BinInputStream
  {
public:
    ///Constructor
    explicit Bzip2InputStream(const String& file_name);

    explicit Bzip2InputStream(const char* const file_name);


    ///Destructor
    ~Bzip2InputStream() override;

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
    XMLSize_t readBytes(XMLByte* const  to_fill, const XMLSize_t max_to_read) override;

    /**
      * @brief returns 0
      *
      * @note Implementation of the xerces-c input stream interface
      *
      * If no content type is provided for the data, 0 is returned (as is the
      * case here, see xerces docs).
      *
    */
    const XMLCh* getContentType() const override;


private:
    ///pointer to an compression stream
    Bzip2Ifstream* bzip2_;
    ///current index of the actual file
    XMLSize_t       file_current_index_;

    //not implemented
    Bzip2InputStream();
    Bzip2InputStream(const Bzip2InputStream& stream);
    Bzip2InputStream& operator=(const Bzip2InputStream& stream);
  };

  inline XMLFilePos Bzip2InputStream::curPos() const
  {
    return file_current_index_;
  }

  inline bool Bzip2InputStream::getIsOpen() const
  {
    return bzip2_->isOpen();
  }

} // namespace OpenMS

