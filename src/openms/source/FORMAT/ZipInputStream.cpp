// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZipInputStream.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/ZipIfstream.h>

using namespace xercesc;

namespace OpenMS
{
  ZipInputStream::ZipInputStream(const String& file_name) :
    zip_(new ZipIfstream(file_name.c_str())), file_current_index_(0)
  {
  }

  ZipInputStream::ZipInputStream(const char* file_name) :
    zip_(new ZipIfstream(file_name)), file_current_index_(0)
  {
  }

  ZipInputStream::~ZipInputStream()
  {
    delete zip_;
  }

  bool ZipInputStream::getIsOpen() const
  {
    return zip_->isOpen();
  }

  XMLSize_t ZipInputStream::readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read)
  {
    if (zip_->streamEnd())
    {
      return 0;
    }

    unsigned char* fill_it = static_cast<unsigned char*>(to_fill);
    XMLSize_t actual_read = (XMLSize_t)zip_->read((char*)fill_it, static_cast<size_t>(max_to_read));
    file_current_index_ += actual_read;
    return actual_read;
  }

  const XMLCh* ZipInputStream::getContentType() const
  {
    return nullptr;
  }

} // namespace OpenMS
