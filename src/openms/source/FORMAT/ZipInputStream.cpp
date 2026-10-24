// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZipInputStream.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ZipIfstream.h>
#include <algorithm>
#include <limits>


namespace OpenMS
{
  ZipInputStream::ZipInputStream(const std::string& file_name) :
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

  Size ZipInputStream::readBytes(unsigned char* const to_fill, const Size max_to_read)
  {
    if (max_to_read == 0) { return 0; }
    if (zip_->streamEnd())
    {
      return 0;
    }

    const Size limit = std::min(max_to_read, static_cast<Size>(std::numeric_limits<int>::max()));
    const Size actual_read = zip_->read(reinterpret_cast<char*>(to_fill), limit);
    file_current_index_ += actual_read;
    return actual_read;
  }


} // namespace OpenMS
