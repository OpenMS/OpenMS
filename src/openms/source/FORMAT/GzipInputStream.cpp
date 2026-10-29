// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/GzipInputStream.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/GzipIfstream.h>
#include <algorithm>
#include <limits>


namespace OpenMS
{
  GzipInputStream::GzipInputStream(const std::string & file_name) :
    gzip_(new GzipIfstream(file_name.c_str())), file_current_index_(0)
  {
  }

  GzipInputStream::GzipInputStream(const char * file_name) :
    gzip_(new GzipIfstream(file_name)), file_current_index_(0)
  {
  }

  GzipInputStream::~GzipInputStream()
  {
    delete gzip_;
  }

  bool GzipInputStream::getIsOpen() const
  {
    return gzip_->isOpen();
  }

  Size GzipInputStream::readBytes(unsigned char * const to_fill, const Size max_to_read)
  {
    if (max_to_read == 0) { return 0; }
    // Figure out whether we can really read.
    if (gzip_->streamEnd())
    {
      return 0;
    }

    const Size limit = std::min(max_to_read, static_cast<Size>(std::numeric_limits<int>::max()));
    const Size actual_read = gzip_->read(reinterpret_cast<char*>(to_fill), limit);
    file_current_index_ += actual_read;
    return actual_read;
  }


} // namespace OpenMS
