// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/GzipInputStream.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/GzipIfstream.h>

using namespace xercesc;

namespace OpenMS
{
  GzipInputStream::GzipInputStream(const String & file_name) :
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

  XMLSize_t GzipInputStream::readBytes(XMLByte * const to_fill, const XMLSize_t max_to_read)
  {
    // Figure out whether we can really read.
    if (gzip_->streamEnd())
    {
      return 0;
    }

    unsigned char * fill_it = static_cast<unsigned char *>(to_fill);
    XMLSize_t actual_read = (XMLSize_t) gzip_->read((char *)fill_it, static_cast<size_t>(max_to_read));
    file_current_index_ += actual_read;
    return actual_read;
  }

  const XMLCh * GzipInputStream::getContentType() const
  {
    return nullptr;
  }

} // namespace OpenMS
