// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace xercesc;

namespace OpenMS
{
  Bzip2InputStream::Bzip2InputStream(const String & file_name) :
    bzip2_(new Bzip2Ifstream(file_name.c_str())), file_current_index_(0)
  {
  }

  Bzip2InputStream::Bzip2InputStream(const char * file_name) :
    bzip2_(new Bzip2Ifstream(file_name)), file_current_index_(0)
  {
  }

/*	Bzip2InputStream::Bzip2InputStream()
    :bzip2_(NULL)
    {

    }*/

  Bzip2InputStream::~Bzip2InputStream()
  {
    delete bzip2_;
  }

  XMLSize_t Bzip2InputStream::readBytes(XMLByte * const to_fill, const XMLSize_t max_to_read)
  {
    //  Figure out whether we can really read.
    if (bzip2_->streamEnd())
    {
      return 0;
    }

    unsigned char * fill_it = static_cast<unsigned char *>(to_fill);
    XMLSize_t actual_read = (XMLSize_t) bzip2_->read((char *)fill_it, static_cast<size_t>(max_to_read));
    file_current_index_ += actual_read;
    return actual_read;
  }

  const XMLCh * Bzip2InputStream::getContentType() const
  {
    return nullptr;
  }

} // namespace OpenMS
