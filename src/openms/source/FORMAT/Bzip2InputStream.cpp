// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <algorithm>
#include <limits>


namespace OpenMS
{
  Bzip2InputStream::Bzip2InputStream(const std::string & file_name) :
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

  Size Bzip2InputStream::readBytes(unsigned char * const to_fill, const Size max_to_read)
  {
    if (max_to_read == 0) { return 0; }
    //  Figure out whether we can really read.
    if (bzip2_->streamEnd())
    {
      return 0;
    }

    const Size limit = std::min(max_to_read, static_cast<Size>(std::numeric_limits<int>::max()));
    const Size actual_read = bzip2_->read(reinterpret_cast<char*>(to_fill), limit);
    file_current_index_ += actual_read;
    return actual_read;
  }


} // namespace OpenMS
