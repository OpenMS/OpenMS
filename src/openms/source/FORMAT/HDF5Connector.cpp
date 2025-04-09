// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HDF5Connector.h>

#include <OpenMS/CONCEPT/Exception.h>

#include "H5Cpp.h"

using namespace H5;

namespace OpenMS
{

  HDF5Connector::~HDF5Connector()
  {
    close();
  }

  void HDF5Connector::close()
  {
    if (file_)
    {
      file_->flush(H5F_SCOPE_LOCAL);
      file_->close();
      delete file_;
      file_ = nullptr;
    }
  }

  HDF5Connector::HDF5Connector(const String& filename, bool createNewFile)
  {
    // H5F_ACC_TRUNC - Truncate file, if it already exists, erasing all data previously stored in the file.
    // H5F_ACC_EXCL - Fail if file already exists. H5F_ACC_TRUNC and H5F_ACC_EXCL are mutually exclusive
    // H5F_ACC_RDONLY - Open file as read-only, if it already exists, and fail, otherwise
    // H5F_ACC_RDWR - Open file for read/write, if it already exists, and fail, otherwise
    unsigned int openFlag = H5F_ACC_RDWR;
    if (createNewFile)
    {
      openFlag = H5F_ACC_TRUNC;
    }
    FileCreatPropList fcparm = FileCreatPropList::DEFAULT;
    FileAccPropList faparm = FileAccPropList::DEFAULT;
    file_ = new H5::H5File(filename, openFlag, fcparm, faparm);
  }

}


