// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

// forward declarations
namespace H5
{
  struct H5File;
}

namespace OpenMS
{
  /**
    @brief File adapter for HDF5 files

    This class contains certain helper functions to deal with HDF5 files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI HDF5Connector
  {
public:

    /// Constructor
    HDF5Connector(const String& filename, bool createNewFile = false);

    /// Destructor
    ~HDF5Connector();

    void close();

protected:

    H5::H5File* file_ = nullptr;

  };


} // namespace OpenMS

