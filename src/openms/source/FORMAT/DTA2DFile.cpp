// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DTA2DFile.h>

using namespace std;

namespace OpenMS
{

  DTA2DFile::DTA2DFile() = default;

  DTA2DFile::~DTA2DFile() = default;

  PeakFileOptions & DTA2DFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & DTA2DFile::getOptions() const
  {
    return options_;
  }

} // namespace OpenMS
