// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <QColor>

namespace OpenMS
{
  namespace
  {
    Annotation1DPeakItem<Peak1D> p(Peak1D(0, 0), "test", QColor());
  }
}

