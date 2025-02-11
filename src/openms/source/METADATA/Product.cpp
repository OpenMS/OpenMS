// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/Product.h>

using namespace std;

namespace OpenMS
{

  bool Product::operator==(const Product & rhs) const
  {
    return mz_ == rhs.mz_ &&
           window_low_ == rhs.window_low_ &&
           window_up_ == rhs.window_up_ &&
           CVTermList::operator==(rhs);
  }

  bool Product::operator!=(const Product & rhs) const
  {
    return !(operator==(rhs));
  }

  double Product::getMZ() const
  {
    return mz_;
  }

  void Product::setMZ(double mz)
  {
    mz_ = mz;
  }

  double Product::getIsolationWindowLowerOffset() const
  {
    return window_low_;
  }

  void Product::setIsolationWindowLowerOffset(double bound)
  {
    window_low_ = bound;
  }

  double Product::getIsolationWindowUpperOffset() const
  {
    return window_up_;
  }

  void Product::setIsolationWindowUpperOffset(double bound)
  {
    window_up_ = bound;
  }

}

