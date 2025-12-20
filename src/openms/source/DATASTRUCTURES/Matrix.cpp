// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS
{
  // Explicit template instantiations to ensure symbols exist in the library.
  // This compiles the template code for these types into libOpenMS.so.
  template class Matrix<int>;
  template class Matrix<double>;
  template class Matrix<float>;
}
