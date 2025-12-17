// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Matrix.h>

// Note: Matrix<T> is now a header-only template class.
// No explicit template instantiations needed since all code is in the header.

namespace OpenMS
{
  // Provide default instantiations to ensure symbols exist in the library
  // (mainly for backwards compatibility with existing binaries)
  template class Matrix<int>;
  template class Matrix<double>;
  template class Matrix<float>;
}
