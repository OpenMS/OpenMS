// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

int pch_probe()
{
  std::vector<int> values {3, 1, 2};
  std::sort(values.begin(), values.end());
  return std::accumulate(values.begin(), values.end(), 0) + std::string("PCH").size();
}
