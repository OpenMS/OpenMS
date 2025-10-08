// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>

using namespace std;

namespace OpenMS
{
  MultiplexSatelliteCentroided::MultiplexSatelliteCentroided(size_t rt_idx, size_t mz_idx) :
    rt_idx_(rt_idx), mz_idx_(mz_idx)
  {
  }

  size_t MultiplexSatelliteCentroided::getMZidx() const
  {
    return mz_idx_;
  }

  size_t MultiplexSatelliteCentroided::getRTidx() const
  {
    return rt_idx_;
  }
  
}
