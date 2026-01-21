// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenSwath
{
  struct LightTargetedExperiment;
}

namespace OpenMS
{

  /**
      @brief Read OpenSwath Parquet library input (.pqp_parquet) into LightTargetedExperiment.

      The Parquet library format is a directory container with separate tables for
      precursors and transitions. The MVP reader materializes all rows into
      OpenSwath::LightTargetedExperiment.
  */
  class OPENMS_DLLAPI TransitionParquetFile
  {
  public:
    /// Default constructor
    TransitionParquetFile() = default;

    /// Default destructor
    ~TransitionParquetFile() = default;

    /// Read a .pqp_parquet library directory and populate a LightTargetedExperiment
    void convertParquetToTargetedExperiment(const String& pqp_parquet_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const;
  };

} // namespace OpenMS
