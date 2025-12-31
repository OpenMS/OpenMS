// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h> // for WITH_PARQUET and OPENMS_DLLAPI

#ifdef WITH_PARQUET

#include <memory>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{
  class MSExperiment;

  /**
    @brief Export MSExperiment to Apache Arrow Table
  */
  OPENMS_DLLAPI std::shared_ptr<arrow::Table> exportMSExperimentToArrow(const MSExperiment& exp);

} // namespace OpenMS

#endif 


