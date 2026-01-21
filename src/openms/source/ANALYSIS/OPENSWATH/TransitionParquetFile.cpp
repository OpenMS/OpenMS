// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  void TransitionParquetFile::convertParquetToTargetedExperiment(
    const String& /*pqp_parquet_dir*/, OpenSwath::LightTargetedExperiment& /*targeted_exp*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }
} // namespace OpenMS
