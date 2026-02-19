// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureMapArrowIO.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>

#include <arrow/api.h>
#include <arrow/builder.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>
#include <parquet/properties.h>

#include <filesystem>

namespace OpenMS
{

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportFeaturesToArrow(
  const FeatureMap& /*feature_map*/)
{
  return nullptr;
}

std::shared_ptr<arrow::Table> FeatureMapArrowIO::exportPSMsToArrow(
  const FeatureMap& /*feature_map*/)
{
  return nullptr;
}

bool FeatureMapArrowIO::exportToParquet(
  const FeatureMap& /*feature_map*/,
  const String& /*directory*/,
  const ParquetWriteConfig& /*config*/)
{
  return false;
}

bool FeatureMapArrowIO::importFeaturesFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

bool FeatureMapArrowIO::importPSMsFromArrow(
  const std::shared_ptr<arrow::Table>& /*table*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

bool FeatureMapArrowIO::importFromParquet(
  const String& /*directory*/,
  FeatureMap& /*feature_map*/)
{
  return false;
}

} // namespace OpenMS

#endif // WITH_PARQUET
