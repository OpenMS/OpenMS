// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/FORMAT/ArrowExport.h>

#ifdef WITH_PARQUET

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <arrow/api.h>

namespace OpenMS
{

  std::shared_ptr<arrow::Table> exportMSExperimentToArrow(const MSExperiment& exp)
  {
    // Handle empty experiment: return empty but valid table
    if (exp.empty())
    {
      auto schema = arrow::schema({
        arrow::field("mz", arrow::float64()),
        arrow::field("intensity", arrow::float64()),
        arrow::field("rt", arrow::float64())
      });
      return arrow::Table::Make(schema, {});
    }

    // Define memory pool (default)
    arrow::MemoryPool* pool = arrow::default_memory_pool();

    // Builders for columns: mz, intensity, rt (all float64)
    arrow::DoubleBuilder mz_builder(pool);
    arrow::DoubleBuilder int_builder(pool);
    arrow::DoubleBuilder rt_builder(pool);

    arrow::Status st;

    // Iterate over spectra
    for (const auto& spectrum : exp)
    {
      double rt = spectrum.getRT();

      // Iterate over peaks in spectrum
      for (const auto& peak : spectrum)
      {
        st = mz_builder.Append(peak.getMZ());
        if (!st.ok()) return nullptr;

        st = int_builder.Append(peak.getIntensity());
        if (!st.ok()) return nullptr;

        st = rt_builder.Append(rt);
        if (!st.ok()) return nullptr;
      }
    }

    // Finalize arrays
    std::shared_ptr<arrow::Array> mz_array;
    std::shared_ptr<arrow::Array> int_array;
    std::shared_ptr<arrow::Array> rt_array;

    if (!mz_builder.Finish(&mz_array).ok()) return nullptr;
    if (!int_builder.Finish(&int_array).ok()) return nullptr;
    if (!rt_builder.Finish(&rt_array).ok()) return nullptr;

    // Create Schema
    auto schema = arrow::schema({
      arrow::field("mz", arrow::float64()),
      arrow::field("intensity", arrow::float64()),
      arrow::field("rt", arrow::float64())
    });

    // Create Table
    return arrow::Table::Make(schema, {mz_array, int_array, rt_array});
  }

} // namespace OpenMS

#endif // WITH_PARQUET
