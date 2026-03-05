// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <memory>

namespace OpenMS
{
  // Forward-declare ExperimentalSettings to avoid heavy includes in header
  class ExperimentalSettings;

  class MobilogramParquetConsumerImpl;

  /// Writes Mobilograms (ion-mobility traces) to a Parquet file using a PyProphet-compatible-like schema.
  class OPENMS_DLLAPI MobilogramParquetConsumer
  {
  public:
    /**
      @brief Construct a parquet consumer for mobilogram export.

      @param[in] filename Output parquet filename.
      @param[in] run_id Run identifier to store with each mobilogram.
      @param[in] source_file Source mzML filename to store with each mobilogram.
      @param[in] transition_exp Optional transition experiment used to annotate mobilograms (nullable use-case)
    */
    MobilogramParquetConsumer(const String& filename,
                             UInt64 run_id,
                             const String& source_file,
                             const OpenSwath::LightTargetedExperiment& transition_exp);

    /// Destructor flushes pending data and closes the parquet writer.
    ~MobilogramParquetConsumer();

    /// Consume a Mobilogram and append it to the parquet output.
    void consumeMobilogram(Mobilogram& m);

    /// Finalize and write the parquet file. Call to surface write errors.
    void finalize();

    /// Reserve storage for expected number of mobilograms
    void setExpectedSize(Size expectedMobilograms);

    /// Set experimental settings (currently unused)
    void setExperimentalSettings(const ExperimentalSettings& exp);

  private:
    std::unique_ptr<MobilogramParquetConsumerImpl> impl_;
  };

} // namespace OpenMS
