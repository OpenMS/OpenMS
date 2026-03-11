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

#include <limits>
#include <memory>

namespace OpenMS
{
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
      @param[in] transition_exp Optional transition experiment used to annotate mobilograms. Note: the
                                signature takes a const reference, so a caller that does not want any
                                transition annotations must pass an empty `OpenSwath::LightTargetedExperiment`
                                (i.e. one with no compounds/transitions). A null pointer is not accepted.
    */
    MobilogramParquetConsumer(const String& filename,
                             UInt64 run_id,
                             const String& source_file,
                             const OpenSwath::LightTargetedExperiment& transition_exp);

    /// Destructor flushes pending data and closes the parquet writer.
    ~MobilogramParquetConsumer();

    /**
      @brief Consume a mobilogram and write it to the parquet file.

      @param[in] m The mobilogram to consume (read-only; the caller retains ownership and the mobilogram is not modified)
      @param[in] mobilogram_type The type of the mobilogram (e.g. "precursor" or "fragment")
      @param[in] ms_level The MS level of the mobilogram (e.g. 1 for precursor, 2 for fragment)
      @param[in] transition_id The id of the corresponding transition in the transition experiment (nullable: -1 means not set)
      @param[in] transition_native_id The native id of the corresponding transition in the transition experiment (nullable: empty string means not set)
      @param[in] feature_rt Optional retention time apex that the mobilogram corresponds to (nullable: NaN means not set)
      @param[in] feature_id Optional feature id associated with the mobilogram (nullable: -1 means not set)
    */
    void consumeMobilogram(const Mobilogram& m,
                 const String& mobilogram_type = "",
                 Int64 ms_level = -1,
                 Int64 transition_id = -1,
                 const String& transition_native_id = "",
                 double feature_rt = std::numeric_limits<double>::quiet_NaN(),
                 Int64 feature_id = -1);

    /// Finalize and write the parquet file. Call to surface write errors.
    void finalize();

    /**
      @brief Reserve storage for expected number of mobilograms.

      Must be called before any concurrent consumeMobilogram() calls (e.g. from
      an OpenMP parallel loop). Calling it after concurrent writes have started
      is undefined behaviour.

      @param[in] expectedMobilograms The expected number of mobilograms
    */
    void setExpectedSize(Size expectedMobilograms);

  private:
    std::unique_ptr<MobilogramParquetConsumerImpl> impl_;
  };

} // namespace OpenMS
