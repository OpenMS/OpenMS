// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportConfig.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathInferenceData.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/config.h>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace OpenSwath
{
struct LightTargetedExperiment;
}

namespace OpenMS
{
/**
  @brief Read and update OpenSWATH results in an OSWPQ directory or archive.

  Provides the result-storage operations used by OpenDIA without exposing
  Arrow types. The prepared library must use the same canonical numeric IDs
  as the library used to write the workflow. Its lookup data is owned by this
  object; the supplied experiment need not outlive the constructor.

  An archive is unpacked once into an owned temporary directory. Updates to
  that workspace are visible to subsequent reads, but only commit() repacks
  the archive. Destruction releases temporary files without committing. For
  directory inputs, writes take effect immediately and cannot be rolled back.
  Feature tables are processed one run at a time; only the run table is cached.

  @ingroup FileIO
*/
class OPENMS_DLLAPI OSWParquetFile : public ProgressLogger
{
public:
  /// Results for one inference level/context, including an explicitly empty result.
  struct OPENMS_DLLAPI InferenceResults
  {
    InferenceLevel level = InferenceLevel::Peptide;
    InferenceContext context = InferenceContext::Global;
    std::vector<LevelContextResultRow> rows;
  };

  /**
    @brief Open a workflow and build its prepared-library lookup.
    @param[in] filename OSWPQ directory or ZIP archive
    @param[in] library Prepared library with canonical numeric IDs
    @param[in] load_transition_metadata Retain transition metadata for feature/transition
               exports and transition quantification; may be false for inference-only use
  */
  OSWParquetFile(const std::string& filename, const OpenSwath::LightTargetedExperiment& library, bool load_transition_metadata = true);
  /// Release the workspace without implicitly committing archive changes.
  ~OSWParquetFile() override;
  OSWParquetFile(const OSWParquetFile&) = delete;
  OSWParquetFile& operator=(const OSWParquetFile&) = delete;
  /// Transfer ownership of the workspace.
  OSWParquetFile(OSWParquetFile&&);
  /// Transfer ownership, releasing the old workspace without committing it.
  OSWParquetFile& operator=(OSWParquetFile&&);

  /// Read compact peptide-, protein-, or gene-level rows for context inference.
  std::vector<LevelContextInputRow> readLevelContextData(InferenceLevel level, InferenceContext context) const;

  /**
    @brief Write inference tables and synchronize the corresponding feature score columns.
    @param[in] results Results for all requested level/context combinations
    All contexts for each supplied level replace that level's inference table.
    Only the requested feature score columns are replaced; an empty result
    writes null scores. Each run's feature table is updated once for the batch.
    Pass the batch with std::move to avoid copying it at the API boundary.
    Peptidoform inference is not supported. Changes to archives require commit().
  */
  void writeLevelContextResults(std::vector<InferenceResults> results);

  /// Read filtered rows for result and matrix exporters.
  std::vector<OpenSwathExportRow> readOpenSwathExportRows(const OpenSwathExportFilterConfig& config) const;
  /// Read feature score rows for the Parquet exporter.
  OpenSwathFeatureScoreTable readOpenSwathFeatureScoreTable(const OpenSwathParquetExportConfig& config) const;
  /// Read transition score rows for the Parquet exporter.
  OpenSwathTransitionScoreTable readOpenSwathTransitionScoreTable(const OpenSwathParquetExportConfig& config) const;
  /// Read run IDs and user-facing filename stems.
  std::map<Int64, std::string> readRunBasenames() const;

  /// Persist pending archive changes and refresh the sidecar index; a no-op if unchanged.
  void commit();

private:
  class Impl;
  std::unique_ptr<Impl> impl_;
};
} // namespace OpenMS
