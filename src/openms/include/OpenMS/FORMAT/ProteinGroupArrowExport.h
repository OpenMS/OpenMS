// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <memory>
#include <string>
#include <vector>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{
  class ExperimentalDesign;
  class ProteinIdentification;
}

namespace OpenMS
{

/**
  @brief Export protein group data to Apache Arrow format following QPX pg schema

  This class provides static methods to export protein group quantification data
  from a ConsensusMap to Apache Arrow Tables and Parquet files. The schema follows
  the QPX (Quantitative Proteomics Exchange) protein group format.

  Protein groups must have quantification annotated via
  PeptideAndProteinQuant::annotateQuantificationsToProteins() before export.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ProteinGroupArrowExport
{
public:
  /**
    @brief Export protein group data to Apache Arrow Table

    Exports indistinguishable protein groups following the active QPX 1.1 pg schema. One row is
    emitted per protein group, experimental-design @c fraction_group, and label. @c grouped_runs
    lists that fraction group's raw files, while scalar @c label and @c intensity carry one
    quantity. Together with @c anchor_protein these columns form the primary key.

    A fraction group must be <b>rectangular</b>: every label it publishes has to exist in every one
    of its runs, which makes "(any file in grouped_runs, label) -> run.samples[]" resolvable.
    Designs that are ragged, assign one run to several fraction groups, or disagree about a label
    are refused. Reusing one sample in several fraction groups is valid and does not join them.

    PeptideAndProteinQuant calculates protein quantities at @c (fraction_group, label) grain before
    protein aggregation and annotates them in named parallel data arrays. The exporter refuses a
    legacy sample-only annotation because an experiment-wide sample abundance cannot be split back
    into independent fraction-group quantities.

    Groups without quantification are melted onto the fraction groups in which a member accession
    has peptide evidence, with null @c label and @c intensity.

    @param[in] cmap The ConsensusMap with annotated protein group quantification
    @param[in] design The experimental design that was used to quantify @p cmap. It defines the
               exact fraction-group/label keys and their grouped runs. The exporter checks those
               keys against every quantified protein group.
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap,
                                                     const ExperimentalDesign& design);

  /**
    @brief Export protein group data to Apache Arrow Table, deriving the design from @p cmap

    Convenience overload for callers that do not hold the design used for quantification. It
    reconstructs one with ExperimentalDesign::fromConsensusMap(), which recovers the fraction
    grouping from the column headers. Prefer the overload taking an explicit design whenever the
    caller has one — a reconstructed design only reproduces the original sample numbering when
    quantification was itself driven by a reconstructed design.

    @param[in] cmap The ConsensusMap with annotated protein group quantification
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap);

  /**
    @brief Export protein group data to Parquet file

    @param[in] cmap The ConsensusMap with annotated protein group quantification
    @param[in] design The experimental design used to quantify @p cmap (see exportToArrow)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const ExperimentalDesign& design,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Export protein group data to Parquet file, deriving the design from @p cmap

    @param[in] cmap The ConsensusMap with annotated protein group quantification
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Export protein group data to Arrow table from identification data (no quantification)

    For search-engine output where no ConsensusMap is available. Populates required
    QPX pg fields (pg_accessions, anchor_protein, grouped_runs, is_decoy, peptides)
    and sets quantification columns (label, intensity, additional_intensities) to null.

    This emits <b>one row per protein group per run</b>: the runs are those in which a member
    accession has peptide evidence, resolved per PSM through IdentifierMSRunMapper. A group in a
    single-file run is keyed on that file even without evidence; a group in a merged run with no
    evidence anywhere is skipped with a diagnostic, since @c grouped_runs is a QPX primary-key
    component that must not be empty.

    Identification input carries no experimental design, so there is nothing to aggregate and
    every @c grouped_runs list has exactly one element — the single-element form the QPX 1.1
    schema documents for unfractionated input.

    @c peptides and @c peptide_counts are scoped to the emitted run, so the same group in two
    runs reports the peptides seen in each.

    @param[in] protein_identifications Protein identifications with protein groups
    @param[in] peptide_identifications Peptide identifications (for peptide-per-protein counts)
    @return Shared pointer to Arrow Table (empty table if no groups), or nullptr if an Arrow
            builder fails
    @throws Exception::MissingInformation if a PSM belongs to a merged run but carries no usable
           @c id_merge_index, so its origin file -- and with it the row's key -- is undetermined
  */
  static std::shared_ptr<arrow::Table> exportToArrow(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications);

  /**
    @brief Export protein group data to Parquet file from identification data (no quantification)

    @param[in] protein_identifications Protein identifications with protein groups
    @param[in] peptide_identifications Peptide identifications (for peptide-per-protein counts)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Write a pre-built QPX protein group Arrow table to a Parquet file

    The table is expected to follow QPXPgSchema (e.g., from exportToArrow).
    Attaches QPX file metadata (qpx_version, file_type="pg_file", UUID, creation_date)
    before writing. Use this overload when the caller already has the table built
    (e.g., for merged output) to avoid rebuilding it.

    @param[in] table QPX pg Arrow table (must not be null)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});
};

} // namespace OpenMS
