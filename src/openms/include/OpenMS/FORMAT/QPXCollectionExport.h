// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS
{
  class ConsensusMap;
  class ExperimentalDesign;

/**
  @brief Whole-collection preflight for the QPX Parquet export (@c -out_qpx)

  A QPX collection is three files written in a fixed order -- @c quantms.feature.parquet, then
  @c quantms.psm.parquet, then @c quantms.pg.parquet. Each view refuses input it cannot
  represent, but it does so before opening <i>its own</i> file, so a refusal in the second or
  third view leaves the earlier files on disk: a partial collection that looks like a
  successful export.

  The exporters cannot avoid that themselves. Both streaming paths open the output and emit a
  valid 0-row file even for empty input, so by the time a view can refuse, its predecessors
  have already written. Only a check that runs before the first write can leave the output
  directory untouched.

  requireExportable() therefore runs every deterministic refusal of all three views up front.
  It does not replace the in-exporter checks -- those remain the contract for callers using an
  exporter directly -- so a condition is validated twice on the TOPP path. That is deliberate:
  the preflight must not become the only thing standing between a library user and a malformed
  table.

  @note Only conditions that are decidable from the input are covered. I/O failures, and
        anything that can only be discovered while building a row, cannot be. Those are handled
        from the other side, by Transaction: it removes whatever was written when the export does
        not run to completion. Use both -- the preflight so a refusable input never opens a file,
        the transaction so anything that slips past it is not left behind.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI QPXCollectionExport
{
public:
  /**
    @brief Refuse a collection that cannot be written in full, before any file is created

    Performs, on the same input the three views will see:
      - the feature view's divergent-identification and merge-index refusals,
      - the psm view's merge-index refusal over assigned <i>and</i> unassigned identifications,
      - the channel-label check shared by the feature and pg views,
      - the pg view's protein/group prerequisites and fraction-group checks (invalid sample,
        inconsistent cell, ragged group, repeated run, or ambiguous/cross-run label),
      - every quantified group's complete parallel @c (fraction_group, label, abundance)
        annotation contract and its agreement with the design.

    @param[in] cmap The ConsensusMap about to be exported
    @param[in] design The experimental design that was used to quantify @p cmap
    @return true if all three views can be written; false for the conditions the views report
            by returning false/nullptr, having logged the reason
    @throws Exception::IllegalArgument if a feature carries divergent peptide annotations
    @throws Exception::MissingInformation if an identification of a merged run carries no usable
           @c id_merge_index
    @throws Exception::InvalidValue if two identification runs share an identifier
  */
  static bool requireExportable(const ConsensusMap& cmap, const ExperimentalDesign& design);

  /**
    @brief Scope guard that leaves no partial QPX collection behind

    requireExportable() cannot decide everything: a view can still refuse a row it only sees
    while building the table, and a write can fail on I/O. Either way the earlier files of the
    collection are already on disk, and a directory holding @c quantms.feature.parquet but no
    @c quantms.psm.parquet is indistinguishable from a complete export of a smaller dataset.

    Construct this before the first write and call commit() once the last one has succeeded.
    On any other exit -- a view returning false, a writer throwing, an early return -- the
    destructor removes every collection file present in the directory, so the export is
    all-or-nothing. Enumerating the failure sites instead would silently miss the next one.

    @code
    QPXCollectionExport::Transaction qpx(out_qpx);
    if (!ConsensusMapArrowExport::exportToParquet(cmap, out_qpx + "/quantms.feature.parquet")) { throw ...; }
    if (!QPXFile::exportToParquet(prot_ids, pep_ids, out_qpx + "/quantms.psm.parquet"))        { throw ...; }
    if (!ProteinGroupArrowExport::exportToParquet(cmap, design, out_qpx + "/quantms.pg.parquet")) { throw ...; }
    qpx.commit();
    @endcode

    @experimental This API is experimental and may change in future versions.
  */
  class OPENMS_DLLAPI Transaction
  {
  public:
    /// Start guarding the collection directory @p directory
    explicit Transaction(const std::string& directory);

    /// Removes the collection files unless commit() was called. Never throws.
    ~Transaction();

    /// Keep the collection: the export completed
    void commit();

    Transaction(const Transaction&) = delete;
    Transaction& operator=(const Transaction&) = delete;
    Transaction(Transaction&&) = delete;
    Transaction& operator=(Transaction&&) = delete;

  private:
    std::string directory_;
    bool committed_ {false};
  };
};

} // namespace OpenMS
