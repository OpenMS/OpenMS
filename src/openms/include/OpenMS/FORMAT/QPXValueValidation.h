// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <memory>
#include <string>
#include <vector>

namespace arrow
{
  class Table;
}

namespace OpenMS
{
  /**
    @brief Stateful value and primary-key validator for QPX Arrow tables

    Arrow schema validation checks declarations such as field type and nullability. It does not
    reject an empty string in an identity-bearing key, physical nulls under a non-nullable field,
    duplicate primary keys, or invalid nested values. This class supplies those write-time checks
    for the OpenMS QPX psm, feature, and pg views. QPX-permitted unmapped features retain OpenMS'
    empty peptidoform representation and are still distinguished by their remaining key fields.

    A validator retains primary keys between validate() calls. Streaming writers can therefore
    validate consecutive batches with one instance and detect duplicates that cross a batch
    boundary. Use reset() before starting a different output file.

    <b>The contract this class enforces</b>

    This is the single place the QPX input contract is written down; the tools that offer
    @c -out_qpx (ProteinQuantifier, ProteomicsLFQ, IsobaricWorkflow, ProSE) point here rather than
    restating it, so it cannot drift from the code. A table that violates it is refused outright:
    the export writes nothing, rather than emitting rows that will not join.

    Consensus maps produced by @c ProteomicsLFQ and @c IsobaricWorkflow satisfy the contract by
    construction, and they are the supported producers. A map assembled by a different pipeline
    may not, and the requirements it can fail are:

    - <b>Every identification needs a scan.</b> The @c scan primary-key component is derived from
      the identification's spectrum reference. Identifications carrying none - protein-inference
      output, some ID converters, transferred (match-between-runs) identifications - cannot be
      exported.
    - <b>Every identification needs a resolvable origin file.</b> @c run_file_name is the key the
      psm view is joined to the feature and pg views on. Protein inference commonly drops the
      path; restore it with ProteinIdentification::setPrimaryMSRunPath() before exporting.
    - <b>Identifications must be unique per (peptidoform, charge, run, scan).</b> IDMapper maps one
      MS2 identification onto every consensus feature whose RT/m-z window contains it, so the same
      identification can end up on two features. Resolve conflicts first, e.g. with
      IDConflictResolver.
    - <b>Every channel needs a QPX label.</b> The vocabulary is label-free, the TMT and iTRAQ
      plexes, and SILAC 2-/3-plex. Those labels are join keys against the run sidecar, so OpenMS
      refuses rather than guesses one. Maps from FeatureFinderMultiplex carrying @c ICPL*,
      @c Leu3, @c "label N" or @c no_label fall outside the vocabulary and cannot be exported.

    Feature rows are keyed on (peptidoform, charge, run_file_name, rt, observed_mz). @c observed_mz
    is part of that key for the sake of unmapped features: QPX permits them, OpenMS writes them
    with an empty peptidoform, and @c rt is narrowed to float32 on write - so without the m/z two
    co-eluting unmapped features of one charge in one run would share a key and the export would
    be refused. Rows agreeing on all five are genuinely indistinguishable and are still refused.

    @experimental This API is experimental and may change in future versions.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI QPXValueValidation
  {
  public:
    /// QPX table whose value contract is checked.
    enum class View
    {
      PSM,
      FEATURE,
      PROTEIN_GROUP
    };

    /// Result of validating one table or streaming batch.
    struct OPENMS_DLLAPI Result
    {
      bool valid = true;
      std::vector<std::string> errors;

      /**
        @brief Join all diagnostics into one log-friendly string

        @return Diagnostics separated by semicolons
      */
      std::string toString() const;
    };

    /**
      @brief Construct a validator for one QPX view

      @param[in] view View whose schema, values, and primary key are checked
    */
    explicit QPXValueValidation(View view);
    ~QPXValueValidation();

    QPXValueValidation(QPXValueValidation&&) noexcept;
    QPXValueValidation& operator=(QPXValueValidation&&) noexcept;

    QPXValueValidation(const QPXValueValidation&) = delete;
    QPXValueValidation& operator=(const QPXValueValidation&) = delete;

    /**
      @brief Validate one whole table or the next streaming batch

      Checks the exact OpenMS QPX schema, physical nulls in required columns, non-empty
      identity-bearing primary-key values, primary-key uniqueness, canonical quantitative labels,
      and view-specific nested invariants. For pg data this includes non-empty, duplicate-free
      @c grouped_runs and run-disjointness within one @c (anchor_protein, label).

      State is advanced only when the complete table is valid, so a rejected batch can be fixed and
      submitted again without first calling reset().

      @param[in] table Table or streaming batch to validate
      @return Structured validation result; @c valid is false if writing must be refused
    */
    Result validate(const std::shared_ptr<arrow::Table>& table);

    /// Forget primary keys accumulated from earlier validate() calls.
    void reset();

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };
} // namespace OpenMS
