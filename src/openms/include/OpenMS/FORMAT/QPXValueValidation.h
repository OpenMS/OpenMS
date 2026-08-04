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
