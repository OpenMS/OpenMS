// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXValueValidation.h>

#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>

#include <arrow/api.h>

#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>

namespace OpenMS
{
  namespace
  {
    void addError(QPXValueValidation::Result& result, std::string error)
    {
      result.valid = false;
      result.errors.push_back(std::move(error));
    }

    std::shared_ptr<arrow::Array> combinedColumn(
      const std::shared_ptr<arrow::Table>& table,
      const std::string& name,
      QPXValueValidation::Result& result)
    {
      const auto column = table->GetColumnByName(name);
      if (!column)
      {
        addError(result, "missing required column '" + name + "'");
        return nullptr;
      }
      if (column->num_chunks() == 0)
      {
        addError(result, "column '" + name + "' has no chunks");
        return nullptr;
      }
      if (column->num_chunks() == 1) { return column->chunk(0); }

      auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
      if (!combined.ok())
      {
        addError(result, "could not combine chunks of column '" + name + "': "
                         + combined.status().ToString());
        return nullptr;
      }
      return combined.ValueOrDie();
    }

    void appendString(std::string& key, const std::string& value)
    {
      key += std::to_string(value.size());
      key += ':';
      key += value;
      key += '|';
    }

    void appendNullableString(std::string& key, const std::optional<std::string>& value)
    {
      if (!value)
      {
        key += "N|";
        return;
      }
      key += "S";
      appendString(key, *value);
    }

    template <typename Integer>
    void appendInteger(std::string& key, Integer value)
    {
      key += std::to_string(value);
      key += '|';
    }

    void appendFloat(std::string& key, float value)
    {
      // QPX/Arrow equality does not distinguish +0 from -0. Normalize before encoding the bits so
      // the in-memory key comparison has the same behaviour.
      if (value == 0.0f) { value = 0.0f; }
      appendInteger(key, std::bit_cast<std::uint32_t>(value));
    }

    /// Append a nullable float key component; a null is its own value, distinct from any float.
    void appendNullableFloat(
      std::string& key, const std::shared_ptr<arrow::FloatArray>& values, int64_t row)
    {
      if (values->IsNull(row)) { key += "N|"; return; }
      key += "F|";
      appendFloat(key, values->Value(row));
    }

    std::shared_ptr<arrow::Schema> schemaFor(QPXValueValidation::View view)
    {
      switch (view)
      {
        case QPXValueValidation::View::PSM:           return QPXPSMSchema::schema();
        case QPXValueValidation::View::FEATURE:       return QPXFeatureSchema::schema();
        case QPXValueValidation::View::PROTEIN_GROUP: return QPXPgSchema::schema();
      }
      return nullptr;
    }

    std::string viewName(QPXValueValidation::View view)
    {
      switch (view)
      {
        case QPXValueValidation::View::PSM:           return "psm";
        case QPXValueValidation::View::FEATURE:       return "feature";
        case QPXValueValidation::View::PROTEIN_GROUP: return "pg";
      }
      return "unknown";
    }

    void requiredColumnsHaveValues(
      const std::shared_ptr<arrow::Table>& table,
      const std::shared_ptr<arrow::Schema>& expected,
      QPXValueValidation::Result& result)
    {
      for (const auto& field : expected->fields())
      {
        if (field->nullable()) { continue; }
        const auto column = table->GetColumnByName(field->name());
        if (column && column->null_count() > 0)
        {
          addError(result, "required column '" + field->name() + "' contains "
                           + std::to_string(column->null_count()) + " physical null value(s)");
        }
      }
    }

    bool isBlank(const std::string& value)
    {
      return std::all_of(value.begin(), value.end(),
                         [](unsigned char c) { return std::isspace(c); });
    }

    bool nonEmptyString(
      const std::shared_ptr<arrow::StringArray>& values,
      int64_t row,
      const std::string& name,
      QPXValueValidation::Result& result)
    {
      if (values->IsNull(row))
      {
        addError(result, "row " + std::to_string(row) + " has a null '" + name
                         + "' primary-key value");
        return false;
      }
      if (isBlank(values->GetString(row)))
      {
        addError(result, "row " + std::to_string(row) + " has an empty or blank '" + name
                         + "' primary-key value");
        return false;
      }
      return true;
    }

    void validateAdditionalIntensities(
      const std::shared_ptr<arrow::Table>& table,
      const std::string& column_name,
      const std::shared_ptr<arrow::StringArray>& scalar_labels,
      QPXValueValidation::Result& result)
    {
      auto outer = std::static_pointer_cast<arrow::ListArray>(
        combinedColumn(table, column_name, result));
      if (!outer) { return; }

      const auto entries = std::static_pointer_cast<arrow::StructArray>(outer->values());
      const auto labels = std::static_pointer_cast<arrow::StringArray>(entries->field(0));
      const auto intensity_lists = std::static_pointer_cast<arrow::ListArray>(entries->field(1));
      const auto intensity_entries = std::static_pointer_cast<arrow::StructArray>(intensity_lists->values());
      const auto names = std::static_pointer_cast<arrow::StringArray>(intensity_entries->field(0));
      const auto values = std::static_pointer_cast<arrow::FloatArray>(intensity_entries->field(1));

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        if (outer->IsNull(row)) { continue; }
        std::unordered_set<std::string> row_labels;
        for (int64_t i = 0; i < outer->value_length(row); ++i)
        {
          const int64_t index = outer->value_offset(row) + i;
          if (entries->IsNull(index))
          {
            addError(result, "row " + std::to_string(row) + " contains a null "
                             + column_name + " entry");
            continue;
          }

          std::optional<std::string> label;
          if (labels->IsNull(index) || isBlank(labels->GetString(index)))
          {
            addError(result, "row " + std::to_string(row) + " contains an empty "
                             + column_name + "[].label");
          }
          else
          {
            label = labels->GetString(index);
            if (!ArrowIOHelpers::qpxIsCanonicalIntensityLabel(*label))
            {
              addError(result, "row " + std::to_string(row) + " has non-canonical "
                               + column_name + " label '" + *label + "'");
            }
            if (!row_labels.insert(*label).second)
            {
              addError(result, "row " + std::to_string(row) + " repeats " + column_name
                               + " label '" + *label + "'");
            }
            if (scalar_labels
                && (scalar_labels->IsNull(row) || scalar_labels->GetString(row) != *label))
            {
              addError(result, "row " + std::to_string(row) + " has " + column_name
                               + " label '" + *label + "' that differs from its pg label");
            }
          }

          if (intensity_lists->IsNull(index))
          {
            addError(result, "row " + std::to_string(row) + " has a null " + column_name
                             + "[].intensities list");
            continue;
          }
          std::unordered_set<std::string> row_names;
          for (int64_t j = 0; j < intensity_lists->value_length(index); ++j)
          {
            const int64_t value_index = intensity_lists->value_offset(index) + j;
            if (intensity_entries->IsNull(value_index))
            {
              addError(result, "row " + std::to_string(row) + " contains a null "
                               + column_name + "[].intensities entry");
              continue;
            }
            if (names->IsNull(value_index) || isBlank(names->GetString(value_index)))
            {
              addError(result, "row " + std::to_string(row) + " contains an empty "
                               + column_name + "[].intensities[].intensity_name");
            }
            else if (!row_names.insert(names->GetString(value_index)).second)
            {
              addError(result, "row " + std::to_string(row) + " repeats " + column_name
                               + " intensity name '" + names->GetString(value_index) + "'");
            }
            if (values->IsNull(value_index) || !std::isfinite(values->Value(value_index)))
            {
              addError(result, "row " + std::to_string(row) + " has a null or non-finite "
                               + column_name + " intensity value");
            }
          }
        }
      }
    }
  } // namespace

  struct QPXValueValidation::Impl
  {
    explicit Impl(View view) : view(view) {}

    View view;
    std::unordered_set<std::string> primary_keys;
    std::unordered_set<std::string> pg_run_keys;
  };

  std::string QPXValueValidation::Result::toString() const
  {
    std::ostringstream stream;
    for (size_t i = 0; i < errors.size(); ++i)
    {
      if (i != 0) { stream << "; "; }
      stream << errors[i];
    }
    return stream.str();
  }

  QPXValueValidation::QPXValueValidation(View view) : impl_(std::make_unique<Impl>(view)) {}

  QPXValueValidation::~QPXValueValidation() = default;
  QPXValueValidation::QPXValueValidation(QPXValueValidation&&) noexcept = default;
  QPXValueValidation& QPXValueValidation::operator=(QPXValueValidation&&) noexcept = default;

  void QPXValueValidation::reset()
  {
    impl_->primary_keys.clear();
    impl_->pg_run_keys.clear();
  }

  QPXValueValidation::Result QPXValueValidation::validate(
    const std::shared_ptr<arrow::Table>& table)
  {
    Result result;
    if (!table)
    {
      addError(result, viewName(impl_->view) + " table is null");
      return result;
    }

    const auto expected = schemaFor(impl_->view);
    const auto schema_validation = ArrowSchemaValidation::validate(
      table, expected, ArrowSchemaValidation::Mode::Strict);
    if (!schema_validation.valid)
    {
      for (const auto& error : schema_validation.errors)
      {
        addError(result, "schema: " + error);
      }
      return result;
    }

    const auto arrow_validation = table->ValidateFull();
    if (!arrow_validation.ok())
    {
      addError(result, "Arrow table is inconsistent: " + arrow_validation.ToString());
      return result;
    }
    requiredColumnsHaveValues(table, expected, result);

    std::unordered_set<std::string> new_primary_keys;
    std::unordered_set<std::string> new_pg_run_keys;

    if (impl_->view == View::PSM)
    {
      auto peptidoform = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXPSMSchema::PEPTIDOFORM, result));
      auto charge = std::static_pointer_cast<arrow::Int16Array>(
        combinedColumn(table, QPXPSMSchema::CHARGE, result));
      auto run = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXPSMSchema::RUN_FILE_NAME, result));
      auto scan = std::static_pointer_cast<arrow::ListArray>(
        combinedColumn(table, QPXPSMSchema::SCAN, result));
      if (!peptidoform || !charge || !run || !scan) { return result; }
      const auto scan_values = std::static_pointer_cast<arrow::Int32Array>(scan->values());

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        const bool peptidoform_valid = nonEmptyString(
          peptidoform, row, QPXPSMSchema::PEPTIDOFORM, result);
        const bool run_valid = nonEmptyString(
          run, row, QPXPSMSchema::RUN_FILE_NAME, result);
        bool charge_valid = true;
        if (charge->IsNull(row))
        {
          addError(result, "row " + std::to_string(row)
                           + " has a null 'charge' primary-key value");
          charge_valid = false;
        }
        bool key_valid = peptidoform_valid && run_valid && charge_valid;
        if (scan->IsNull(row) || scan->value_length(row) == 0)
        {
          addError(result, "row " + std::to_string(row)
                           + " has an empty 'scan' primary-key value");
          key_valid = false;
        }

        std::vector<int32_t> scans;
        if (!scan->IsNull(row))
        {
          scans.reserve(static_cast<size_t>(scan->value_length(row)));
          for (int64_t i = 0; i < scan->value_length(row); ++i)
          {
            const int64_t index = scan->value_offset(row) + i;
            if (scan_values->IsNull(index))
            {
              addError(result, "row " + std::to_string(row)
                               + " has a null scan component");
              key_valid = false;
            }
            else
            {
              scans.push_back(scan_values->Value(index));
            }
          }
        }
        if (!key_valid) { continue; }

        // QPX canonicalizes list-valued primary-key fields by sorting before comparison.
        std::sort(scans.begin(), scans.end());
        std::string key;
        appendString(key, peptidoform->GetString(row));
        appendInteger(key, charge->Value(row));
        appendString(key, run->GetString(row));
        appendInteger(key, scans.size());
        for (const int32_t value : scans) { appendInteger(key, value); }

        if (impl_->primary_keys.contains(key) || !new_primary_keys.insert(key).second)
        {
          addError(result, "row " + std::to_string(row)
                           + " repeats the QPX psm primary key "
                             "(peptidoform, charge, run_file_name, scan)");
        }
      }
    }
    else if (impl_->view == View::FEATURE)
    {
      auto peptidoform = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXFeatureSchema::PEPTIDOFORM, result));
      auto charge = std::static_pointer_cast<arrow::Int16Array>(
        combinedColumn(table, QPXFeatureSchema::CHARGE, result));
      auto run = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXFeatureSchema::RUN_FILE_NAME, result));
      auto rt = std::static_pointer_cast<arrow::FloatArray>(
        combinedColumn(table, QPXFeatureSchema::RT, result));
      auto observed_mz = std::static_pointer_cast<arrow::FloatArray>(
        combinedColumn(table, QPXFeatureSchema::OBSERVED_MZ, result));
      auto intensities = std::static_pointer_cast<arrow::ListArray>(
        combinedColumn(table, QPXFeatureSchema::INTENSITIES, result));
      if (!peptidoform || !charge || !run || !rt || !observed_mz || !intensities) { return result; }

      const auto intensity_values = std::static_pointer_cast<arrow::StructArray>(intensities->values());
      const auto labels = std::static_pointer_cast<arrow::StringArray>(intensity_values->field(0));
      const auto values = std::static_pointer_cast<arrow::FloatArray>(intensity_values->field(1));

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        // QPX explicitly permits unmapped features; OpenMS represents their unknown peptidoform
        // as an empty string. The run identity, unlike the optional mapping, must always be
        // present.
        //
        // 'observed_mz' is part of the key precisely because of those unmapped rows. With an
        // empty peptidoform the key collapses to (charge, run_file_name, rt), and 'rt' is
        // narrowed to float32 on write - one ULP is 244 us at a retention time of 3000 s. Two
        // co-eluting unmapped features of the same charge in one run then share a key even
        // though they are plainly different features, and the whole export is refused. That is
        // not a corner case for the producers OpenMS ships: ProteomicsLFQ with seeds leaves the
        // large majority of its consensus features unidentified (79 % in the
        // ProteomicsLFQ_3 test fixture), and a real run carries tens of thousands of them, where
        // the birthday bound over ~2.2e7 representable float32 retention times makes a collision
        // the expected outcome rather than a remote one. The m/z is written one column over and
        // separates them.
        const bool strings_valid = nonEmptyString(
          run, row, QPXFeatureSchema::RUN_FILE_NAME, result);
        const bool peptidoform_valid = !peptidoform->IsNull(row);
        if (!peptidoform_valid)
        {
          addError(result, "row " + std::to_string(row)
                           + " has a null 'peptidoform' primary-key value");
        }
        const bool charge_valid = !charge->IsNull(row);
        if (!charge_valid)
        {
          addError(result, "row " + std::to_string(row)
                           + " has a null 'charge' primary-key value");
        }
        bool rt_valid = true;
        if (!rt->IsNull(row) && !std::isfinite(rt->Value(row)))
        {
          addError(result, "row " + std::to_string(row)
                           + " has a non-finite 'rt' primary-key value");
          rt_valid = false;
        }
        bool observed_mz_valid = true;
        if (!observed_mz->IsNull(row) && !std::isfinite(observed_mz->Value(row)))
        {
          addError(result, "row " + std::to_string(row)
                           + " has a non-finite 'observed_mz' primary-key value");
          observed_mz_valid = false;
        }
        if (strings_valid && peptidoform_valid && charge_valid && rt_valid && observed_mz_valid)
        {
          std::string key;
          appendString(key, peptidoform->GetString(row));
          appendInteger(key, charge->Value(row));
          appendString(key, run->GetString(row));
          appendNullableFloat(key, rt, row);
          appendNullableFloat(key, observed_mz, row);
          if (impl_->primary_keys.contains(key) || !new_primary_keys.insert(key).second)
          {
            addError(result, "row " + std::to_string(row)
                             + " repeats the QPX feature primary key "
                               "(peptidoform, charge, run_file_name, rt, observed_mz)");
          }
        }

        if (intensities->IsNull(row)) { continue; }
        std::unordered_set<std::string> row_labels;
        for (int64_t i = 0; i < intensities->value_length(row); ++i)
        {
          const int64_t index = intensities->value_offset(row) + i;
          if (intensity_values->IsNull(index))
          {
            addError(result, "row " + std::to_string(row)
                             + " contains a null intensities entry");
            continue;
          }
          if (labels->IsNull(index) || isBlank(labels->GetString(index)))
          {
            addError(result, "row " + std::to_string(row)
                             + " contains an empty intensities[].label");
            continue;
          }
          const std::string label = labels->GetString(index);
          if (!ArrowIOHelpers::qpxIsCanonicalIntensityLabel(label))
          {
            addError(result, "row " + std::to_string(row) + " has non-canonical intensity label '"
                             + label + "'");
          }
          if (!row_labels.insert(label).second)
          {
            addError(result, "row " + std::to_string(row) + " repeats intensity label '"
                             + label + "'");
          }
          if (values->IsNull(index) || !std::isfinite(values->Value(index)))
          {
            addError(result, "row " + std::to_string(row) + " has a null or non-finite intensity for label '"
                             + label + "'");
          }
        }
      }
      validateAdditionalIntensities(
        table, QPXFeatureSchema::ADDITIONAL_INTENSITIES, nullptr, result);
    }
    else
    {
      auto anchor = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXPgSchema::ANCHOR_PROTEIN, result));
      auto grouped_runs = std::static_pointer_cast<arrow::ListArray>(
        combinedColumn(table, QPXPgSchema::GROUPED_RUNS, result));
      auto label = std::static_pointer_cast<arrow::StringArray>(
        combinedColumn(table, QPXPgSchema::LABEL, result));
      auto intensity = std::static_pointer_cast<arrow::FloatArray>(
        combinedColumn(table, QPXPgSchema::INTENSITY, result));
      if (!anchor || !grouped_runs || !label || !intensity) { return result; }
      const auto run_values = std::static_pointer_cast<arrow::StringArray>(grouped_runs->values());

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        bool key_valid = nonEmptyString(anchor, row, QPXPgSchema::ANCHOR_PROTEIN, result);
        std::vector<std::string> runs;
        std::unordered_set<std::string> row_runs;
        if (grouped_runs->IsNull(row) || grouped_runs->value_length(row) == 0)
        {
          addError(result, "row " + std::to_string(row)
                           + " has an empty 'grouped_runs' primary-key value");
          key_valid = false;
        }
        else
        {
          runs.reserve(static_cast<size_t>(grouped_runs->value_length(row)));
          for (int64_t i = 0; i < grouped_runs->value_length(row); ++i)
          {
            const int64_t index = grouped_runs->value_offset(row) + i;
            if (run_values->IsNull(index) || isBlank(run_values->GetString(index)))
            {
              addError(result, "row " + std::to_string(row)
                               + " has an empty member in grouped_runs");
              key_valid = false;
              continue;
            }
            const std::string run_name = run_values->GetString(index);
            if (!row_runs.insert(run_name).second)
            {
              addError(result, "row " + std::to_string(row) + " repeats run '" + run_name
                               + "' inside grouped_runs");
              key_valid = false;
              continue;
            }
            runs.push_back(run_name);
          }
        }

        const bool has_label = !label->IsNull(row);
        const bool has_intensity = !intensity->IsNull(row);
        if (has_label != has_intensity)
        {
          addError(result, "row " + std::to_string(row)
                           + " must set label and intensity together, or leave both null");
        }
        std::optional<std::string> label_value;
        if (has_label)
        {
          label_value = label->GetString(row);
          if (isBlank(*label_value))
          {
            addError(result, "row " + std::to_string(row)
                             + " has an empty 'label' primary-key value");
            key_valid = false;
          }
          else if (!ArrowIOHelpers::qpxIsCanonicalIntensityLabel(*label_value))
          {
            addError(result, "row " + std::to_string(row) + " has non-canonical label '"
                             + *label_value + "'");
          }
        }
        if (has_intensity && !std::isfinite(intensity->Value(row)))
        {
          addError(result, "row " + std::to_string(row) + " has a non-finite intensity");
        }

        std::sort(runs.begin(), runs.end());
        if (!key_valid) { continue; }

        std::string key;
        appendString(key, anchor->GetString(row));
        appendInteger(key, runs.size());
        for (const auto& run_name : runs) { appendString(key, run_name); }
        appendNullableString(key, label_value);
        if (impl_->primary_keys.contains(key) || !new_primary_keys.insert(key).second)
        {
          addError(result, "row " + std::to_string(row)
                           + " repeats the QPX pg primary key "
                             "(anchor_protein, grouped_runs, label)");
        }

        // One raw file may contribute to at most one protein quantity at this
        // (anchor_protein, label). Otherwise the same measurement is double-counted.
        for (const auto& run_name : runs)
        {
          std::string run_key;
          appendString(run_key, anchor->GetString(row));
          appendNullableString(run_key, label_value);
          appendString(run_key, run_name);
          if (impl_->pg_run_keys.contains(run_key) || !new_pg_run_keys.insert(run_key).second)
          {
            addError(result, "row " + std::to_string(row) + " reuses run '" + run_name
                             + "' in another pg row with the same (anchor_protein, label)");
          }
        }
      }
      validateAdditionalIntensities(
        table, QPXPgSchema::ADDITIONAL_INTENSITIES, label, result);
    }

    if (result.valid)
    {
      impl_->primary_keys.insert(new_primary_keys.begin(), new_primary_keys.end());
      impl_->pg_run_keys.insert(new_pg_run_keys.begin(), new_pg_run_keys.end());
    }
    return result;
  }
} // namespace OpenMS
