// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParquetTableComparator.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace OpenMS
{
  namespace
  {
    /// Read a Parquet file into a single-chunk Arrow table. Returns nullptr and logs on failure.
    std::shared_ptr<arrow::Table> readTable_(const std::string& filename)
    {
      auto infile_result = arrow::io::ReadableFile::Open(filename);
      if (!infile_result.ok())
      {
        OPENMS_LOG_ERROR << "ParquetTableComparator: cannot open '" << filename
                         << "': " << infile_result.status().ToString() << std::endl;
        return nullptr;
      }

      auto reader_result = parquet::arrow::OpenFile(*infile_result, arrow::default_memory_pool());
      if (!reader_result.ok())
      {
        OPENMS_LOG_ERROR << "ParquetTableComparator: '" << filename
                         << "' is not readable as Parquet: " << reader_result.status().ToString() << std::endl;
        return nullptr;
      }
      auto reader = std::move(reader_result.ValueOrDie());

      // The Result-returning overload is newer than the Arrow in contrib, so use the
      // pointer-out form that ConsensusMapArrowIO uses. It is deprecated in recent Arrow,
      // which is why -Wdeprecated-declarations fires here on newer system installs.
      std::shared_ptr<arrow::Table> table;
      const auto read_status = reader->ReadTable(&table);
      if (!read_status.ok())
      {
        OPENMS_LOG_ERROR << "ParquetTableComparator: cannot read table from '" << filename
                         << "': " << read_status.ToString() << std::endl;
        return nullptr;
      }

      // Deliberately NOT CombineChunks(): it documents that "binary columns may be combined
      // into multiple chunks" to avoid offset overflow, so a single chunk is not guaranteed
      // for large string columns. cellAt_ indexes the ChunkedArray globally instead.
      return table;
    }

    /// Value of a key in the table's schema metadata, or "" when absent.
    std::string metadataValue_(const std::shared_ptr<arrow::Table>& table, const std::string& key)
    {
      const auto& md = table->schema()->metadata();
      if (!md) { return ""; }
      const int idx = md->FindKey(key);
      return (idx < 0) ? std::string() : md->value(idx);
    }

    /**
      @brief Element at @p row of a column.

      Indexes the ChunkedArray globally rather than assuming a single chunk: a large string
      column can legitimately span several chunks. @p ok distinguishes "this cell is null"
      from "this cell could not be read" -- conflating the two would let an unreadable cell
      on both sides compare equal.
    */
    std::shared_ptr<arrow::Scalar> cellAt_(const std::shared_ptr<arrow::Table>& table,
                                           int col, int64_t row, bool& ok)
    {
      ok = false;
      if (col < 0) { return nullptr; }
      auto scalar_result = table->column(col)->GetScalar(row);
      if (!scalar_result.ok()) { return nullptr; }
      ok = true;
      return *scalar_result;
    }

    /// True for the Arrow scalar types that derive from BaseListScalar.
    bool isListLike_(arrow::Type::type id)
    {
      return id == arrow::Type::LIST || id == arrow::Type::LARGE_LIST
             || id == arrow::Type::FIXED_SIZE_LIST || id == arrow::Type::MAP;
    }

    /**
      @brief Canonical, self-delimiting encoding of a scalar, used to form primary keys.

      Every part is length-prefixed, so concatenating the encodings of several key columns is
      unambiguous: a value containing any byte -- including a separator character or the text
      of the null marker -- cannot be confused with a structural boundary. A plain delimiter
      would not be safe, because Arrow string values may contain arbitrary bytes.

      List elements are sorted so that a set-valued key column (QPX @c grouped_runs, and the
      @c scan component of the psm key) compares equal regardless of the order the producer
      happened to emit. This mirrors what the QPX reference implementation does before it
      checks primary-key uniqueness.
    */
    std::string canonicalKey_(const std::shared_ptr<arrow::Scalar>& scalar)
    {
      if (!scalar) { return "E;"; }        // unreadable
      if (!scalar->is_valid) { return "N;"; } // null

      if (isListLike_(scalar->type->id()))
      {
        const auto& values = static_cast<const arrow::BaseListScalar&>(*scalar).value;
        std::vector<std::string> parts;
        parts.reserve(static_cast<size_t>(values->length()));
        for (int64_t i = 0; i < values->length(); ++i)
        {
          auto element = values->GetScalar(i);
          parts.push_back(element.ok() ? canonicalKey_(*element) : "E;");
        }
        std::sort(parts.begin(), parts.end());
        std::string joined = "L" + std::to_string(parts.size()) + ";";
        for (const auto& p : parts) { joined += p; }
        return joined;
      }

      const std::string text = scalar->ToString();
      return "V" + std::to_string(text.size()) + ":" + text;
    }

    /// True when both are NaN, which regression comparison should treat as equal.
    bool bothNaN_(double a, double b)
    {
      return std::isnan(a) && std::isnan(b);
    }

    /**
      @brief Tolerant numeric comparison.

      Mirrors FuzzyStringComparator: either the absolute difference or the ratio has to be
      acceptable. A ratio is only meaningful for two non-zero values of the same sign; the
      "zero vs. epsilon" case is what @p acceptable_absdiff exists for.
    */
    bool numbersAcceptable_(double a, double b, const ParquetDiffSettings& settings,
                            double& ratio, double& absdiff)
    {
      ratio = 1.0;
      absdiff = 0.0;
      if (bothNaN_(a, b)) { return true; }
      if (std::isnan(a) || std::isnan(b)) { return false; }
      if (a == b) { return true; }

      absdiff = std::fabs(a - b);

      // Both metrics are computed unconditionally, even when the first one already accepts,
      // so that ParquetDiffResult::max_ratio / max_absdiff really are the largest observed.
      if (a != 0.0 && b != 0.0 && ((a > 0.0) == (b > 0.0)))
      {
        const double aa = std::fabs(a);
        const double bb = std::fabs(b);
        ratio = (aa > bb) ? (aa / bb) : (bb / aa);
      }
      else
      {
        // No ratio bridges zero to non-zero, or opposite signs; only absdiff can accept those.
        ratio = std::numeric_limits<double>::infinity();
      }
      return absdiff <= settings.acceptable_absdiff || ratio <= settings.acceptable_ratio;
    }

    /**
      @brief Numeric value of a scalar, if it holds one.

      Integers are included: "numeric" in the tolerance contract means every number, not only
      the floating-point ones. Values wider than the 53-bit mantissa of a double lose precision
      here, which the caller compensates for by never accepting a zero difference between two
      scalars that are not exactly equal.
    */
    bool asDouble_(const arrow::Scalar& scalar, double& out)
    {
      const auto id = scalar.type->id();
      if (!arrow::is_integer(id) && !arrow::is_floating(id)) { return false; }
      auto casted = scalar.CastTo(arrow::float64());
      if (!casted.ok()) { return false; }
      out = static_cast<const arrow::DoubleScalar&>(**casted).value;
      return true;
    }

    /**
      @brief Recursive, tolerance-aware scalar comparison.

      Descends into list and struct values so that a float nested inside
      @c intensities: list<struct{label, intensity}> is compared with the same tolerance as a
      top-level float. @p path accumulates a human-readable location for the message.
    */
    bool scalarsEqual_(const std::shared_ptr<arrow::Scalar>& a,
                       const std::shared_ptr<arrow::Scalar>& b,
                       const ParquetDiffSettings& settings,
                       const std::string& path,
                       std::string& message,
                       double& max_ratio,
                       double& max_absdiff)
    {
      const bool a_valid = a && a->is_valid;
      const bool b_valid = b && b->is_valid;
      if (!a_valid && !b_valid) { return true; }
      if (a_valid != b_valid)
      {
        message = path + ": " + (a_valid ? a->ToString() : "<null>") + " vs. "
                  + (b_valid ? b->ToString() : "<null>");
        return false;
      }

      if (!a->type->Equals(*b->type))
      {
        message = path + ": type " + a->type->ToString() + " vs. " + b->type->ToString();
        return false;
      }

      double av = 0.0;
      double bv = 0.0;
      if (asDouble_(*a, av) && asDouble_(*b, bv))
      {
        double ratio = 1.0;
        double absdiff = 0.0;
        bool ok = numbersAcceptable_(av, bv, settings, ratio, absdiff);
        // A 64-bit integer difference can vanish when both sides are widened to a double.
        // Two scalars that are not exactly equal must never be accepted on a zero difference.
        if (ok && absdiff == 0.0 && !a->Equals(*b)) { ok = false; }
        if (std::isfinite(ratio)) { max_ratio = std::max(max_ratio, ratio); }
        max_absdiff = std::max(max_absdiff, absdiff);
        if (!ok)
        {
          message = path + ": " + a->ToString() + " vs. " + b->ToString()
                    + " (absdiff " + std::to_string(absdiff) + ", ratio "
                    + (std::isfinite(ratio) ? std::to_string(ratio) : std::string("inf")) + ")";
        }
        return ok;
      }

      const auto id = a->type->id();
      if (isListLike_(id))
      {
        const auto& av_list = static_cast<const arrow::BaseListScalar&>(*a).value;
        const auto& bv_list = static_cast<const arrow::BaseListScalar&>(*b).value;
        if (av_list->length() != bv_list->length())
        {
          message = path + ": list length " + std::to_string(av_list->length()) + " vs. "
                    + std::to_string(bv_list->length());
          return false;
        }

        // Multiset semantics: pair elements by canonical string before comparing, so that a
        // reordered list is not reported as a difference. Exact numeric tolerance still applies
        // within each matched pair.
        std::vector<int64_t> order_a(static_cast<size_t>(av_list->length()));
        std::vector<int64_t> order_b(order_a.size());
        for (size_t i = 0; i < order_a.size(); ++i)
        {
          order_a[i] = static_cast<int64_t>(i);
          order_b[i] = static_cast<int64_t>(i);
        }
        if (settings.unordered_lists)
        {
          auto key_of = [](const std::shared_ptr<arrow::Array>& arr, int64_t i) {
            auto s = arr->GetScalar(i);
            return s.ok() ? canonicalKey_(*s) : std::string("E;");
          };
          std::stable_sort(order_a.begin(), order_a.end(),
                           [&](int64_t l, int64_t r) { return key_of(av_list, l) < key_of(av_list, r); });
          std::stable_sort(order_b.begin(), order_b.end(),
                           [&](int64_t l, int64_t r) { return key_of(bv_list, l) < key_of(bv_list, r); });
        }

        for (size_t i = 0; i < order_a.size(); ++i)
        {
          auto ea = av_list->GetScalar(order_a[i]);
          auto eb = bv_list->GetScalar(order_b[i]);
          if (!ea.ok() || !eb.ok())
          {
            message = path + "[" + std::to_string(i) + "]: unreadable list element";
            return false;
          }
          if (!scalarsEqual_(*ea, *eb, settings, path + "[" + std::to_string(i) + "]",
                             message, max_ratio, max_absdiff))
          {
            return false;
          }
        }
        return true;
      }

      if (id == arrow::Type::STRUCT)
      {
        const auto& fields_a = static_cast<const arrow::StructScalar&>(*a).value;
        const auto& fields_b = static_cast<const arrow::StructScalar&>(*b).value;
        if (fields_a.size() != fields_b.size())
        {
          message = path + ": struct field count " + std::to_string(fields_a.size()) + " vs. "
                    + std::to_string(fields_b.size());
          return false;
        }
        const auto& struct_type = static_cast<const arrow::StructType&>(*a->type);
        for (size_t i = 0; i < fields_a.size(); ++i)
        {
          const std::string child = path + "." + struct_type.field(static_cast<int>(i))->name();
          if (!scalarsEqual_(fields_a[i], fields_b[i], settings, child, message,
                             max_ratio, max_absdiff))
          {
            return false;
          }
        }
        return true;
      }

      if (!a->Equals(*b))
      {
        message = path + ": " + a->ToString() + " vs. " + b->ToString();
        return false;
      }
      return true;
    }

    /// Append @p msg to @p sink unless the cap is reached, in which case count it as suppressed.
    void report_(std::vector<std::string>& sink, Size max_reported, Size& suppressed,
                 const std::string& msg)
    {
      if (max_reported != 0 && sink.size() >= max_reported) { ++suppressed; return; }
      sink.push_back(msg);
    }

    /// Compare two Arrow schemas by column name, type and nullability.
    void compareSchemas_(const std::shared_ptr<arrow::Table>& t1,
                         const std::shared_ptr<arrow::Table>& t2,
                         const ParquetDiffSettings& settings,
                         ParquetDiffResult& result)
    {
      const auto& s1 = t1->schema();
      const auto& s2 = t2->schema();

      for (const auto& field : s1->fields())
      {
        const int idx = s2->GetFieldIndex(field->name());
        if (idx < 0)
        {
          report_(result.schema_errors, settings.max_reported, result.suppressed,
                  "column '" + field->name() + "' only in file 1");
          continue;
        }
        const auto& other = s2->field(idx);
        if (!field->type()->Equals(*other->type()))
        {
          report_(result.schema_errors, settings.max_reported, result.suppressed,
                  "column '" + field->name() + "': type " + field->type()->ToString()
                    + " vs. " + other->type()->ToString());
        }
        if (field->nullable() != other->nullable())
        {
          report_(result.schema_errors, settings.max_reported, result.suppressed,
                  "column '" + field->name() + "': nullable " + (field->nullable() ? "true" : "false")
                    + " vs. " + (other->nullable() ? "true" : "false"));
        }
      }
      for (const auto& field : s2->fields())
      {
        if (s1->GetFieldIndex(field->name()) < 0)
        {
          report_(result.schema_errors, settings.max_reported, result.suppressed,
                  "column '" + field->name() + "' only in file 2");
        }
      }
    }

    /**
      @brief Build primary keys for every row.

      Returns key -> all rows carrying it, ordered, so that a duplicated key is visible as a
      multiplicity rather than silently collapsing to its first occurrence. A std::map keeps
      iteration deterministic, which matters because ParquetDiffSettings::max_reported decides
      which differences are shown.
    */
    std::map<std::string, std::vector<int64_t>> buildKeys_(const std::shared_ptr<arrow::Table>& table,
                                                           const std::vector<std::string>& key_columns,
                                                           const std::string& which,
                                                           const ParquetDiffSettings& settings,
                                                           ParquetDiffResult& result,
                                                           bool& ok)
    {
      ok = true;
      std::map<std::string, std::vector<int64_t>> keys;

      std::vector<int> indices;
      indices.reserve(key_columns.size());
      for (const auto& name : key_columns)
      {
        const int idx = table->schema()->GetFieldIndex(name);
        if (idx < 0)
        {
          report_(result.key_errors, settings.max_reported, result.suppressed,
                  which + ": primary-key column '" + name + "' does not exist");
          ok = false;
        }
        indices.push_back(idx);
      }
      if (!ok) { return keys; }

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        // canonicalKey_ is self-delimiting, so plain concatenation is unambiguous.
        std::string key;
        for (const int index : indices)
        {
          bool cell_ok = false;
          auto cell = cellAt_(table, index, row, cell_ok);
          if (!cell_ok)
          {
            report_(result.key_errors, settings.max_reported, result.suppressed,
                    which + ": row " + std::to_string(row) + " has an unreadable key cell in '"
                      + table->schema()->field(index)->name() + "'");
            ok = false;
            return keys;
          }
          key += canonicalKey_(cell);
        }
        keys[key].push_back(row);
      }

      for (const auto& [key, rows] : keys)
      {
        if (rows.size() > 1)
        {
          std::string row_list;
          for (size_t i = 0; i < rows.size(); ++i)
          {
            if (i != 0) { row_list += ", "; }
            row_list += std::to_string(rows[i]);
          }
          report_(result.key_errors, settings.max_reported, result.suppressed,
                  which + ": primary key occurs " + std::to_string(rows.size())
                    + " times (rows " + row_list + "): " + key);
        }
      }
      return keys;
    }

  } // anonymous namespace

  std::string ParquetTableComparator::viewFromFileType(const std::string& file_type)
  {
    if (file_type == "psm_file" || file_type == "psm") { return "psm"; }
    if (file_type == "feature_file" || file_type == "feature") { return "feature"; }
    if (file_type == "pg_file" || file_type == "pg") { return "pg"; }
    return "";
  }

  std::vector<std::string> ParquetTableComparator::qpxPrimaryKey(const std::string& view)
  {
    // The keys OpenMS currently emits. QPX >= the identity-id change additionally defines a
    // single surrogate key per view; pass -pk explicitly when comparing such files.
    if (view == "psm")
    {
      return {QPXPSMSchema::PEPTIDOFORM, QPXPSMSchema::CHARGE,
              QPXPSMSchema::RUN_FILE_NAME, QPXPSMSchema::SCAN};
    }
    if (view == "feature")
    {
      return {QPXFeatureSchema::PEPTIDOFORM, QPXFeatureSchema::CHARGE,
              QPXFeatureSchema::RUN_FILE_NAME, QPXFeatureSchema::RT};
    }
    if (view == "pg")
    {
      return {QPXPgSchema::ANCHOR_PROTEIN, QPXPgSchema::GROUPED_RUNS, QPXPgSchema::LABEL};
    }
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "unknown QPX view '" + view + "' (expected psm, feature or pg)");
  }

  ParquetDiffResult ParquetTableComparator::compare(const std::string& file_1,
                                                    const std::string& file_2,
                                                    const ParquetDiffSettings& settings)
  {
    ParquetDiffResult result;

    auto t1 = readTable_(file_1);
    auto t2 = readTable_(file_2);
    if (!t1 || !t2)
    {
      result.schema_errors.push_back("one or both files could not be read as Parquet");
      return result;
    }
    result.rows_1 = static_cast<Size>(t1->num_rows());
    result.rows_2 = static_cast<Size>(t2->num_rows());

    compareSchemas_(t1, t2, settings, result);

    if (settings.schema_only)
    {
      result.equal = result.schema_errors.empty();
      return result;
    }

    // Resolve the primary key: explicit setting first, then the file's QPX file_type.
    std::vector<std::string> key_columns = settings.primary_key;
    if (key_columns.empty())
    {
      const std::string view = viewFromFileType(metadataValue_(t1, "file_type"));
      if (view.empty())
      {
        result.key_errors.push_back(
          "no primary key given and file 1 carries no recognised QPX 'file_type' metadata; "
          "pass the key columns explicitly");
        return result;
      }
      key_columns = qpxPrimaryKey(view);
    }
    result.primary_key_used = key_columns;

    // A type mismatch on a key column makes every key differ, which would bury the real cause.
    if (!result.schema_errors.empty())
    {
      for (const auto& key : key_columns)
      {
        const int i1 = t1->schema()->GetFieldIndex(key);
        const int i2 = t2->schema()->GetFieldIndex(key);
        if (i1 < 0 || i2 < 0 || !t1->schema()->field(i1)->type()->Equals(*t2->schema()->field(i2)->type()))
        {
          result.key_errors.push_back(
            "primary-key column '" + key + "' differs in schema; skipping value comparison");
          return result;
        }
      }
    }

    bool ok_1 = false;
    bool ok_2 = false;
    const auto keys_1 = buildKeys_(t1, key_columns, "file 1", settings, result, ok_1);
    const auto keys_2 = buildKeys_(t2, key_columns, "file 2", settings, result, ok_2);
    if (!ok_1 || !ok_2) { return result; }

    // Key-set drift, including multiplicity: a key present twice on one side and once on the
    // other is a real difference, not something to collapse.
    for (const auto& [key, rows] : keys_1)
    {
      const auto it = keys_2.find(key);
      if (it == keys_2.end())
      {
        report_(result.key_errors, settings.max_reported, result.suppressed,
                "row only in file 1: " + key);
      }
      else if (it->second.size() != rows.size())
      {
        report_(result.key_errors, settings.max_reported, result.suppressed,
                "key occurs " + std::to_string(rows.size()) + " times in file 1 but "
                  + std::to_string(it->second.size()) + " times in file 2: " + key);
      }
    }
    for (const auto& [key, rows] : keys_2)
    {
      (void)rows;
      if (keys_1.find(key) == keys_1.end())
      {
        report_(result.key_errors, settings.max_reported, result.suppressed,
                "row only in file 2: " + key);
      }
    }

    // Columns to compare: those both tables agree on, minus the ignored ones, and minus the
    // primary-key columns. The key columns already matched the rows, and they matched under
    // canonical (order-insensitive) semantics -- re-comparing them as ordinary values would
    // contradict that by reporting a reordered grouped_runs as a difference.
    std::vector<std::pair<int, int>> columns;
    std::vector<std::string> column_names;
    for (const auto& field : t1->schema()->fields())
    {
      if (std::find(settings.ignore_columns.begin(), settings.ignore_columns.end(), field->name())
          != settings.ignore_columns.end())
      {
        continue;
      }
      if (std::find(key_columns.begin(), key_columns.end(), field->name()) != key_columns.end())
      {
        continue;
      }
      const int idx2 = t2->schema()->GetFieldIndex(field->name());
      if (idx2 < 0) { continue; }
      if (!field->type()->Equals(*t2->schema()->field(idx2)->type())) { continue; }
      columns.emplace_back(t1->schema()->GetFieldIndex(field->name()), idx2);
      column_names.push_back(field->name());
    }

    for (const auto& [key, rows_1] : keys_1)
    {
      const auto it = keys_2.find(key);
      if (it == keys_2.end()) { continue; }
      // With a well-formed key both sides hold exactly one row; duplicates were reported above
      // and the first occurrence is compared so that the value report is still useful.
      const int64_t row_1 = rows_1.front();
      const int64_t row_2 = it->second.front();
      ++result.rows_compared;

      for (size_t c = 0; c < columns.size(); ++c)
      {
        bool ok_a = false;
        bool ok_b = false;
        auto cell_a = cellAt_(t1, columns[c].first, row_1, ok_a);
        auto cell_b = cellAt_(t2, columns[c].second, row_2, ok_b);
        if (!ok_a || !ok_b)
        {
          report_(result.value_errors, settings.max_reported, result.suppressed,
                  key + " " + column_names[c] + ": cell could not be read");
          continue;
        }

        std::string message;
        if (!scalarsEqual_(cell_a, cell_b, settings, column_names[c], message,
                           result.max_ratio, result.max_absdiff))
        {
          report_(result.value_errors, settings.max_reported, result.suppressed,
                  key + " " + message);
        }
      }
    }

    result.equal = result.schema_errors.empty() && result.key_errors.empty()
                   && result.value_errors.empty();
    return result;
  }
  ParquetDiffResult ParquetTableComparator::validate(const std::string& file,
                                                     const std::string& view,
                                                     const ParquetDiffSettings& settings)
  {
    ParquetDiffResult result;

    std::shared_ptr<arrow::Schema> expected;
    if (view == "psm") { expected = QPXPSMSchema::schema(); }
    else if (view == "feature") { expected = QPXFeatureSchema::schema(); }
    else if (view == "pg") { expected = QPXPgSchema::schema(); }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "unknown QPX view '" + view + "' (expected psm, feature or pg)");
    }

    auto table = readTable_(file);
    if (!table)
    {
      result.schema_errors.push_back("file could not be read as Parquet");
      return result;
    }
    result.rows_1 = static_cast<Size>(table->num_rows());

    const auto& actual = table->schema();
    for (const auto& field : expected->fields())
    {
      const int idx = actual->GetFieldIndex(field->name());
      if (idx < 0)
      {
        report_(result.schema_errors, settings.max_reported, result.suppressed,
                "missing required column '" + field->name() + "' ("
                  + field->type()->ToString() + ")");
        continue;
      }
      const auto& present = actual->field(idx);
      if (!field->type()->Equals(*present->type()))
      {
        report_(result.schema_errors, settings.max_reported, result.suppressed,
                "column '" + field->name() + "': expected " + field->type()->ToString()
                  + ", found " + present->type()->ToString());
      }
      if (field->nullable() != present->nullable())
      {
        report_(result.schema_errors, settings.max_reported, result.suppressed,
                "column '" + field->name() + "': expected nullable "
                  + (field->nullable() ? "true" : "false") + ", found "
                  + (present->nullable() ? "true" : "false"));
      }
    }
    for (const auto& field : actual->fields())
    {
      if (expected->GetFieldIndex(field->name()) < 0)
      {
        report_(result.schema_errors, settings.max_reported, result.suppressed,
                "column '" + field->name() + "' is not part of the " + view + " view");
      }
    }

    // A QPX primary key must be unique; buildKeys_ reports any duplicate it finds.
    const auto key_columns = settings.primary_key.empty() ? qpxPrimaryKey(view) : settings.primary_key;
    result.primary_key_used = key_columns;
    bool ok = false;
    (void)buildKeys_(table, key_columns, File::basename(file), settings, result, ok);

    result.equal = result.schema_errors.empty() && result.key_errors.empty();
    return result;
  }

} // namespace OpenMS
