// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <arrow/api.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <set>

namespace OpenMS
{
  namespace
  {
    void throwIfArrowAppendFailed_(const arrow::Status& status, const String& context)
    {
      if (!status.ok())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Failed to build parquet matrix column for " + context + ": " + String(status.ToString()));
      }
    }

    struct PeptideSummaryRow
    {
      String sequence;
      String full_peptide_name;
      String protein_name;
      String gene_name;
      std::map<Int64, double> sample_values;
    };

    double median_(std::vector<double> values)
    {
      if (values.empty())
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
      std::sort(values.begin(), values.end());
      const Size mid = values.size() / 2;
      if (values.size() % 2 == 1)
      {
        return values[mid];
      }
      return (values[mid - 1] + values[mid]) / 2.0;
    }

    std::map<Int64, String> buildSampleNames_(const std::vector<OpenSwathExportRow>& rows)
    {
      std::map<Int64, String> names;
      std::set<String> used;
      for (const auto& row : rows)
      {
        if (names.find(row.run_id) != names.end())
        {
          continue;
        }
        String label = row.run_name.empty() ? (String("RUN_ID_") + String(row.run_id)) : row.run_name;
        if (used.count(label))
        {
          label += "_RUN" + String(row.run_id);
        }
        used.insert(label);
        names[row.run_id] = std::move(label);
      }
      return names;
    }

    std::vector<OpenSwathExportRow> selectBestPeakGroups_(std::vector<OpenSwathExportRow> rows)
    {
      std::sort(rows.begin(), rows.end(),
                [](const auto& lhs, const auto& rhs)
                {
                  if (lhs.run_id != rhs.run_id) return lhs.run_id < rhs.run_id;
                  if (lhs.transition_group_id != rhs.transition_group_id) return lhs.transition_group_id < rhs.transition_group_id;
                  if (lhs.m_score != rhs.m_score) return lhs.m_score < rhs.m_score;
                  if (lhs.intensity != rhs.intensity) return lhs.intensity > rhs.intensity;
                  return lhs.feature_id < rhs.feature_id;
                });

      std::vector<OpenSwathExportRow> best_rows;
      String current_group;
      Int64 current_run = std::numeric_limits<Int64>::min();
      for (const auto& row : rows)
      {
        if (row.run_id != current_run || row.transition_group_id != current_group)
        {
          best_rows.push_back(row);
          current_run = row.run_id;
          current_group = row.transition_group_id;
        }
      }
      return best_rows;
    }

    std::vector<PeptideSummaryRow> summarizePeptides_(const std::vector<OpenSwathExportRow>& best_rows,
                                                      const Size top_n,
                                                      const bool consistent_top)
    {
      using PeptideKey = std::pair<String, String>;
      std::map<PeptideKey, String> protein_by_key;
      std::map<PeptideKey, String> gene_by_key;
      for (const auto& row : best_rows)
      {
        const PeptideKey key{row.sequence, row.full_peptide_name};
        protein_by_key.emplace(key, row.protein_name);
        gene_by_key.emplace(key, row.gene_name);
      }

      std::map<PeptideKey, std::map<Int64, double>> per_run_mean;
      if (consistent_top)
      {
        std::map<PeptideKey, std::map<String, std::vector<double>>> intensities_by_precursor;
        for (const auto& row : best_rows)
        {
          intensities_by_precursor[{row.sequence, row.full_peptide_name}][row.transition_group_id].push_back(row.intensity);
        }

        std::map<PeptideKey, std::set<String>> selected_precursors;
        for (const auto& entry : intensities_by_precursor)
        {
          std::vector<std::pair<String, double>> medians;
          for (const auto& precursor : entry.second)
          {
            medians.emplace_back(precursor.first, median_(precursor.second));
          }
          std::sort(medians.begin(), medians.end(),
                    [](const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
          for (Size i = 0; i < std::min<Size>(top_n, medians.size()); ++i)
          {
            selected_precursors[entry.first].insert(medians[i].first);
          }
        }

        std::map<std::tuple<PeptideKey, Int64>, std::vector<double>> grouped;
        for (const auto& row : best_rows)
        {
          const PeptideKey key{row.sequence, row.full_peptide_name};
          const auto selected_it = selected_precursors.find(key);
          if (selected_it == selected_precursors.end() || !selected_it->second.count(row.transition_group_id))
          {
            continue;
          }
          grouped[{key, row.run_id}].push_back(row.intensity);
        }
        for (const auto& entry : grouped)
        {
          per_run_mean[std::get<0>(entry.first)][std::get<1>(entry.first)] =
            std::accumulate(entry.second.begin(), entry.second.end(), 0.0) / static_cast<double>(entry.second.size());
        }
      }
      else
      {
        std::map<std::tuple<PeptideKey, Int64>, std::vector<OpenSwathExportRow>> grouped_rows;
        for (const auto& row : best_rows)
        {
          grouped_rows[{PeptideKey{row.sequence, row.full_peptide_name}, row.run_id}].push_back(row);
        }

        for (auto& entry : grouped_rows)
        {
          auto& group_rows = entry.second;
          std::sort(group_rows.begin(), group_rows.end(),
                    [](const auto& lhs, const auto& rhs) { return lhs.intensity > rhs.intensity; });
          const Size keep = std::min<Size>(top_n, group_rows.size());
          double sum = 0.0;
          for (Size i = 0; i < keep; ++i)
          {
            sum += group_rows[i].intensity;
          }
          per_run_mean[std::get<0>(entry.first)][std::get<1>(entry.first)] = sum / static_cast<double>(keep);
        }
      }

      std::vector<PeptideSummaryRow> rows;
      rows.reserve(per_run_mean.size());
      for (const auto& entry : per_run_mean)
      {
        PeptideSummaryRow row;
        row.sequence = entry.first.first;
        row.full_peptide_name = entry.first.second;
        row.protein_name = protein_by_key[entry.first];
        row.gene_name = gene_by_key[entry.first];
        row.sample_values = entry.second;
        rows.push_back(std::move(row));
      }
      return rows;
    }

    void validateMatrixDimensions_(const OpenSwathQuantMatrix& matrix)
    {
      if (matrix.identifier_rows.size() != matrix.values.size())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Matrix identifiers and values must have identical row counts.");
      }

      for (Size row_idx = 0; row_idx < matrix.identifier_rows.size(); ++row_idx)
      {
        if (matrix.identifier_rows[row_idx].size() != matrix.identifier_column_names.size())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Identifier row width does not match identifier columns.");
        }
        if (matrix.values[row_idx].size() != matrix.sample_column_names.size())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Value row width does not match sample columns.");
        }
      }
    }

    void applyNormalization_(OpenSwathQuantMatrix& matrix, const OpenSwathMatrixNormalization normalization)
    {
      if (normalization == OpenSwathMatrixNormalization::None)
      {
        return;
      }

      std::vector<double> sample_medians(matrix.sample_column_names.size(), std::numeric_limits<double>::quiet_NaN());
      for (Size col = 0; col < matrix.sample_column_names.size(); ++col)
      {
        std::vector<double> values;
        values.reserve(matrix.values.size());
        for (const auto& row : matrix.values)
        {
          if (row[col].has_value())
          {
            values.push_back(*row[col]);
          }
        }
        sample_medians[col] = median_(values);
      }

      double global_median = std::numeric_limits<double>::quiet_NaN();
      if (normalization == OpenSwathMatrixNormalization::MedianMedian)
      {
        std::vector<double> medians;
        for (const double med : sample_medians)
        {
          if (std::isfinite(med))
          {
            medians.push_back(med);
          }
        }
        global_median = median_(medians);
      }

      for (auto& row : matrix.values)
      {
        for (Size col = 0; col < row.size(); ++col)
        {
          if (!row[col].has_value() || !std::isfinite(sample_medians[col]) || sample_medians[col] == 0.0)
          {
            continue;
          }
          double value = *row[col] / sample_medians[col];
          if (normalization == OpenSwathMatrixNormalization::MedianMedian && std::isfinite(global_median))
          {
            value *= global_median;
          }
          row[col] = value;
        }
      }
    }
  } // namespace

  OpenSwathQuantMatrix OpenSwathMatrixExporter::buildMatrix(const std::vector<OpenSwathExportRow>& rows,
                                                            const OpenSwathMatrixExportConfig& config)
  {
    if (rows.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No OpenSWATH rows were available for matrix export.");
    }
    if (config.level != OpenSwathMatrixLevel::Precursor && config.top_n == 0)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Matrix export requires top_n > 0 for peptide/protein/gene levels.");
    }

    const auto sample_names = buildSampleNames_(rows);
    const auto best_rows = selectBestPeakGroups_(rows);
    OpenSwathQuantMatrix matrix;

    auto emitMatrix = [&](const auto& grouped_rows)
    {
      matrix.identifier_rows.clear();
      matrix.values.clear();
      matrix.sample_column_names.clear();
      for (const auto& name : sample_names)
      {
        matrix.sample_column_names.push_back(name.second);
      }

      for (const auto& grouped_row : grouped_rows)
      {
        matrix.identifier_rows.push_back(grouped_row.first);
        std::vector<std::optional<double>> values(sample_names.size(), std::nullopt);
        Size sample_index = 0;
        for (const auto& sample : sample_names)
        {
          const auto value_it = grouped_row.second.find(sample.first);
          if (value_it != grouped_row.second.end())
          {
            values[sample_index] = value_it->second;
          }
          ++sample_index;
        }
        matrix.values.push_back(std::move(values));
      }
      applyNormalization_(matrix, config.normalization);
    };

    if (config.level == OpenSwathMatrixLevel::Precursor)
    {
      matrix.identifier_column_names = {"transition_group_id", "Sequence", "FullPeptideName", "Charge", "ProteinName"};
      std::map<std::vector<String>, std::map<Int64, double>> grouped;
      for (const auto& row : best_rows)
      {
        grouped[{
          row.transition_group_id,
          row.sequence,
          row.full_peptide_name,
          String(row.charge),
          row.protein_name
        }][row.run_id] = row.intensity;
      }
      emitMatrix(grouped);
      return matrix;
    }

    const auto peptide_rows = summarizePeptides_(best_rows, config.top_n, config.consistent_top);
    if (peptide_rows.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No OpenSWATH rows passed peptide-level matrix summarization.");
    }

    if (config.level == OpenSwathMatrixLevel::Peptide)
    {
      matrix.identifier_column_names = {"Sequence", "FullPeptideName"};
      std::map<std::vector<String>, std::map<Int64, double>> grouped;
      for (const auto& row : peptide_rows)
      {
        grouped[{
          row.sequence,
          row.full_peptide_name
        }] = row.sample_values;
      }
      emitMatrix(grouped);
      return matrix;
    }

    auto summarizeHigher = [&](const bool gene_level)
    {
      matrix.identifier_column_names = {gene_level ? "GeneName" : "ProteinName"};
      std::map<String, std::vector<const PeptideSummaryRow*>> grouped_peptides;
      for (const auto& row : peptide_rows)
      {
        const String& key = gene_level ? row.gene_name : row.protein_name;
        if (key.empty())
        {
          continue;
        }
        grouped_peptides[key].push_back(&row);
      }
      if (grouped_peptides.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      gene_level ? "No gene annotations were available for matrix export."
                                                 : "No protein annotations were available for matrix export.");
      }

      std::map<std::vector<String>, std::map<Int64, double>> grouped;
      for (const auto& entry : grouped_peptides)
      {
        std::vector<const PeptideSummaryRow*> selected = entry.second;
        if (config.consistent_top)
        {
          std::vector<std::pair<const PeptideSummaryRow*, double>> medians;
          for (const auto* peptide_row : entry.second)
          {
            std::vector<double> values;
            for (const auto& sample_value : peptide_row->sample_values)
            {
              values.push_back(sample_value.second);
            }
            medians.emplace_back(peptide_row, median_(values));
          }
          std::sort(medians.begin(), medians.end(),
                    [](const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
          selected.clear();
          for (Size i = 0; i < std::min<Size>(config.top_n, medians.size()); ++i)
          {
            selected.push_back(medians[i].first);
          }
        }

        std::map<Int64, std::vector<double>> per_sample;
        for (const auto* peptide_row : selected)
        {
          for (const auto& sample_value : peptide_row->sample_values)
          {
            per_sample[sample_value.first].push_back(sample_value.second);
          }
        }

        std::map<Int64, double> sample_means;
        for (const auto& sample : per_sample)
        {
          sample_means[sample.first] =
            std::accumulate(sample.second.begin(), sample.second.end(), 0.0) / static_cast<double>(sample.second.size());
        }
        grouped[{entry.first}] = std::move(sample_means);
      }
      emitMatrix(grouped);
    };

    if (config.level == OpenSwathMatrixLevel::Protein)
    {
      summarizeHigher(false);
    }
    else
    {
      summarizeHigher(true);
    }
    return matrix;
  }

  void OpenSwathMatrixExporter::writeMatrix(const String& filename,
                                            const OpenSwathQuantMatrix& matrix,
                                            const OpenSwathMatrixExportConfig& config)
  {
    if (matrix.identifier_rows.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Cannot write an empty OpenSWATH quantification matrix.");
    }
    validateMatrixDimensions_(matrix);

    if (config.format == OpenSwathExportFileFormat::TSV)
    {
      std::ofstream os(filename.c_str());
      if (!os)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      for (Size i = 0; i < matrix.identifier_column_names.size(); ++i)
      {
        if (i != 0) os << '\t';
        os << matrix.identifier_column_names[i];
      }
      for (const auto& sample : matrix.sample_column_names)
      {
        os << '\t' << sample;
      }
      os << '\n';

      for (Size row_idx = 0; row_idx < matrix.identifier_rows.size(); ++row_idx)
      {
        const auto& ids = matrix.identifier_rows[row_idx];
        for (Size i = 0; i < ids.size(); ++i)
        {
          if (i != 0) os << '\t';
          os << ids[i];
        }
        for (const auto& value : matrix.values[row_idx])
        {
          os << '\t';
          if (value.has_value())
          {
            os << *value;
          }
        }
        os << '\n';
      }
      return;
    }

    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;

    for (Size col = 0; col < matrix.identifier_column_names.size(); ++col)
    {
      arrow::StringBuilder builder;
      for (const auto& row : matrix.identifier_rows)
      {
        throwIfArrowAppendFailed_(builder.Append(row[col].c_str()), matrix.identifier_column_names[col]);
      }
      fields.push_back(arrow::field(matrix.identifier_column_names[col].c_str(), arrow::utf8()));
      arrays.push_back(builder.Finish().ValueOrDie());
    }

    for (Size col = 0; col < matrix.sample_column_names.size(); ++col)
    {
      arrow::DoubleBuilder builder;
      for (const auto& row : matrix.values)
      {
        if (row[col].has_value())
        {
          throwIfArrowAppendFailed_(builder.Append(*row[col]), matrix.sample_column_names[col]);
        }
        else
        {
          throwIfArrowAppendFailed_(builder.AppendNull(), matrix.sample_column_names[col]);
        }
      }
      fields.push_back(arrow::field(matrix.sample_column_names[col].c_str(), arrow::float64()));
      arrays.push_back(builder.Finish().ValueOrDie());
    }

    const auto table = arrow::Table::Make(arrow::schema(fields), arrays);
    if (!ArrowIOHelpers::writeTableToParquet(table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }
} // namespace OpenMS
