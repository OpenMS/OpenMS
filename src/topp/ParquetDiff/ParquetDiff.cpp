// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include "ParquetTableComparator.h"

#include <iostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ParquetDiff ParquetDiff

@brief Compares two Parquet tables by primary key, tolerating numeric differences.

The table counterpart of @ref TOPP_FuzzyDiff. Two Parquet files holding the same logical result
may legitimately differ in row order and in the low-order bits of floating-point columns, so
neither a byte comparison nor a text diff answers "is this the same result?". ParquetDiff matches
rows by a primary key, compares matched rows cell by cell with a numeric tolerance, and reports
schema drift separately from value drift.

Only one of 'ratio' or 'absdiff' has to be satisfied. Use "absdiff" to deal with cases like
"zero vs. epsilon".

Rows are matched on 'pk'. For a QPX file (psm, feature or pg) the key is derived from the file's
own @c file_type metadata when 'pk' is not given. List-valued key columns are canonicalised by
sorting their elements, so a set-valued key such as @c grouped_runs matches regardless of the
order the producer emitted.

With 'schema' the tool takes a single input and checks it against the built-in QPX schema for
that view instead of comparing two files; this reports missing or extra columns, wrong Arrow
types, wrong nullability, and duplicate primary keys.

'min_rows' requires every input to hold at least that many rows, and is checked in both modes.
It closes a gap neither of them covers: two empty tables compare equal, and an empty table
satisfies any schema, so a producer that wrote a correctly shaped nothing passes both. Combining
@p schema with @p min_rows therefore asserts something about a single file without needing a
reference at all - useful where a reference would only add churn. '0' is accepted and always
passes; use it to record deliberately that a table is expected to be empty.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ParquetDiff.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ParquetDiff.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPParquetDiff :
  public TOPPBase
{
public:
  TOPPParquetDiff() :
    TOPPBase("ParquetDiff", "Compares two Parquet tables by primary key, tolerating numeric differences.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    addEmptyLine_();
    registerInputFile_("in1", "<file>", "", "first input file", true, false);
    setValidFormats_("in1", ListUtils::create<std::string>("parquet"));
    registerInputFile_("in2", "<file>", "", "second input file (omit when 'schema' is given)", false, false);
    setValidFormats_("in2", ListUtils::create<std::string>("parquet"));
    addEmptyLine_();

    registerStringList_("pk", "<string list>", ListUtils::create<std::string>(""),
                        "primary-key columns used to match rows. If empty, derived from the QPX "
                        "'file_type' metadata of the first file.", false, false);
    registerStringOption_("schema", "<view>", "",
                          "instead of comparing two files, check 'in1' against the built-in QPX "
                          "schema of this view", false, false);
    setValidStrings_("schema", ListUtils::create<std::string>("psm,feature,pg"));
    registerOutputFile_("out_tsv", "<file>", "",
                        "instead of comparing, write 'in1' out as TSV sorted by primary key. A "
                        "Parquet reference can only be regenerated wholesale, which hides "
                        "unrelated drift; a text dump can be reviewed in a diff and patched line "
                        "by line, and compared with FuzzyDiff like every other reference.",
                        false, false);
    setValidFormats_("out_tsv", ListUtils::create<std::string>("tsv"));
    addEmptyLine_();

    registerDoubleOption_("ratio", "<double>", 1, R"(acceptable relative error. Only one of 'ratio' or 'absdiff' has to be satisfied.  Use "absdiff" to deal with cases like "zero vs. epsilon".)", false, false);
    setMinFloat_("ratio", 1);
    registerDoubleOption_("absdiff", "<double>", 0, "acceptable absolute difference. Only one of 'ratio' or 'absdiff' has to be satisfied. ", false, false);
    setMinFloat_("absdiff", 0);
    addEmptyLine_();

    registerIntOption_("min_rows", "<int>", -1,
                       "require every input to hold at least this many rows ('-1' does not check). "
                       "A schema check alone is satisfied by an empty table, so this is what "
                       "distinguishes 'the tool wrote a correct table' from 'the tool wrote a "
                       "correctly shaped nothing'. '0' is accepted and always passes: use it to "
                       "record that a table is expected to be empty.", false, false);
    setMinInt_("min_rows", -1);
    addEmptyLine_();

    registerStringList_("ignore", "<string list>", ListUtils::create<std::string>(""),
                        "columns excluded from value comparison. They are still schema-checked.", false, true);
    registerFlag_("schema_only", "compare schemas only; do not compare values", true);
    registerFlag_("with_ids",
                  "also compare the QPX identity columns (feature_id/psm_id/pg_id and the "
                  "cross-references between them). Off by default: those values are derived, and "
                  "the spec states that identity is meaningful within a file only, so comparing "
                  "them across two files is a join it tells you not to make. 'schema' mode checks "
                  "instead that they are present, non-null and unique.", true);
    registerFlag_("unordered_lists", "compare list-valued cells as multisets rather than sequences", true);
    registerIntOption_("max_reported", "<int>", 25,
                       "stop listing differences of one kind after this many (0 = unlimited)", false, true);
    setMinInt_("max_reported", 0);

    registerIntOption_("verbose", "<int>", 2, "set verbose level:\n"
                                              "0 = very quiet mode (absolutely no output)\n"
                                              "1 = quiet mode (no output unless differences detected)\n"
                                              "2 = default (include summary at end)\n",
                       false, false);
    setMinInt_("verbose", 0);
    setMaxInt_("verbose", 2);
  }

  /// Print one category of messages, honouring the verbose level.
  static void printSection_(int verbose, const std::string& title,
                            const std::vector<std::string>& messages)
  {
    if (verbose < 1 || messages.empty()) { return; }
    std::cout << title << " (" << messages.size() << "):" << std::endl;
    for (const auto& m : messages) { std::cout << "  " << m << std::endl; }
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in1 = getStringOption_("in1");
    const std::string in2 = getStringOption_("in2");
    const std::string view = getStringOption_("schema");
    const std::string out_tsv = getStringOption_("out_tsv");
    const int verbose = getIntOption_("verbose");

    if (view.empty() && in2.empty() && out_tsv.empty())
    {
      writeLogError_("Error: 'in2' is required unless 'schema' or 'out_tsv' is given.");
      return ILLEGAL_PARAMETERS;
    }

    ParquetDiffSettings settings;
    settings.acceptable_ratio = getDoubleOption_("ratio");
    settings.acceptable_absdiff = getDoubleOption_("absdiff");
    settings.schema_only = getFlag_("schema_only");
    settings.unordered_lists = getFlag_("unordered_lists");
    settings.compare_identity_columns = getFlag_("with_ids");
    settings.max_reported = static_cast<Size>(getIntOption_("max_reported"));

    // Drop empty entries a user may have supplied explicitly; the default is already empty.
    for (const auto& c : getStringList_("pk"))
    {
      if (!c.empty()) { settings.primary_key.push_back(c); }
    }
    for (const auto& c : getStringList_("ignore"))
    {
      if (!c.empty()) { settings.ignore_columns.push_back(c); }
    }

    // Dump mode is a distinct job from comparing, so it returns on its own rather than falling
    // through to a comparison that has no second file to make.
    if (!out_tsv.empty())
    {
      if (!ParquetTableComparator::dumpToTsv(in1, out_tsv, settings)) { return INPUT_FILE_CORRUPT; }
      if (verbose >= 2) { std::cout << "Wrote " << out_tsv << std::endl; }
      return EXECUTION_OK;
    }

    const ParquetDiffResult result = view.empty()
      ? ParquetTableComparator::compare(in1, in2, settings)
      : ParquetTableComparator::validate(in1, view, settings);

    // Row-count floor. Deliberately independent of 'equal': two empty tables compare equal, and an
    // empty table satisfies any schema, so neither mode notices that the producer wrote nothing.
    const int min_rows = getIntOption_("min_rows");
    bool too_few_rows = false;
    if (min_rows >= 0)
    {
      const auto report = [&](const std::string& file, Size rows)
      {
        if (static_cast<Int64>(rows) >= min_rows) { return; }
        too_few_rows = true;
        if (verbose >= 1)
        {
          std::cout << "Row count: '" << file << "' holds " << rows << " row(s), at least "
                    << min_rows << " required" << std::endl;
        }
      };
      report(in1, result.rows_1);
      if (view.empty()) { report(in2, result.rows_2); }
    }

    if (!result.equal || verbose >= 2)
    {
      printSection_(verbose, "Schema differences", result.schema_errors);
      printSection_(verbose, "Key differences", result.key_errors);
      printSection_(verbose, "Value differences", result.value_errors);
    }

    if (verbose >= 2)
    {
      std::cout << "Summary:" << std::endl;
      if (view.empty())
      {
        std::cout << "  rows: " << result.rows_1 << " vs. " << result.rows_2
                  << ", compared: " << result.rows_compared << std::endl;
        std::cout << "  max ratio observed:   " << result.max_ratio << std::endl;
        std::cout << "  max absdiff observed: " << result.max_absdiff << std::endl;
      }
      else
      {
        std::cout << "  view: " << view << ", rows: " << result.rows_1 << std::endl;
      }
      if (!result.primary_key_used.empty())
      {
        std::cout << "  primary key: "
                  << ListUtils::concatenate(result.primary_key_used, ", ") << std::endl;
      }
      if (result.suppressed != 0)
      {
        std::cout << "  " << result.suppressed
                  << " further difference(s) not shown (see 'max_reported')" << std::endl;
      }
      std::cout << ((result.equal && !too_few_rows) ? "  RESULT: equal" : "  RESULT: differ")
                << std::endl;
    }

    return (result.equal && !too_few_rows) ? EXECUTION_OK : INCOMPATIBLE_INPUT_DATA;
  }

};

int main(int argc, const char** argv)
{
  TOPPParquetDiff tool;
  return tool.main(argc, argv);
}

/// @endcond
