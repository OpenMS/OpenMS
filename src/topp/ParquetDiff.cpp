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
#include <OpenMS/FORMAT/ParquetTableComparator.h>

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
    addEmptyLine_();

    registerDoubleOption_("ratio", "<double>", 1, R"(acceptable relative error. Only one of 'ratio' or 'absdiff' has to be satisfied.  Use "absdiff" to deal with cases like "zero vs. epsilon".)", false, false);
    setMinFloat_("ratio", 1);
    registerDoubleOption_("absdiff", "<double>", 0, "acceptable absolute difference. Only one of 'ratio' or 'absdiff' has to be satisfied. ", false, false);
    setMinFloat_("absdiff", 0);
    addEmptyLine_();

    registerStringList_("ignore", "<string list>", ListUtils::create<std::string>(""),
                        "columns excluded from value comparison. They are still schema-checked.", false, true);
    registerFlag_("schema_only", "compare schemas only; do not compare values", true);
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
    const int verbose = getIntOption_("verbose");

    if (view.empty() && in2.empty())
    {
      writeLogError_("Error: 'in2' is required unless 'schema' is given.");
      return ILLEGAL_PARAMETERS;
    }

    ParquetDiffSettings settings;
    settings.acceptable_ratio = getDoubleOption_("ratio");
    settings.acceptable_absdiff = getDoubleOption_("absdiff");
    settings.schema_only = getFlag_("schema_only");
    settings.unordered_lists = getFlag_("unordered_lists");
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

    const ParquetDiffResult result = view.empty()
      ? ParquetTableComparator::compare(in1, in2, settings)
      : ParquetTableComparator::validate(in1, view, settings);

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
      std::cout << (result.equal ? "  RESULT: equal" : "  RESULT: differ") << std::endl;
    }

    return result.equal ? EXECUTION_OK : INCOMPATIBLE_INPUT_DATA;
  }

};

int main(int argc, const char** argv)
{
  TOPPParquetDiff tool;
  return tool.main(argc, argv);
}

/// @endcond
