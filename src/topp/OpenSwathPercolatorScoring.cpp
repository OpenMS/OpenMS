// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_OpenSwathPercolatorScoring OpenSwathPercolatorScoring

@brief Rescore OpenSWATH OSW / OSWPQ peak groups and transitions with the in-process OpenMS Percolator implementation.

OpenSwathPercolatorScoring operates directly on existing OpenSWATH result containers:
- SQLite OSW files (`.osw`)
- OpenSWATH Parquet archives or directories (`.oswpq`)

One or more scoring levels can be requested sequentially via @p level:
- @c ms1 rescoring writes @c SCORE_MS1 or @c score_ms1_*
- @c ms2 rescoring writes @c SCORE_MS2 or @c score_ms2_*
- @c ms1ms2 combines MS2 and MS1 features and writes the historical MS2 score table/columns
- @c transition rescoring writes @c SCORE_TRANSITION or @c score_transition_*

If @p out is empty, the input is updated in place. Otherwise the input is first copied to @p out and all requested rescoring passes are applied there.

Transition rescoring requires prior precursor-level scores in the same file. Starting from an unscored result, run either `-level ms2 transition` or `-level ms1ms2 transition` in one invocation.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathPercolatorScoring.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathPercolatorScoring.html
*/

/// @cond TOPPCLASSES
class TOPPOpenSwathPercolatorScoring final :
  public TOPPBase
{
public:
  TOPPOpenSwathPercolatorScoring() :
    TOPPBase("OpenSwathPercolatorScoring", "Rescore OpenSWATH OSW / OSWPQ results with in-process Percolator.")
  {
  }

protected:
  using Scorer = OpenMS::OpenSwathPercolatorScoring;
  using Level = Scorer::Level;

  static StringList validLevels_()
  {
    StringList levels;
    levels.reserve(Scorer::names_of_level.size());
    for (const auto& name : Scorer::names_of_level)
    {
      levels.push_back(name);
    }
    return levels;
  }

  static Level toLevel_(const String& value)
  {
    const auto it = std::find(Scorer::names_of_level.begin(),
                              Scorer::names_of_level.end(),
                              static_cast<std::string>(value));
    if (it == Scorer::names_of_level.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unsupported OpenSWATH Percolator level.",
                                    value);
    }
    return static_cast<Level>(std::distance(Scorer::names_of_level.begin(), it));
  }

  static vector<Level> parseLevels_(const StringList& values)
  {
    vector<Level> levels;
    levels.reserve(values.size());
    for (const auto& value : values)
    {
      const Level level = toLevel_(value);
      if (std::find(levels.begin(), levels.end(), level) == levels.end())
      {
        levels.push_back(level);
      }
    }
    if (levels.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "At least one scoring level must be requested.", "");
    }
    return levels;
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input OpenSWATH result file (.osw) or OSWPQ parquet archive/directory (.oswpq).");
    setValidFormats_("in", {"osw", "oswpq"});
    registerOutputFile_("out", "<file>", "", "Optional output path. If empty, the input is updated in place.", false);
    setValidFormats_("out", {"osw", "oswpq"});

    registerStringList_("level", "<levels>", StringList{String("ms2")},
                        "One or more scoring levels to run sequentially. Use 'ms2 transition' or 'ms1ms2 transition' to derive transition scores from an unscored input.",
                        false);
    setValidStrings_("level", validLevels_());

    registerFullParam_(scorer_.getDefaults());
  }

  ExitCodes main_(int, const char**) override
  {
    const String input_file = getStringOption_("in");
    const String output_file = getStringOption_("out");
    const vector<Level> levels = parseLevels_(getStringList_("level"));

    scorer_.setParameters(getParam_().copySubset(scorer_.getDefaults()));

    const String working_file = output_file.empty() ? input_file : output_file;
    bool first_pass = true;
    for (const Level level : levels)
    {
      const auto summary = scorer_.score(first_pass ? input_file : working_file,
                                         level,
                                         first_pass ? output_file : "");
      OPENMS_LOG_INFO << "Level '" << Scorer::names_of_level[static_cast<Size>(level)]
                      << "': rescored " << summary.total_rows << " rows ("
                      << summary.target_rows << " target, "
                      << summary.decoy_rows << " decoy) across "
                      << summary.feature_count << " features.\n";
      first_pass = false;
    }

    OPENMS_LOG_INFO << "OpenSWATH Percolator rescoring finished successfully.\n";
    return EXECUTION_OK;
  }

private:
  Scorer scorer_;
};

int main(int argc, const char** argv)
{
  TOPPOpenSwathPercolatorScoring tool;
  return tool.main(argc, argv);
}

/// @endcond
