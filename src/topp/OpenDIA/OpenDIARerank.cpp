// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIARerank.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <string>

using namespace OpenMS;

namespace
{
  class TOPPOpenDIARerank final : public TOPPBase
  {
  public:
    TOPPOpenDIARerank() : TOPPBase("OpenDIARerank", "Train and apply an experimental target-decoy peak-group model", {}, false)
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFile_("train", "<file>", "", "Training-run OpenDIA candidate feature table", true);
      setValidFormats_("train", {"tsv"});
      registerInputFile_("score", "<file>", "", "Held-out-run OpenDIA candidate feature table", true);
      setValidFormats_("score", {"tsv"});
      registerOutputFile_("out", "<file>", "", "Held-out candidate scores and assay-level q-values", true);
      setValidFormats_("out", {"tsv"});
      registerIntOption_("seed", "<int>", 1, "Percolator random seed", false);
      registerIntOption_("threads", "<int>", 1, "Percolator training threads", false);
      setMinInt_("threads", 1);
    }

    ExitCodes main_(int, const char**) override
    {
      std::ifstream training_file(getStringOption_("train"));
      std::ifstream scoring_file(getStringOption_("score"));
      std::ofstream output_file(getStringOption_("out"));
      if (!training_file || !scoring_file || !output_file)
      {
        writeLogError_("Could not open one of the candidate feature input/output files.");
        return CANNOT_WRITE_OUTPUT_FILE;
      }

      try
      {
        const auto training = OpenDIA::readCandidateFeatures(training_file);
        const auto scoring = OpenDIA::readCandidateFeatures(scoring_file);
        const auto result = OpenDIA::rescoreCandidates(training, scoring,
                                                       getIntOption_("seed"), getIntOption_("threads"));
        OpenDIA::writeScoredCandidates(output_file, result);
        output_file.flush();
        if (!output_file)
        {
          writeLogError_("Failed while writing rescored candidate table.");
          return CANNOT_WRITE_OUTPUT_FILE;
        }
        writeLogInfo_("Rescored " + std::to_string(result.size()) + " held-out candidate windows.");
      }
      catch (const std::exception& error)
      {
        writeLogError_(error.what());
        return ILLEGAL_PARAMETERS;
      }
      return EXECUTION_OK;
    }
  };
}

int main(int argc, const char** argv)
{
  TOPPOpenDIARerank tool;
  return tool.main(argc, argv);
}
