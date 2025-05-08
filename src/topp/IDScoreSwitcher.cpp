// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDScoreSwitcher IDScoreSwitcher

@brief Switches between different scores of peptide hits (PSMs) or protein hits in identification data.

In the idXML file format and in OpenMS' internal representation of identification data, every peptide spectrum match (PSM, "peptide hit") and every protein hit can be associated with a single numeric (quality) score of an arbitrary type. However, database search engines that generate PSMs or tools for post-processing of identification data may assign multiple scores of different types to each PSM/protein. These scores can be captured as meta data associated with the PSMs/proteins (in idXML: "UserParam" elements), but they are typically not considered by TOPP tools that utilize the scores. This utility allows to switch between "primary" scores and scores stored as meta values.

By default this tool operates on PSM scores; to consider protein scores instead, set the @p proteins flag. The meta value that is supposed to replace the PSM/protein score - given by parameter @p new_score - has to be numeric (type "float") and exist for every peptide or protein hit, respectively. The old score will be stored as a meta value, the name for which is given by the parameter @p old_score. It is an error if a meta value with this name already exists for any hit, unless that meta value already stores the same score.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDScoreSwitcher.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDScoreSwitcher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDScoreSwitcher :
  public TOPPBase
{
public:

  TOPPIDScoreSwitcher() :
    TOPPBase("IDScoreSwitcher", "Switches between different scores of peptide or protein hits in identification data")
  {
  }

protected:

  IDScoreSwitcherAlgorithm switcher_{};

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerFullParam_(switcher_.getParameters());
  }

  ExitCodes main_(int, const char**) override
  {
    switcher_.setParameters(getParam_().copySubset(switcher_.getParameters()));
    String in = getStringOption_("in"), out = getStringOption_("out");
    bool do_proteins_ = getFlag_("proteins"); // from full param of IDScoreSwitcherAlgorithm

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    FileHandler().loadIdentifications(in, proteins, peptides, {FileTypes::IDXML});

    Size counter = 0;
    if (do_proteins_)
    {
      for (auto& pid : proteins)
      {
        switcher_.switchScores<ProteinIdentification>(pid, counter);
      }
    }
    else
    {
      for (auto& pepid : peptides)
      {
        switcher_.switchScores<PeptideIdentification>(pepid, counter);
      }
    }

    FileHandler().storeIdentifications(out, proteins, peptides, {FileTypes::IDXML});

    OPENMS_LOG_INFO << "Successfully switched " << counter << " "
             << (do_proteins_ ? "protein" : "PSM") << " scores." << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDScoreSwitcher tool;
  return tool.main(argc, argv);
}

/// @endcond
