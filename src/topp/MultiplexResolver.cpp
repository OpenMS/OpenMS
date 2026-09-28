// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>

using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MultiplexResolver MultiplexResolver

@brief Completes peptide multiplets and resolves conflicts within them.

<CENTER>
  <table>
    <tr>
      <th ALIGN = "center"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> &rarr; MultiplexResolver &rarr;</td>
      <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDConflictResolver </td>
    </tr>
  </table>
</CENTER>

Tools such as FeatureFinderMultiplex can detect peptide feature multiplets in labeled experimental data. The multiplets can then be annotated with peptide sequences
using the IDMapper tool (*). The MultiplexResolver tool is consolidating these results in two steps.
- Any multiplets with conflicting quantitative and sequence information are filtered out. As example, let us consider a triple SILAC analyis. Let us assume a sequence
"LDNLVAIFDINR(Label:13C(6)15N(4))" with a single Arg10 label is mapped to the light feature in a SILAC triplet. Either peptide feature detection or sequence information
must be incorrect und the triplet is removed.
- In a second step, any incomplete peptide feature groups are completed with dummy features of zero intensity. As example, let us stay with the triple SILAC analysis.
But let us now assume the sequence "LDNLVAIFDINR(Label:13C(6)15N(4))" is mapped to the heavy partner of a peptide feature pair. This is no conflict. Medium and heavy
peptides have been correctly detected. The MultiplexResolver adds a dummy peptide feature of zero intensity at the light position and thereby completes the triplet.

(*) Note that the MultiplexResolver tool takes only a single (the first) peptide sequence annotation into account. By running IDConflictResolver first, it is assured that
each multiplet has only one peptide sequence annotation, the best one. Multiplets without sequence annotation are passed to the optional out_conflicts output.

The same steps are available as a library class (@ref OpenMS::MultiplexResolverAlgorithm), which the @ref TOPP_MS1LabeledWorkflow tool
uses to run the complete MS1-labeled quantification in one go.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MultiplexResolver.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MultiplexResolver.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMultiplexResolver :
  public TOPPBase
{
private:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Peptide multiplets with assigned sequence information");
    setValidFormats_("in", ListUtils::create<std::string>("consensusXML"));
    registerInputFile_("in_blacklist", "<file>", "", "Optional input containing spectral peaks blacklisted during feature detection. Needed for generation of dummy features.", false);
    setValidFormats_("in_blacklist", ListUtils::create<std::string>("mzML"));
    registerOutputFile_("out", "<file>", "", "Complete peptide multiplets.");
    setValidFormats_("out", ListUtils::create<std::string>("consensusXML"));
    registerOutputFile_("out_conflicts", "<file>", "", "Optional output containing peptide multiplets without ID annotation or with conflicting quant/ID information.", false);
    setValidFormats_("out_conflicts", ListUtils::create<std::string>("consensusXML"));

    // sections 'algorithm' and 'labels', as defined by the algorithm
    registerFullParam_(MultiplexResolverAlgorithm().getDefaults());
  }

public:
  TOPPMultiplexResolver() :
    TOPPBase("MultiplexResolver", "Completes peptide multiplets and resolves conflicts within them.")
  {
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string in = getStringOption_("in");
    const std::string in_blacklist = getStringOption_("in_blacklist");
    const std::string out = getStringOption_("out");
    const std::string out_conflicts = getStringOption_("out_conflicts");

    // load consensus map
    ConsensusMap map_in;
    FileHandler().loadConsensusFeatures(in, map_in, {FileTypes::CONSENSUSXML}, log_type_);

    // load (optional) blacklist
    MSExperiment exp_blacklist;
    if (!in_blacklist.empty())
    {
      FileHandler().loadExperiment(in_blacklist, exp_blacklist, {FileTypes::MZML}, log_type_);
    }

    // the algorithm's parameters are the 'algorithm' and 'labels' sections of this tool
    MultiplexResolverAlgorithm resolver;
    Param resolver_param = resolver.getDefaults();
    resolver_param.insert("algorithm:", getParam_().copy("algorithm:", true));
    resolver_param.insert("labels:", getParam_().copy("labels:", true));
    writeDebug_("Parameters passed to MultiplexResolverAlgorithm", resolver_param, 3);
    resolver.setParameters(resolver_param);

    // construct the new consensus maps
    ConsensusMap map_out;
    ConsensusMap map_conflicts;
    resolver.resolve(map_in, map_out, map_conflicts, exp_blacklist);

    // store consensus maps
    FileHandler().storeConsensusFeatures(out, map_out, {FileTypes::CONSENSUSXML}, log_type_);
    if (!out_conflicts.empty())
    {
      FileHandler().storeConsensusFeatures(out_conflicts, map_conflicts, {FileTypes::CONSENSUSXML}, log_type_);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPMultiplexResolver tool;
  return tool.main(argc, argv);
}

///@endcond
