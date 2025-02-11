// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDConflictResolver IDConflictResolver

@brief Resolves ambiguous annotations of features with peptide identifications.

<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> potential predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> &rarr; IDConflictResolver &rarr;</td>
        <th ALIGN = "center"> potential successor tools </td>
    </tr>
    <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping algorithm) </td>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
    </tr>
</table>
</CENTER>

The peptide identifications are filtered so that only one identification
with a single hit (with the best score) is associated to each feature. (If
two IDs have the same best score, either one of them may be selected.)

The the filtered identifications are added to the vector of unassigned peptides
and also reduced to a single best hit.

This step may be useful before applying @ref TOPP_ProteinQuantifier
"ProteinQuantifier", because features with ambiguous annotation are not
considered for the quantification.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDConflictResolver.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDConflictResolver.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDConflictResolver :
  public TOPPBase
{
public:

  TOPPIDConflictResolver() :
    TOPPBase("IDConflictResolver", "Resolves ambiguous annotations of features with peptide identifications")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (data with one peptide identification per feature)");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));
    registerStringOption_("resolve_between_features", "<resolve_between_features>", "off", "A map may contain multiple features with both identical (possibly modified i.e. not stripped) sequence and charge state. The feature with the 'highest intensity' is very likely the most reliable one. When switched on, the filter removes the sequence annotation from the lower intensity features, thereby resolving the multiplicity. Only the most reliable features for each (possibly modified i.e. not stripped) sequence maintain annotated with this peptide sequence.", false);
    setValidStrings_("resolve_between_features", ListUtils::create<String>("off,highest_intensity"));
  }

  ExitCodes main_(int, const char **) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");
    String resolve_between_features = getStringOption_("resolve_between_features");
    
    FileTypes::Type in_type = FileHandler::getType(in);
    
    if (in_type == FileTypes::FEATUREXML) // featureXML
    {
      FeatureMap features;
      FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
      
      IDConflictResolverAlgorithm::resolve(features);
      
      if (resolve_between_features=="highest_intensity")
      {
        IDConflictResolverAlgorithm::resolveBetweenFeatures(features);
      }
      
      addDataProcessing_(features, getProcessingInfo_(DataProcessing::FILTERING));
      FileHandler().storeFeatures(out, features, {FileTypes::FEATUREXML});
    }
    else // consensusXML
    {
      ConsensusMap consensus;
      FileHandler().loadConsensusFeatures(in, consensus, {FileTypes::CONSENSUSXML});
      
      IDConflictResolverAlgorithm::resolve(consensus);
      
      if (resolve_between_features=="highest_intensity")
      {
        IDConflictResolverAlgorithm::resolveBetweenFeatures(consensus);
      }
      
      addDataProcessing_(consensus, getProcessingInfo_(DataProcessing::FILTERING));
      FileHandler().storeConsensusFeatures(out, consensus, {FileTypes::CONSENSUSXML});
    }
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDConflictResolver tool;
  return tool.main(argc, argv);
}

/// @endcond
