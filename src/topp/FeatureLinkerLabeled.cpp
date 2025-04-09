// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>

#include "FeatureLinkerBase.cpp"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureLinkerLabeled FeatureLinkerLabeled

@brief Groups corresponding isotope-labeled features in a feature map.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureLinkerLabeled &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureFinderCentroided @n (or another feature detection algorithm) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
    </table>
</CENTER>


This tool provides an algorithm for grouping corresponding features in isotope-labeled experiments. For more details and algorithm-specific parameters (set in the ini file) see "Detailed Description" in the @ref OpenMS::FeatureGroupingAlgorithmLabeled "algorithm documentation".

FeatureLinkerLabeled takes one feature map (featureXML file) and stores the corresponding features in a consensus map (consensusXML file). Feature maps can be created from MS experiments (peak data) using one of the FeatureFinder TOPP tools.

@see @ref TOPP_FeatureLinkerUnlabeled @ref TOPP_FeatureLinkerUnlabeledQT

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureLinkerLabeled.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureLinkerLabeled.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerLabeled :
  public TOPPFeatureLinkerBase
{

public:
  TOPPFeatureLinkerLabeled() :
    TOPPFeatureLinkerBase("FeatureLinkerLabeled", "Groups corresponding isotope-labeled features in a feature map.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file", true);
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file", true);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    FeatureGroupingAlgorithmLabeled algo;
    Param p = algo.getParameters();
    return p;
  }

  ExitCodes main_(int, const char **) override
  {
    FeatureGroupingAlgorithmLabeled algo;
    return TOPPFeatureLinkerBase::common_main_(&algo, true);
  }

};


int main(int argc, const char ** argv)
{
  TOPPFeatureLinkerLabeled tool;
  return tool.main(argc, argv);
}

/// @endcond
