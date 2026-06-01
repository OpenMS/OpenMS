// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h>

#include "FeatureLinkerBase.cpp"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureLinkerWNet FeatureLinkerWNet

@brief Groups corresponding features from multiple maps using Wasserstein optimal transport.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; FeatureLinkerWNet &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature detection algorithm) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapAlignerPoseClustering @n (or another map alignment algorithm) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
    </table>
</CENTER>

This tool provides an algorithm for grouping corresponding features in
multiple runs of label-free experiments using Wasserstein (optimal transport)
distance. Features are matched by solving a minimum-cost network flow problem
over their (m/z, RT) positions.

For more details and algorithm-specific parameters (set in the ini file)
see "Detailed Description" in the @ref OpenMS::FeatureGroupingAlgorithmWNet
"algorithm documentation".

FeatureLinkerWNet takes several feature maps (featureXML files) and stores
the corresponding features in a consensus map (consensusXML file). Feature
maps can be created from MS experiments (peak data) using one of the
FeatureFinder TOPP tools.

@see @ref TOPP_FeatureLinkerUnlabeledQT @ref TOPP_FeatureLinkerUnlabeledKD

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureLinkerWNet.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureLinkerWNet.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerWNet :
  public TOPPFeatureLinkerBase
{
public:
  TOPPFeatureLinkerWNet() :
    TOPPFeatureLinkerBase("FeatureLinkerWNet", "Groups corresponding features from multiple maps using Wasserstein optimal transport.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    FeatureGroupingAlgorithmWNet algo;
    Param p = algo.getParameters();
    return p;
  }

  ExitCodes main_(int, const char **) override
  {
    FeatureGroupingAlgorithmWNet algo;
    return TOPPFeatureLinkerBase::common_main_(&algo);
  }
};


int main(int argc, const char ** argv)
{
  TOPPFeatureLinkerWNet tool;
  return tool.main(argc, argv);
}

/// @endcond
