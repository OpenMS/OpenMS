// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

#include "FeatureLinkerBase.cpp"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_FeatureLinkerUnlabeledQT FeatureLinkerUnlabeledQT

@brief Groups corresponding features from multiple maps using a QT clustering approach.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> &rarr; FeatureLinkerUnlabeledQT &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature detection algorithm) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MapAlignerPoseClustering @n (or another map alignment algorithm) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
        </tr>
    </table>
</CENTER>

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

This tool provides an algorithm for grouping corresponding features in
multiple runs of label-free experiments. For more details and
algorithm-specific parameters (set in the ini file) see "Detailed
Description" in the @ref OpenMS::FeatureGroupingAlgorithmQT "algorithm
documentation".

FeatureLinkerUnlabeledQT takes several feature maps (featureXML files) and
stores the corresponding features in a consensus map (consensusXML file).
Feature maps can be created from MS experiments (peak data) using one of
the FeatureFinder TOPP tools.

@see @ref TOPP_FeatureLinkerUnlabeled @ref TOPP_FeatureLinkerLabeled

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_FeatureLinkerUnlabeledQT.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_FeatureLinkerUnlabeledQT.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerUnlabeledQT :
  public TOPPFeatureLinkerBase
{

public:
  TOPPFeatureLinkerUnlabeledQT() :
    TOPPFeatureLinkerBase("FeatureLinkerUnlabeledQT", "Groups corresponding features from multiple maps.")
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
    FeatureGroupingAlgorithmQT algo;
    Param p = algo.getParameters();
    return p;
  }

  ExitCodes main_(int, const char **) override
  {
    FeatureGroupingAlgorithmQT algo;
    return TOPPFeatureLinkerBase::common_main_(&algo);
  }

};


int main(int argc, const char ** argv)
{
  TOPPFeatureLinkerUnlabeledQT tool;
  return tool.main(argc, argv);
}

/// @endcond
