// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmRTLocal.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ConsensusMapNormalizer ConsensusMapNormalizer

@brief Normalization of intensities in a set of maps using robust regression.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; ConsensusMapNormalizer &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
        </tr>
    </table>
</CENTER>

The tool normalizes the intensities of a set of maps (consensusXML file). The following normalization algorithms are available:

- Robust regression: Maps are normalized pair-wise relative to the map with the most features. Given two maps, peptide featues are classified as non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) or outliers. From the non-outliers an average intensity ratio is calculated and used for normalization.

- Median correction: The median of all maps is set to the median of the map with the most features.

- Quantile normalization: Performs an exact quantile normalization if the number of features is equal across all maps. Otherwise, an approximate quantile normalization using resampling is applied.

- RT-local median correction ('rt_local'): Every map gets a log2 offset that may vary with retention time, so that both
  different loading and a drift of the ion signal along the gradient are removed. The offset of a feature is the median,
  over the @p rt_local:neighbors features of the same map before it and after it in retention time, of their log2 deviation
  from their mean over the maps. Features quantified in at least @p rt_local:min_maps maps (and passing the accession and
  description filters) serve as references; the offsets are refined in @p rt_local:rounds backfitting rounds.
  See @ref OpenMS::ConsensusMapNormalizerAlgorithmRTLocal for details.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ConsensusMapNormalizer.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ConsensusMapNormalizer.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPConsensusMapNormalizer :
  public TOPPBase
{

public:
  TOPPConsensusMapNormalizer() :
    TOPPBase("ConsensusMapNormalizer", "Normalizes maps of one consensusXML file")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<std::string>("consensusXML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<std::string>("consensusXML"));
    addEmptyLine_();
    registerStringOption_("algorithm_type", "<type>", "robust_regression", "The normalization algorithm that is applied. 'robust_regression' scales each map by a fator computed from the ratios of non-differential background features (as determined by the ratio_threshold parameter), 'quantile' performs quantile normalization, 'median' scales all maps to the same median intensity, 'median_shift' shifts the median instead of scaling (WARNING: if you have regular, log-normal MS data, 'median_shift' is probably the wrong choice. Use only if you know what you're doing!), 'rt_local' scales each map by a factor that may vary with retention time, estimated from the features quantified in several maps (see 'rt_local:*')", false, false);
    setValidStrings_("algorithm_type", ListUtils::create<std::string>("robust_regression,median,median_shift,quantile,rt_local"));
    registerDoubleOption_("ratio_threshold", "<ratio>", 0.67, "Only for 'robust_regression': the parameter is used to distinguish between non-outliers (ratio_threshold < intensity ratio < 1/ratio_threshold) and outliers.", false);
    setMinFloat_("ratio_threshold", 0.001);
    setMaxFloat_("ratio_threshold", 1.0);
    registerStringOption_("accession_filter", "<regexp>", "", "Use only features with accessions (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);
    registerStringOption_("description_filter", "<regexp>", "", "Use only features with description (partially) matching this regular expression for computing the normalization factors. Useful, e.g., if you have known house keeping proteins in your samples. When this parameter is empty or the regular expression matches the empty string, all features are used (even those without an ID). No effect if quantile normalization is used.", false, true);

    registerTOPPSubsection_("rt_local", "Options for the 'rt_local' normalization");
    const ConsensusMapNormalizerAlgorithmRTLocal::Settings rt_local_defaults;
    registerIntOption_("rt_local:neighbors", "<number>", int(rt_local_defaults.neighbors), "Number of reference features on each side (before and after in retention time) whose median log2 deviation is a feature's offset. Fewer: follows faster drifts, but noisier offsets.", false, true);
    setMinInt_("rt_local:neighbors", 1);
    registerIntOption_("rt_local:min_maps", "<number>", int(rt_local_defaults.min_maps), "A consensus feature is a reference feature if it has an intensity in at least this many maps.", false, true);
    setMinInt_("rt_local:min_maps", 2);
    registerIntOption_("rt_local:min_reference", "<number>", int(rt_local_defaults.min_reference), "Maps with fewer reference features are not normalized.", false, true);
    setMinInt_("rt_local:min_reference", 1);
    registerIntOption_("rt_local:rounds", "<number>", int(rt_local_defaults.rounds), "Backfitting rounds: each round recomputes the features' mean levels from the offset-corrected intensities and then the offsets.", false, true);
    setMinInt_("rt_local:rounds", 1);
  }

  ExitCodes main_(int, const char **) override
  {
    std::string in = getStringOption_("in");
    std::string out = getStringOption_("out");
    std::string algo_type = getStringOption_("algorithm_type");
    std::string acc_filter = getStringOption_("accession_filter");
    std::string desc_filter = getStringOption_("description_filter");
    double ratio_threshold = getDoubleOption_("ratio_threshold");

    FileHandler infile;
    ConsensusMap map;
    infile.loadConsensusFeatures(in, map, {FileTypes::CONSENSUSXML}, log_type_);

    //map normalization
    if (algo_type == "robust_regression")
    {
      map.sortBySize();
      vector<double> results = ConsensusMapNormalizerAlgorithmThreshold::computeCorrelation(map, ratio_threshold, acc_filter, desc_filter);
      ConsensusMapNormalizerAlgorithmThreshold::normalizeMaps(map, results);
    }
    else if (algo_type == "median")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SCALE, acc_filter, desc_filter);
    }
    else if (algo_type == "median_shift")
    {
      ConsensusMapNormalizerAlgorithmMedian::normalizeMaps(map, ConsensusMapNormalizerAlgorithmMedian::NM_SHIFT, acc_filter, desc_filter);
    }
    else if (algo_type == "quantile")
    {
      if (!acc_filter.empty() || !desc_filter.empty())
      {
        OPENMS_LOG_WARN << endl << "NOTE: Accession / description filtering is not supported in quantile normalization mode. Ignoring filters." << endl << endl;
      }
      ConsensusMapNormalizerAlgorithmQuantile::normalizeMaps(map);
    }
    else if (algo_type == "rt_local")
    {
      ConsensusMapNormalizerAlgorithmRTLocal::Settings settings;
      settings.neighbors = getIntOption_("rt_local:neighbors");
      settings.min_maps = getIntOption_("rt_local:min_maps");
      settings.min_reference = getIntOption_("rt_local:min_reference");
      settings.rounds = getIntOption_("rt_local:rounds");
      ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(map, settings, acc_filter, desc_filter);
    }
    else
    {
      cerr << "Unknown algorithm type  '" << algo_type.c_str() << "'." << endl;
      return ILLEGAL_PARAMETERS;
    }

    //annotate output with data processing info and save output file
    addDataProcessing_(map, getProcessingInfo_(DataProcessing::NORMALIZATION));
    infile.storeConsensusFeatures(out, map, {FileTypes::CONSENSUSXML});

    return EXECUTION_OK;
  }
};


int main(int argc, const char ** argv)
{
  TOPPConsensusMapNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
