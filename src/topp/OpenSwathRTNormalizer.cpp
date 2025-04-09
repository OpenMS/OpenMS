// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, George Rosenberger $
// $Authors: Hannes Roest, George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>

using namespace OpenMS;

//
//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathRTNormalizer OpenSwathRTNormalizer

@brief The OpenSwathRTNormalizer will find retention time peptides in data.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; OpenSwathRTNormalizer &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2></td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathAnalyzer </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenSwathWorkflow </td>
        </tr>
    </table>
</CENTER>

This tool will find retention time normalization peptides in data and use them
to generate a transformation between the experimental RT space and the
normalized RT space. The output is a transformation file on how to transform
the RT space into the normalized space.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathRTNormalizer.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathRTNormalizer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathRTNormalizer : public TOPPBase
{
public:

  TOPPOpenSwathRTNormalizer() :
  TOPPBase("OpenSwathRTNormalizer",
           "Generate a transformation file on how to transform the RT space into the normalized space given a description of RT peptides and their normalized retention time.",
           true)
  {
  }

protected:

  typedef PeakMap MapType;

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file with the RT peptides ('TraML' or 'csv')");
    setValidFormats_("tr", ListUtils::create<String>("csv,traML"));
    
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("trafoXML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library)", false);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false);

    registerFlag_("estimateBestPeptides", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).", false);

    registerSubsection_("algorithm", "Algorithm parameters section");

    registerSubsection_("peptideEstimation", "Parameters for the peptide estimation (use -estimateBestPeptides to enable).");

    registerSubsection_("RTNormalization", "Parameters for the RTNormalization. RT normalization and outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
  }

  Param getSubsectionDefaults_(const String & section) const override
  {
    if (section == "algorithm")
    {
      return MRMFeatureFinderScoring().getDefaults();
    }
    else if (section == "peptideEstimation")
    {
      Param p;
      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (section == "RTNormalization")
    {
      Param p;
      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", {"iter_residual","iter_jackknife","ransac","none"});

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", {"true","false"});

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      return p;
    }
    return Param();
  }

  ExitCodes main_(int, const char **) override
  {

    ///////////////////////////////////
    // Read input files and parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file_str = getStringOption_("tr");
    String out = getStringOption_("out");
    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");
    bool estimateBestPeptides = getFlag_("estimateBestPeptides");
    const char * tr_file  = tr_file_str.c_str();

    MapType all_xic_maps; // all XICs from all files
    OpenSwath::LightTargetedExperiment targeted_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      TargetedExperiment transition_exp_;
      FileHandler().loadTransitions(tr_file, transition_exp_, {FileTypes::TRAML});
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);
    }

    Param pepEstimationParams = getParam_().copy("peptideEstimation:", true);
    Param RTNormParams = getParam_().copy("RTNormalization:", true);
    String outlier_method = RTNormParams.getValue("outlierMethod").toString();

    // 1. Estimate the retention time range of the whole experiment
    std::pair<double,double> RTRange = OpenSwathHelper::estimateRTRange(targeted_exp);
    std::cout << "Detected retention time range from " << RTRange.first << " to " << RTRange.second << std::endl;

    // 2. Store the peptide retention times in an intermediate map
    std::map<std::string, double> PeptideRTMap;
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      PeptideRTMap[targeted_exp.getCompounds()[i].id] = targeted_exp.getCompounds()[i].rt; 
    }

    TransformationDescription trafo;

    // If we have a transformation file, trafo will transform the RT in the
    // scoring according to the model. If we don't have one, it will apply the
    // null transformation.
    if (!getStringOption_("rt_norm").empty())
    {
      String trafo_in = getStringOption_("rt_norm");
      String model_type = "linear"; //getStringOption_("model:type");
      FileHandler().loadTransformations(trafo_in, trafo, true, {FileTypes::TRANSFORMATIONXML});
    }

    ///////////////////////////////////
    // Start computation
    ///////////////////////////////////

    // 3. Extract the RT pairs from the input data
    std::vector<std::pair<double, double> > pairs;
    for (Size i = 0; i < file_list.size(); ++i)
    {
      boost::shared_ptr<MapType> swath_map (new MapType()); // the map with the extracted ion chromatograms
      boost::shared_ptr<MapType> xic_map (new MapType());
      FeatureMap featureFile;
      std::cout << "RT Normalization working on " << file_list[i] << std::endl;
      FileHandler().loadExperiment(file_list[i], *xic_map.get(), {FileTypes::MZML}, log_type_);

      // Initialize the featureFile and set its parameters (disable for example
      // the RT score since here do not know the RT transformation) 
      MRMFeatureFinderScoring featureFinder;
      Param scoring_params = getParam_().copy("algorithm:", true);
      scoring_params.setValue("Scores:use_rt_score", "false");
      scoring_params.setValue("Scores:use_elution_model_score", "false");
      if (estimateBestPeptides)
      {
        scoring_params.setValue("TransitionGroupPicker:compute_peak_quality", "true");
        scoring_params.setValue("TransitionGroupPicker:minimal_quality", pepEstimationParams.getValue("InitialQualityCutoff"));
      }
      featureFinder.setParameters(scoring_params);
      featureFinder.setStrictFlag(false);
      
      std::vector< OpenSwath::SwathMap > swath_maps(1);
      swath_maps[0].sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(xic_map);
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;
      featureFinder.pickExperiment(chromatogram_ptr, featureFile, targeted_exp, trafo, swath_maps, transition_group_map);

      // add all the chromatograms to the output
      for (Size k = 0; k < xic_map->getChromatograms().size(); k++)
      {
        all_xic_maps.addChromatogram(xic_map->getChromatograms()[k]);
      }

      // find most likely correct feature for each group and add it to the
      // "pairs" vector by computing pairs of iRT and real RT
      std::map<std::string, double> res = OpenSwathHelper::simpleFindBestFeature(transition_group_map, 
        estimateBestPeptides, pepEstimationParams.getValue("OverallQualityCutoff"));
      for (std::map<std::string, double>::iterator it = res.begin(); it != res.end(); ++it)
      {
        pairs.push_back(std::make_pair(it->second, PeptideRTMap[it->first])); // pair<exp_rt, theor_rt>
      }
    }

    // 4. Perform the outlier detection
    std::vector<std::pair<double, double> > pairs_corrected;
    if (outlier_method == "iter_residual" || outlier_method == "iter_jackknife")
    {
      pairs_corrected = MRMRTNormalizer::removeOutliersIterative(pairs, min_rsq, min_coverage,
        RTNormParams.getValue("useIterativeChauvenet").toBool(), outlier_method);
    }
    else if (outlier_method == "ransac")
    {
      // First, estimate of the maximum deviation from RT that is tolerated:
      //   Because 120 min gradient can have around 4 min elution shift, we use
      //   a default value of 3 % of the gradient to find upper RT threshold (3.6 min).
      double pcnt_rt_threshold = RTNormParams.getValue("RANSACMaxPercentRTThreshold");
      double max_rt_threshold = (RTRange.second - RTRange.first) * pcnt_rt_threshold / 100.0;

      pairs_corrected = MRMRTNormalizer::removeOutliersRANSAC(pairs, min_rsq, min_coverage,
        RTNormParams.getValue("RANSACMaxIterations"), max_rt_threshold,
        RTNormParams.getValue("RANSACSamplingSize"));
    }
    else if (outlier_method == "none") 
    {
      pairs_corrected = pairs;
    }
    else 
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("Illegal argument '") + outlier_method + "' used for outlierMethod (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none').");
    }

    // 5. Check whether the found peptides fulfill the binned coverage criteria
    // set by the user.
    bool enoughPeptides = MRMRTNormalizer::computeBinnedCoverage(RTRange, pairs_corrected,
      pepEstimationParams.getValue("NrRTBins"),
      pepEstimationParams.getValue("MinPeptidesPerBin"),
      pepEstimationParams.getValue("MinBinsFilled") );
    if (estimateBestPeptides && !enoughPeptides)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "There were not enough bins with the minimal number of peptides");
    }

    ///////////////////////////////////
    // Write output
    ///////////////////////////////////

    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    String model_type = "linear";
    trafo_out.fitModel(model_type, model_params);
    FileHandler().storeTransformations(out, trafo_out, {FileTypes::TRANSFORMATIONXML});

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathRTNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
